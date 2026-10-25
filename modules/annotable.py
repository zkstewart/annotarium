# Copyright (C) 2026 Zachary Kenneth Stewart

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os, sys
import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, GzCapableWriter
from gff3tarium import GFF3Feature, GFF3Tarium
from validation import DirectoryNotFoundError

class Annotable:
    '''
    Columns in HIT_COLUMNS store one value per BLAST hit, ordered from best to worst and
    formatted like 'hit1 [hit2][hit3]'. Columns in IDMAPPING_COLUMNS are derived from the
    best hit which had an ID mapping (as given in 'Best_hit_with_idmapping').
    '''
    HIT_COLUMNS = ["Target_accession", "Gene_names", "Taxa_names", "Percentage_identity",
                   "Alignment_length", "Mismatches", "Gap_opens", "Query_start", "Query_end",
                   "Target_start", "Target_end", "Expect_value", "Bit_score"]
    IDMAPPING_COLUMNS = ["Best_hit_with_idmapping", "Best_mapped_GOs", "Best_mapped_GOs_+_parents"]
    
    def __init__(self, annotableFile, indexCol="#Query"):
        self.annotableFile = annotableFile
        self.df = None
        self.load_annotable(self.annotableFile, indexCol=indexCol)
        self.isAnnotable = True # flag for easier type checking
    
    @property
    def annotableFile(self):
        return self._annotableFile
    
    @annotableFile.setter
    def annotableFile(self, value):
        if not isinstance(value, str):
            raise ValueError("File location must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"File location '{value}' is not a file")
        
        self._annotableFile = value
    
    @property
    def nrow(self):
        return len(self.df)
    
    @property
    def ncol(self):
        return len(self.df.columns)
    
    @property
    def rowNames(self):
        return self.df.index.tolist()
    
    @property
    def rows(self):
        for index, rowSeries in self.df.iterrows():
            yield rowSeries
    
    @property
    def colNames(self):
        return self.df.columns.tolist()
    
    @property
    def columns(self):
        for colIndex in self.df.columns:
            yield self.df[colIndex]
    
    def load_annotable(self, annotableFile, indexCol="#Query"):
        with read_gz_file(annotableFile) as fileIn:
            self.df = pd.read_csv(fileIn, sep="\t", index_col=indexCol,
                                  dtype=str, keep_default_na=False) # retain values exactly as written
    
    def write(self, outputFileName):
        '''
        Writes the annotable object stored under self.df to file. File overwriting is not allowed.
        
        Parameters:
            outputFileName -- a string pointing to a file location that does not already exist.
                              If full path is given, the parent folder(s) must already exist.
        '''
        outputFileName = os.path.abspath(outputFileName)
        if os.path.exists(outputFileName):
            raise FileExistsError(f"Cannot .write annotable to '{outputFileName}' as it already exists!")
        parentDir = os.path.dirname(outputFileName)
        if not os.path.isdir(parentDir):
            raise DirectoryNotFoundError(f"Cannot .write annotable '{os.path.basename(outputFileName)}' as " + 
                                         f"its parent directory '{parentDir}' does not exist!")
        
        with GzCapableWriter(outputFileName) as fileOut:
            self.df.to_csv(fileOut, sep="\t")
    
    @staticmethod
    def split_hits(value, numHits):
        '''
        Parameters:
            value -- a string formatted like 'hit1 [hit2][hit3]', or '.' if there are no hits
            numHits -- an integer indicating how many hits the value should contain
        Returns:
            hits -- a list of strings, one per hit; or None if the value could not be split
                    into numHits hits
        '''
        if numHits == 0:
            return [] if value == "." else None
        if numHits == 1:
            return [value]
        if not value.endswith("]"):
            return None
        
        # Split from the right, since only the first hit is not enclosed in brackets
        hits = value[:-1].rsplit("][", maxsplit=numHits-2)
        firstAndSecond = hits[0].rsplit(" [", maxsplit=1)
        hits = firstAndSecond + hits[1:]
        if len(hits) != numHits:
            return None
        return hits
    
    @staticmethod
    def join_hits(hits):
        '''
        Parameters:
            hits -- a list of strings, one per hit
        Returns:
            value -- a string formatted like 'hit1 [hit2][hit3]', or '.' if there are no hits
        '''
        if len(hits) == 0:
            return "."
        if len(hits) == 1:
            return hits[0]
        return hits[0] + " " + "".join([ f"[{hit}]" for hit in hits[1:] ])
    
    def filter_hits(self, evalue=None, maxhits=None):
        '''
        Removes hits worse than the evalue cutoff, then retains up to maxhits of the remaining
        hits. If no hits remain, the row's hit columns are set to '.'; the row is not removed.
        If the 'Best_hit_with_idmapping' hit is removed, the ID mapping columns are set to '.'.
        
        Parameters:
            evalue -- OPTIONAL; a float of the maximum E-value for a hit to be retained
            maxhits -- OPTIONAL; an integer of the maximum number of hits to retain
        Returns:
            numModified -- an integer indicating how many rows had hits removed
        '''
        hitColumns = [ col for col in Annotable.HIT_COLUMNS if col in self.df.columns ]
        idmappingColumns = [ col for col in Annotable.IDMAPPING_COLUMNS if col in self.df.columns ]
        colIndex = { col:self.df.columns.get_loc(col) for col in self.df.columns }
        
        numModified = 0
        for rowIndex in range(len(self.df)):
            rowID = self.df.index[rowIndex]
            
            # Split hit columns into their individual hits
            accessions = self.df.iat[rowIndex, colIndex["Target_accession"]]
            numHits = 0 if accessions == "." else len(accessions.replace(" ", "").replace("]", "").split("["))
            rowHits = {}
            for col in hitColumns:
                hits = Annotable.split_hits(self.df.iat[rowIndex, colIndex[col]], numHits)
                if hits == None:
                    raise ValueError(f"'{rowID}' column '{col}' could not be split into the {numHits} " +
                                     "hits indicated by its 'Target_accession' column")
                rowHits[col] = hits
            
            # Identify which hits to keep
            toKeep = list(range(numHits))
            if evalue != None:
                toKeep = [ i for i in toKeep if float(rowHits["Expect_value"][i]) <= evalue ]
            if maxhits != None:
                toKeep = toKeep[:maxhits]
            if len(toKeep) == numHits:
                continue
            numModified += 1
            
            # Update hit columns
            for col in hitColumns:
                self.df.iat[rowIndex, colIndex[col]] = Annotable.join_hits([ rowHits[col][i] for i in toKeep ])
            if toKeep == [] and "Source" in colIndex:
                self.df.iat[rowIndex, colIndex["Source"]] = "."
            
            # Update ID mapping columns if their hit was removed
            if "Best_hit_with_idmapping" in colIndex:
                keptAccessions = [ rowHits["Target_accession"][i] for i in toKeep ]
                if not self.df.iat[rowIndex, colIndex["Best_hit_with_idmapping"]] in keptAccessions:
                    for col in idmappingColumns:
                        self.df.iat[rowIndex, colIndex[col]] = "."
        
        return numModified
    
    def __delitem__(self, key):
        if not key in self:
            raise KeyError(f"'{key}' cannot be deleted as it is not contained within this MSA")
        
        # Drop the sequence row
        self.df.drop(key, inplace=True)
    
    def __getitem__(self, key):
        return self.df.loc[key]
    
    def __len__(self):
        return len(self.df)
    
    def __iter__(self):
        yield from self.df.iterrows()
    
    def __contains__(self, value):
        return value in self.df.index
    
    def __repr__(self):
        return "<Annotable object;file='{0}';num_genes={1};num_annotation_cols={2}>".format(
            self.annotableFile,
            self.nrow,
            self.ncol
        )

def annotable_ids(args):
    annotable = Annotable(args.annotableFile)
    with GzCapableWriter(args.outputFileName) as fileOut:
        for rowID in annotable.rowNames:
            fileOut.write(f"{rowID}\n")

def annotable_tx2gene(args):
    # Parse annotable
    annotable = Annotable(args.annotableFile)
    
    # Parse GFF3
    gff3 = GFF3Tarium(args.gff3File)
    
    # Extract tx2gene and gene2tx out of GFF3
    tx2gene = {}
    gene2tx = {}
    for geneID in gff3.ftypes["gene"]:
        geneFeature = gff3[geneID]
        if hasattr(geneFeature, "mRNA"):
            gene2tx.setdefault(geneID, set())
            for mrnaFeature in geneFeature.mRNA:
                tx2gene[mrnaFeature.ID] = geneFeature.ID
                gene2tx[geneID].add(mrnaFeature.ID)
    
    # Identify the best hit for each gene
    bestGeneTx = {}
    for txID, row in annotable.df.iterrows():
        geneID = tx2gene[txID]
        try:
            bitscore = float(row["Bit_score"].split(" [")[0])
        except:
            bitscore = None
        
        bestGeneTx.setdefault(geneID, [txID, None])
        prevBest = bestGeneTx[geneID][1]
        if (bitscore is not None) and ((prevBest is None) or (bitscore > prevBest)):
            bestGeneTx[geneID] = [txID, bitscore]
    
    # Subset annotable DataFrame down to just the best transcript IDs
    toKeep = [ txID for geneID, (txID, bitscore) in bestGeneTx.items() ]
    annotable.df = annotable.df.loc[toKeep]
    
    # Convert txID to geneID
    annotable.df.rename(index=tx2gene, inplace=True)
    
    # Write to file
    annotable.write(args.outputFileName)

def annotable_goseq(args):
    # Parse annotable
    annotable = Annotable(args.annotableFile)
    if not args.columnHeader in annotable.df:
        raise KeyError(f"--col value '{args.columnHeader}' is not a header in -i '{args.annotableFile}' ; " + 
                       f"available columns are {annotable.df.columns.to_list()}")
    
    # Extract the column and reset nulls
    gos = annotable.df[args.columnHeader]
    gos.replace(".", args.nullCharacter, inplace=True)
    
    # Write to file
    with GzCapableWriter(args.outputFileName) as fileOut:
        gos.to_csv(fileOut, sep="\t", header=False)

def annotable_filter(args):
    '''
    Hits are filtered first by --evalue and --maxhits, such that --include and --exclude
    are only matched against the values that remain.
    
    Parameters:
        args.columns -- a list of column headers to search for --include and --exclude strings
        args.include -- a list of strings; if non-empty, only rows containing at least
                        one of these strings are retained
        args.exclude -- a list of strings; rows containing any of these strings are
                        removed (takes precedence over args.include)
    '''
    annotable = Annotable(args.annotableFile)
    
    # Filter hits within each row
    if args.evalue != None or args.maxhits != None:
        numModified = annotable.filter_hits(evalue=args.evalue, maxhits=args.maxhits)
        print(f"# {numModified} row(s) had hits removed by --evalue and/or --maxhits")
    
    # Filter rows by string matching
    if args.include != [] or args.exclude != []:
        searchIndex = annotable.df.index.name in args.columns
        searchColumns = [ col for col in args.columns if col != annotable.df.index.name ]
        
        toKeep = []
        for rowID, row in annotable.df[searchColumns].iterrows():
            searchValues = ([rowID] if searchIndex else []) + row.tolist()
            if args.include != [] and not any([ x in value for x in args.include for value in searchValues ]):
                continue
            if any([ x in value for x in args.exclude for value in searchValues ]):
                continue
            toKeep.append(rowID)
        
        print(f"# {len(toKeep)} of {len(annotable)} row(s) retained by --include and/or --exclude")
        annotable.df = annotable.df.loc[toKeep]
    
    # Write to file
    annotable.write(args.outputFileName)
