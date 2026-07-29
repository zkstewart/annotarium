#! python3

import os, sys
import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, GzCapableWriter
from gff3tarium import GFF3Feature, GFF3Tarium

class Annotable:
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
            self.df = pd.read_csv(fileIn, sep="\t", index_col=indexCol)
    
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
