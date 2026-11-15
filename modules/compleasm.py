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

import os, random, sys
from collections import Counter
from glob import glob
from itertools import permutations

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file

TABLE_NAME = "full_table_busco_format.tsv"
COMPLETE_STATUSES = ["Complete", "Duplicated"]
FRAGMENT_STATUSES = ["Fragmented", "Interspersed"]
LINEAGE_HEADER = "# The lineage dataset is: "
SUMMARY_NAME = "summary.txt"
SUMMARY_LINEAGE = "## lineage: "
SUMMARY_COLUMNS = ["S", "D", "F", "M", "I", "N"] # 'I' is only reported by genome runs

def locate_compleasm_table(runDirectory):
    '''
    Parameters:
        runDirectory -- a string indicating the location of a compleasm output directory,
                        or the lineage subdirectory within it
    Returns:
        tableFile -- a string indicating the location of the 'full_table_busco_format.tsv' file
    '''
    tableFile = os.path.join(runDirectory, TABLE_NAME)
    if os.path.isfile(tableFile):
        return tableFile

    tableFiles = glob(os.path.join(runDirectory, "*", TABLE_NAME))
    if len(tableFiles) == 0:
        raise FileNotFoundError(f"No '{TABLE_NAME}' file was found in '{runDirectory}' or its " +
                                "lineage subdirectory; is this a compleasm output directory?")
    if len(tableFiles) > 1:
        raise ValueError(f"'{runDirectory}' contains results for more than one lineage; give the " +
                         "location of the lineage subdirectory you want to use instead")
    return tableFiles[0]

def find_compleasm_summaries(directories):
    '''
    Parameters:
        directories -- a list of strings indicating the locations of compleasm output
                       directories, and/or directories containing them
    Returns:
        summaryFiles -- a list of 'summary.txt' file locations, sorted within each directory
    '''
    summaryFiles = []
    for directory in directories:
        summaryFile = os.path.join(directory, SUMMARY_NAME)
        if os.path.isfile(summaryFile):
            summaryFiles.append(summaryFile)
        else:
            summaryFiles += sorted(glob(os.path.join(directory, "*", SUMMARY_NAME)))
    return summaryFiles

def sample_name_from_summary(summaryFile):
    '''
    Names a sample after the directory holding its 'summary.txt' file, minus any
    'compleasm' suffix (e.g., 'hap1_compleasm' becomes 'hap1').
    '''
    sampleName = os.path.basename(os.path.dirname(summaryFile))
    if sampleName.endswith("compleasm") and sampleName != "compleasm":
        sampleName = sampleName[:-len("compleasm")].rstrip("_.")
    return sampleName

def parse_compleasm_summary(summaryFile):
    '''
    Parameters:
        summaryFile -- a string indicating the location of a compleasm 'summary.txt' file
    Returns:
        lineageStats -- a dictionary with structure like:
                        {
                            "eudicotyledons_odb12": { "S": "98.75", "D": "0.89", ..., "N": "2805" },
                            ...
                        }
                        where values are percentages, except for 'N' which is the number of genes
    '''
    lineageStats = {}
    lineage = None
    with read_gz_file(summaryFile) as fileIn:
        for line in fileIn:
            l = line.strip()
            if l.startswith(SUMMARY_LINEAGE):
                lineage = l[len(SUMMARY_LINEAGE):].strip()
                lineageStats[lineage] = {}
                continue
            if l.startswith("#") or l == "":
                continue
            if lineage is None or not ":" in l:
                raise ValueError(f"'{summaryFile}' has an unexpected line ('{l}'); is it a compleasm summary?")
            
            # Lines look like 'S:98.75%, 2770' or 'N:2805'
            key, value = l.split(":", maxsplit=1)
            lineageStats[lineage][key] = value.split("%")[0].strip()
    
    # Make sure each lineage has the expected statistics
    if lineageStats == {}:
        raise ValueError(f"No lineage results were found in '{summaryFile}'; is it a compleasm summary?")
    for lineage, stats in lineageStats.items():
        missingKeys = [ key for key in SUMMARY_COLUMNS if key != "I" and not key in stats ]
        if missingKeys != []:
            raise ValueError(f"'{summaryFile}' lacks {', '.join(missingKeys)} value(s) for lineage '{lineage}'")
    return lineageStats

class CompleasmTable:
    '''
    Holds the gene hits from a compleasm 'full_table_busco_format.tsv' file.

    Attributes:
        genes -- a set of all gene IDs in the lineage, including those which are missing
        completeHits -- a dictionary where keys are contig IDs and values are Counters
                        of gene IDs to their number of complete copies on that contig
        fragmentHits -- a dictionary where keys are contig IDs and values are sets of
                        gene IDs fragmented on that contig
    '''
    def __init__(self, tableFile):
        self.tableFile = tableFile
        self.lineage = os.path.basename(os.path.dirname(tableFile)) # updated from the header if present
        self.genes = set()
        self.completeHits = {}
        self.fragmentHits = {}
        self._parse()

    def _parse(self):
        with read_gz_file(self.tableFile) as fileIn:
            for line in fileIn:
                if line.startswith(LINEAGE_HEADER):
                    self.lineage = line[len(LINEAGE_HEADER):].split()[0]
                if line.startswith("#") or line.strip() == "":
                    continue
                sl = line.rstrip("\r\n").split("\t")
                geneID, status = sl[0], sl[1]
                self.genes.add(geneID)

                if status in COMPLETE_STATUSES:
                    contigID = sl[2]
                    self.completeHits.setdefault(contigID, Counter())[geneID] += 1
                elif status in FRAGMENT_STATUSES:
                    contigID = sl[2]
                    self.fragmentHits.setdefault(contigID, set()).add(geneID)

        if self.genes == set():
            raise ValueError(f"No genes were found in '{self.tableFile}'; is it a compleasm table?")

    @property
    def contigs(self):
        return set(self.completeHits) | set(self.fragmentHits)

class SetSummary:
    '''
    Tallies BUSCO-style completeness statistics for each haplotype set.

    Parameters:
        contigSets -- a dictionary where keys are (fileIndex, contigID) tuples and
                      values are the index of the haplotype set the contig belongs to
        tables -- a list of CompleasmTable objects, one per file index
        genes -- a set of all gene IDs in the lineage
        numSets -- an integer giving the number of haplotype sets
    '''
    def __init__(self, contigSets, tables, genes, numSets):
        self.numGenes = len(genes)
        self.single, self.duplicated, self.fragmented, self.missing = [], [], [], []

        for setIndex in range(numSets):
            copies = Counter()
            fragments = set()
            for (fileIndex, contigID), contigSet in contigSets.items():
                if contigSet != setIndex:
                    continue
                copies.update(tables[fileIndex].completeHits.get(contigID, {}))
                fragments.update(tables[fileIndex].fragmentHits.get(contigID, set()))

            single = sum(1 for count in copies.values() if count == 1)
            duplicated = sum(1 for count in copies.values() if count > 1)
            fragmented = len(fragments - set(copies))
            self.single.append(single)
            self.duplicated.append(duplicated)
            self.fragmented.append(fragmented)
            self.missing.append(self.numGenes - single - duplicated - fragmented)

    def format(self, setIndex):
        complete = self.single[setIndex] + self.duplicated[setIndex]
        pct = lambda x: f"{100 * x / self.numGenes:.2f}%"
        return (f"C:{pct(complete)} [S:{pct(self.single[setIndex])}, D:{pct(self.duplicated[setIndex])}], " +
                f"F:{pct(self.fragmented[setIndex])}, M:{pct(self.missing[setIndex])}, n:{self.numGenes}")

class SetOptimiser:
    '''
    Shuffles homologous contigs between haplotype sets to maximise the number of genes
    which are complete in each set and, secondarily, the number which are single-copy
    (i.e., minimising duplicated genes).

    Contigs sharing an ID across files are treated as homologous, and their copies are
    shuffled among the sets of the files containing them. Contigs whose ID occurs in only
    one file, or which carry no complete genes, stay in their own file's set.

    Parameters:
        contigGenes -- a dictionary where keys are (fileIndex, contigID) tuples and values
                       are Counters of gene IDs to their number of complete copies
        groups -- a dictionary where keys are contig IDs and values are lists of the file
                  indices containing that contig ID, in ascending order
        numSets -- an integer giving the number of haplotype sets (i.e., the number of files)
        restarts -- an integer giving the number of random restarts for hill climbing
        seed -- an integer to seed random restarts with, for reproducibility
    '''
    def __init__(self, contigGenes, groups, numSets, restarts=30, seed=1):
        self.contigGenes = contigGenes
        self.numSets = numSets
        self.restarts = restarts
        self.seed = seed

        # Only groups spanning multiple files and carrying genes are worth shuffling
        self.groups = {
            contigID: fileIndices for contigID, fileIndices in groups.items()
            if len(fileIndices) > 1 and any((fileIndex, contigID) in contigGenes for fileIndex in fileIndices)
        }
        self.options = { contigID: list(permutations(fileIndices)) for contigID, fileIndices in self.groups.items() }

        # Contigs which can't move always contribute their genes to their own file's set
        self.fixedCounts = [ Counter() for _ in range(numSets) ]
        for (fileIndex, contigID), genes in contigGenes.items():
            if not contigID in self.groups:
                self.fixedCounts[fileIndex].update(genes)

    @staticmethod
    def gene_score(copies):
        '''
        Returns:
            score -- a tuple of (isComplete, isSingleCopy) as integers, which sums correctly
                     when compared lexicographically
        '''
        return (int(copies >= 1), int(copies == 1))

    def identity(self):
        '''
        Returns:
            assignment -- a dictionary where keys are contig IDs and values are tuples giving
                          the set index of each file's copy, in the order of self.groups;
                          here, each copy stays in its own file's set
        '''
        return { contigID: tuple(fileIndices) for contigID, fileIndices in self.groups.items() }

    def set_counts(self, assignment):
        '''
        Returns:
            counts -- a list with one Counter per set, of gene IDs to their complete copies
        '''
        counts = [ Counter(fixed) for fixed in self.fixedCounts ]
        for contigID, setIndices in assignment.items():
            for fileIndex, setIndex in zip(self.groups[contigID], setIndices):
                counts[setIndex].update(self.contigGenes.get((fileIndex, contigID), {}))
        return counts

    def score(self, assignment):
        '''
        Returns:
            score -- a tuple of (complete, single) giving the number of complete and
                     single-copy genes summed across all sets
        '''
        complete, single = 0, 0
        for counts in self.set_counts(assignment):
            for copies in counts.values():
                isComplete, isSingle = self.gene_score(copies)
                complete += isComplete
                single += isSingle
        return (complete, single)

    def _move_delta(self, counts, contigID, oldSets, newSets):
        '''
        Computes the change in score from moving a contig group's copies between sets,
        considering only the genes those copies carry.

        Returns:
            delta -- a tuple of (complete, single) score changes
        '''
        # Tally how each set's gene copies change
        changes = {} # (setIndex, geneID): change in copies
        for fileIndex, oldSet, newSet in zip(self.groups[contigID], oldSets, newSets):
            if oldSet == newSet:
                continue
            for geneID, copies in self.contigGenes.get((fileIndex, contigID), {}).items():
                changes[(oldSet, geneID)] = changes.get((oldSet, geneID), 0) - copies
                changes[(newSet, geneID)] = changes.get((newSet, geneID), 0) + copies

        # Sum the per-gene score changes
        deltaComplete, deltaSingle = 0, 0
        for (setIndex, geneID), change in changes.items():
            before = counts[setIndex][geneID]
            oldComplete, oldSingle = self.gene_score(before)
            newComplete, newSingle = self.gene_score(before + change)
            deltaComplete += newComplete - oldComplete
            deltaSingle += newSingle - oldSingle
        return (deltaComplete, deltaSingle)

    def _apply_move(self, counts, contigID, oldSets, newSets):
        for fileIndex, oldSet, newSet in zip(self.groups[contigID], oldSets, newSets):
            genes = self.contigGenes.get((fileIndex, contigID), {})
            counts[oldSet].subtract(genes)
            counts[newSet].update(genes)

    def _hill_climb(self, assignment):
        '''
        Repeatedly applies the best-improving rearrangement of any one contig group until
        no rearrangement improves the score.

        Returns:
            assignment -- the locally optimal assignment
        '''
        assignment = dict(assignment)
        counts = self.set_counts(assignment)
        improved = True
        while improved:
            improved = False
            for contigID in self.groups:
                bestDelta, bestOption = (0, 0), None
                for option in self.options[contigID]:
                    delta = self._move_delta(counts, contigID, assignment[contigID], option)
                    if delta > bestDelta:
                        bestDelta, bestOption = delta, option
                if bestOption is not None:
                    self._apply_move(counts, contigID, assignment[contigID], bestOption)
                    assignment[contigID] = bestOption
                    improved = True
        return assignment

    def _retained(self, assignment):
        'Returns the number of contig copies left in their own file\'s set'
        return sum(fileIndex == setIndex
                   for contigID, setIndices in assignment.items()
                   for fileIndex, setIndex in zip(self.groups[contigID], setIndices))

    def _relabel(self, assignment):
        '''
        Sets are interchangeable when every contig with genes can move, so equally scoring
        optima may differ only by which set is called which. This function renames sets
        to keep as many contigs in their original file's set as possible, without
        changing the score.

        Returns:
            assignment -- the relabelled assignment
        '''
        bestScore = self.score(assignment)
        best = assignment
        for labels in permutations(range(self.numSets)):
            relabelled = { contigID: tuple(labels[s] for s in setIndices) for contigID, setIndices in assignment.items() }

            # Relabelling must not send a copy to a set its group can't occupy
            if any(set(setIndices) != set(self.groups[contigID]) for contigID, setIndices in relabelled.items()):
                continue
            if self.score(relabelled) == bestScore and self._retained(relabelled) > self._retained(best):
                best = relabelled
        return best

    def optimise(self):
        '''
        Returns:
            assignment -- the best assignment found, starting from the original sets
                          then from random restarts; ties keep the earliest found
        '''
        rng = random.Random(self.seed)
        best = self._hill_climb(self.identity())
        bestScore = self.score(best)
        for _ in range(self.restarts):
            start = { contigID: rng.choice(options) for contigID, options in self.options.items() }
            candidate = self._hill_climb(start)
            candidateScore = self.score(candidate)
            if candidateScore > bestScore:
                best, bestScore = candidate, candidateScore
        return self._relabel(best)

    def contig_sets(self, assignment, fileContigs):
        '''
        Parameters:
            assignment -- an assignment as returned by optimise() or identity()
            fileContigs -- a list with one list of contig IDs per file index
        Returns:
            contigSets -- a dictionary where keys are (fileIndex, contigID) tuples and
                          values are the set index the contig belongs to
        '''
        contigSets = {}
        for fileIndex, contigIDs in enumerate(fileContigs):
            for contigID in contigIDs:
                if contigID in assignment:
                    position = self.groups[contigID].index(fileIndex)
                    contigSets[(fileIndex, contigID)] = assignment[contigID][position]
                else:
                    contigSets[(fileIndex, contigID)] = fileIndex
        return contigSets
