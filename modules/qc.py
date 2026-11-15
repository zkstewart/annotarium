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
from collections import Counter

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import GzCapableWriter, FILE_CONTIG_HEADER
from assembly import get_fasta_lengths
from compleasm import CompleasmTable, SetOptimiser, SetSummary

def check_table_contigs(table, contigLengths, fastaFile, runDirectory):
    '''
    Ensures that a compleasm run was performed against the FASTA file it is paired with.

    Parameters:
        table -- a CompleasmTable object
        contigLengths -- a dictionary of the FASTA file's contig IDs to their lengths
        fastaFile -- a string indicating the location of the paired FASTA file
        runDirectory -- a string indicating the location of the compleasm run directory
    '''
    unknownContigs = table.contigs - set(contigLengths)
    if unknownContigs != set():
        examples = ", ".join(f"'{x}'" for x in sorted(unknownContigs)[:5])
        raise ValueError(f"Contig(s) in the compleasm results of '{runDirectory}' are not in " +
                         f"'{fastaFile}' (e.g., {examples}); are -i and -f given in the same order?")

def ordered_contig_ids(fileContigs):
    '''
    Parameters:
        fileContigs -- a list with one list of contig IDs per FASTA file
    Returns:
        contigIDs -- a list of unique contig IDs in order of first appearance
    '''
    contigIDs = {} # dict acts as an ordered set
    for contigs in fileContigs:
        for contigID in contigs:
            contigIDs[contigID] = None
    return list(contigIDs)

def report_copy_shortfall(tables, genes, numSets):
    '''
    Reports how many gene copies are missing from the genome as a whole, regardless of which
    set or contig they are on; these cannot be recovered by shuffling contigs between sets.
    '''
    totalCopies = Counter()
    for table in tables:
        for contigCounts in table.completeHits.values():
            totalCopies.update(contigCounts)

    shortGenes = [ geneID for geneID in genes if totalCopies[geneID] < numSets ]
    shortCopies = sum(numSets - totalCopies[geneID] for geneID in shortGenes)
    absentGenes = sum(1 for geneID in genes if totalCopies[geneID] == 0)
    expectedCopies = numSets * len(genes)
    print(f"# Complete gene copies short of {numSets} per gene across all contigs: {shortCopies} of " +
          f"{expectedCopies} ({100 * shortCopies / expectedCopies:.2f}%) in {len(shortGenes)} genes; " +
          f"{absentGenes} genes are absent from all contigs")

def qc_compleasm_set_optimise(args):
    '''
    Parameters:
        args.runDirectories -- a list of compleasm output directories, one per haplotype
        args.tableFiles -- a list of 'full_table_busco_format.tsv' files within those directories
        args.fastaFiles -- a list of FASTA files compleasm was run on, in the same order
        args.setFileNames -- a list of output TSV file locations, one per haplotype set
        args.restarts -- an integer giving the number of random restarts for optimisation
    '''
    numSets = len(args.fastaFiles)
    fastaNames = [ os.path.basename(fastaFile) for fastaFile in args.fastaFiles ]

    # Parse compleasm results and check they belong to their paired FASTA file
    tables = [ CompleasmTable(tableFile) for tableFile in args.tableFiles ]
    lineages = { table.lineage for table in tables }
    if len(lineages) > 1:
        raise ValueError(f"compleasm results use different lineages ({', '.join(sorted(lineages))}); " +
                         "all runs must use the same lineage")
    genes = set.union(*[ table.genes for table in tables ])

    fileLengths = [ get_fasta_lengths([fastaFile]) for fastaFile in args.fastaFiles ]
    fileContigs = [ list(contigLengths) for contigLengths in fileLengths ]
    for table, contigLengths, fastaFile, runDirectory in zip(tables, fileLengths, args.fastaFiles, args.runDirectories):
        check_table_contigs(table, contigLengths, fastaFile, runDirectory)

    # Group contigs which share an ID across files, as these are treated as homologous
    contigGenes = {}
    groups = {}
    for fileIndex, (table, contigs) in enumerate(zip(tables, fileContigs)):
        for contigID in contigs:
            groups.setdefault(contigID, []).append(fileIndex)
            if contigID in table.completeHits:
                contigGenes[(fileIndex, contigID)] = table.completeHits[contigID]

    # Optimise set membership
    optimiser = SetOptimiser(contigGenes, groups, numSets, restarts=args.restarts)
    if optimiser.groups == {}:
        raise ValueError("No contig IDs carrying complete genes are shared between the -f FASTA files, " +
                         "so there is nothing to shuffle; homologous contigs must have the same ID in each file")
    print(f"# {len(optimiser.groups)} homologous contig group(s) can be shuffled between sets")

    before = optimiser.identity()
    after = optimiser.optimise()
    beforeSets = optimiser.contig_sets(before, fileContigs)
    afterSets = optimiser.contig_sets(after, fileContigs)
    beforeSummary = SetSummary(beforeSets, tables, genes, numSets)
    afterSummary = SetSummary(afterSets, tables, genes, numSets)

    # Report on results
    print("# Completeness per FASTA file (before optimisation):")
    for setIndex in range(numSets):
        print(f"set{setIndex+1} ({fastaNames[setIndex]}): {beforeSummary.format(setIndex)}")
    print("# Completeness per haplotype set (after optimisation):")
    for setIndex in range(numSets):
        print(f"set{setIndex+1}: {afterSummary.format(setIndex)}")

    print("# Contigs moved from another file's set:")
    for setIndex in range(numSets):
        moved = [ f"{contigID} (from {fastaNames[fileIndex]})"
                  for (fileIndex, contigID), contigSet in afterSets.items()
                  if contigSet == setIndex and fileIndex != setIndex ]
        print(f"set{setIndex+1}: {', '.join(moved) if moved != [] else 'none'}")
    report_copy_shortfall(tables, genes, numSets)

    # Write one file+contig TSV per set, ordered by contig ID's first appearance
    contigIDs = ordered_contig_ids(fileContigs)
    for setIndex, setFileName in enumerate(args.setFileNames):
        with GzCapableWriter(setFileName) as fileOut:
            fileOut.write(FILE_CONTIG_HEADER)
            for contigID in contigIDs:
                for fileIndex in groups[contigID]:
                    if afterSets[(fileIndex, contigID)] == setIndex:
                        fileOut.write(f"{fastaNames[fileIndex]}\t{contigID}\n")
