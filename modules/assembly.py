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

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file
from fasta import FASTATarium
from minimap2 import Minimap2, Paf2dotplot, parse_paf
from ragtag import RagTag
from phasing import HaplotypePartitioner, deduplicate_query, build_footprints, unique_coverage, group_coverage

def get_fasta_lengths(fastaFiles):
    '''
    Parameters:
        fastaFiles -- a list of FASTA file locations
    Returns:
        lengthsDict -- a dictionary where keys are sequence IDs and values are integers
                       giving the sequence's length, in file order
    '''
    lengthsDict = {}
    for fastaFile in fastaFiles:
        with read_gz_file(fastaFile) as fileIn:
            for record in FASTATarium.fasta_parser(fileIn):
                if record.id in lengthsDict:
                    raise ValueError(f"'{record.id}' occurs more than once in the FASTA file(s) " +
                                     f"(most recently in '{fastaFile}'); sequence IDs must be unique")
                lengthsDict[record.id] = len(record)
    return lengthsDict

def format_regions(footprint):
    '''
    Parameters:
        footprint -- a dictionary where keys are reference contig IDs and values are lists
                     of merged (start, end) intervals
    Returns:
        regions -- a comma-separated string of 'tid:start-end' spans (1-based, inclusive)
                   ordered by the number of aligned bases, most first
        alignedBp -- a comma-separated string of the aligned bases in the same order
    '''
    alignedBp = { tid: sum([ end - start for start, end in intervals ]) for tid, intervals in footprint.items() }
    order = sorted(footprint, key=lambda tid: (-alignedBp[tid], tid))
    regions = ",".join([ f"{tid}:{footprint[tid][0][0]+1}-{footprint[tid][-1][1]}" for tid in order ])
    return regions, ",".join([ str(alignedBp[tid]) for tid in order ])

# assembly > phase
def assembly_phase(args):
    '''
    Parameters:
        args.inputFiles -- a list of query FASTA file locations
        args.referenceFile -- a string of the reference FASTA file location
        args.ploidy -- an integer of how many haplotype groups to partition into
        args.pafFile -- a string of where the minimap2 PAF file is (or will be) located
        args.outputFileName -- a string of where to write the phasing TSV
        args.coverageFileName -- a string of where to write the coverage TSV
    '''
    # Obtain sequence lengths, also ensuring IDs are unique
    queryLengths = get_fasta_lengths(args.inputFiles)
    refLengths = get_fasta_lengths([args.referenceFile])

    # Align query contigs against the reference
    if os.path.isfile(args.pafFile + ".ok"):
        print(f"# Reusing existing minimap2 alignment '{args.pafFile}'")
    else:
        print("# Running minimap2 alignment...")
        aligner = Minimap2(args.minimap2, preset=args.preset, threads=args.threads)
        aligner.align(args.referenceFile, args.inputFiles, args.pafFile)
        open(args.pafFile + ".ok", "w").close()

    # Summarise the reference regions each contig aligns to
    alignments = parse_paf(args.pafFile, minQ=args.minQ, minAlignLen=args.minAlignLen)
    unknownIDs = set(alignments).difference(queryLengths)
    if unknownIDs != set():
        raise ValueError(f"{len(unknownIDs)} query ID(s) in '{args.pafFile}' are not in the -i FASTA file(s), " +
                         f"such as '{sorted(unknownIDs)[0]}'; was the alignment made from different inputs?")
    footprints = build_footprints(deduplicate_query(alignments))

    # Partition contigs into haplotype groups
    print(f"# Partitioning {len(footprints)} aligned contig(s) into {args.ploidy} group(s)...")
    partitioner = HaplotypePartitioner(args.ploidy, penalty=args.overlapPenalty)
    states = partitioner.partition(footprints)
    uniqueBp = unique_coverage(states, footprints)

    # Write per-contig phasing results
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write("#contig\tlength\tpartition\tstatus\taligned_bp\tunique_bp\treference_regions\treference_bp\n")
        for qid, length in queryLengths.items():
            if not qid in footprints:
                fileOut.write(f"{qid}\t{length}\t.\tunaligned\t0\t0\t.\t.\n")
                continue

            regions, regionsBp = format_regions(footprints[qid])
            alignedBp = sum([ int(x) for x in regionsBp.split(",") ])
            if states[qid] is None:
                fileOut.write(f"{qid}\t{length}\t.\tredundant\t{alignedBp}\t0\t{regions}\t{regionsBp}\n")
            else:
                fileOut.write(f"{qid}\t{length}\t{states[qid]+1}\tassigned\t{alignedBp}\t" +
                              f"{uniqueBp[qid]}\t{regions}\t{regionsBp}\n")

    # Write and display per-reference coverage of each partition
    coverage = group_coverage(states, footprints, args.ploidy)
    partitionNames = [ f"partition{i+1}" for i in range(args.ploidy) ]
    idWidth = max([ len(tid) for tid in refLengths ] + [len("#reference")]) + 2

    print("# Percentage of each reference contig covered by each partition's alignments:")
    print(f"{'#reference':<{idWidth}}" + "".join([ f"{name:>12}" for name in partitionNames ]))
    with open(args.coverageFileName, "w") as fileOut:
        fileOut.write("#reference\tlength\t" + "\t".join([ f"{name}_bp\t{name}_pct" for name in partitionNames ]) + "\n")
        for tid, length in refLengths.items():
            groupBp = coverage.get(tid, [0] * args.ploidy)
            groupPct = [ round(bp / length * 100, 2) for bp in groupBp ]
            fileOut.write(f"{tid}\t{length}\t" + "\t".join([ f"{bp}\t{pct}" for bp, pct in zip(groupBp, groupPct) ]) + "\n")
            print(f"{tid:<{idWidth}}" + "".join([ f"{pct:>12.2f}" for pct in groupPct ]))

    # Report on contig assignments
    numAssigned = [ list(states.values()).count(i) for i in range(args.ploidy) ]
    numRedundant = list(states.values()).count(None)
    numUnaligned = len(queryLengths) - len(footprints)
    print("# Contigs per partition: " + ", ".join([ f"{name}={num}" for name, num in zip(partitionNames, numAssigned) ]))
    print(f"# {numRedundant} redundant contig(s) were left unassigned")
    print(f"# {numUnaligned} contig(s) had no qualifying alignments and were left unassigned")

# assembly > scaffold
def assembly_scaffold(args):
    '''
    Parameters:
        args.partitions -- a dictionary as produced by parsing.parse_phasing_tsv()
        args.partitionPaths -- a dictionary where keys are partition labels and values are
                               dictionaries of file locations used for that partition
        args.fastaFiles -- a list of contig FASTA file locations
        args.referenceFile -- a string of the reference FASTA file location
        args.paf2dotplot -- a string of the paf2dotplot.r file location, or None to skip dotplotting
    '''
    # Make sure all assigned contigs can be found
    contigLengths = get_fasta_lengths(args.fastaFiles)
    missingIDs = [ contigID for contigs in args.partitions.values() for contigID in contigs if not contigID in contigLengths ]
    if missingIDs != []:
        raise ValueError(f"{len(missingIDs)} contig(s) in the phasing TSV were not found in the -f FASTA file(s), " +
                         f"such as '{missingIDs[0]}'; were they the same files given to 'assembly phase'?")

    # Write each partition's contigs to file, skipping those that were already scaffolded
    toWrite = { partition: paths for partition, paths in args.partitionPaths.items() if not os.path.exists(paths["flag"]) }
    contigToPartition = { contigID: partition for partition in toWrite for contigID in args.partitions[partition] }
    if toWrite != {}:
        fileOuts = {}
        try:
            for partition, paths in toWrite.items():
                os.makedirs(paths["workDir"])
                fileOuts[partition] = open(os.path.join(paths["workDir"], "contigs.fasta"), "w")
            for fastaFile in args.fastaFiles:
                with read_gz_file(fastaFile) as fileIn:
                    for record in FASTATarium.fasta_parser(fileIn):
                        if record.id in contigToPartition:
                            fileOuts[contigToPartition[record.id]].write(record.format())
        finally:
            for fileOut in fileOuts.values():
                fileOut.close()

    # Scaffold each partition with RagTag
    scaffolder = RagTag(args.ragtag, minimap2Path=args.minimap2, threads=args.threads)
    for partition, paths in args.partitionPaths.items():
        name = f"partition{partition}"
        ragtagDir = os.path.join(paths["workDir"], "ragtag")
        if os.path.exists(paths["flag"]):
            print(f"# Reusing existing RagTag scaffolding of {name}")
        else:
            print(f"# Running RagTag scaffolding of {name}...")
            scaffolder.scaffold(args.referenceFile, os.path.join(paths["workDir"], "contigs.fasta"), ragtagDir)
            open(paths["flag"], "w").close()

        # Write scaffolds without the RagTag suffix
        RagTag.rename_scaffolds(os.path.join(ragtagDir, "ragtag.scaffold.fasta"), paths["fasta"],
                                lineLength=args.lineLength)

        # Report on how contigs were scaffolded
        scaffolds, placed = RagTag.parse_agp(os.path.join(ragtagDir, "ragtag.scaffold.agp"))
        numPlaced = sum([ len(scaffolds[scaffoldID]) for scaffoldID in placed ])
        unplaced = [ scaffoldID for scaffoldID in scaffolds if not scaffoldID in placed ]
        print(f"# {name}: {numPlaced} contig(s) placed into {len(placed)} scaffold(s); " +
              f"{len(unplaced)} contig(s) left unplaced")
        for scaffoldID in scaffolds:
            if scaffoldID in placed:
                print(f"#   {scaffoldID}: {', '.join(scaffolds[scaffoldID])}")
        if unplaced != []:
            print(f"#   unplaced: {', '.join(unplaced)}")

    # Generate dotplots of each scaffolded partition against the reference
    if args.paf2dotplot != None:
        aligner = Minimap2(args.minimap2, preset=args.preset, threads=args.threads)
        plotter = Paf2dotplot(args.paf2dotplot, plotSize=args.plotSize, minRefLen=args.minRefLen,
                              minQueryLen=args.minQueryLen, minAlignLen=args.minPlotAlignLen)
        for partition, paths in args.partitionPaths.items():
            print(f"# Generating dotplot of partition{partition}...")
            aligner.align(args.referenceFile, [paths["fasta"]], paths["paf"])
            plotter.plot(paths["paf"], paths["plotPrefix"])
