#! python3

import os, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, GzCapableWriter

def rnammer_reformat(args):
    ongoingCounts = {}
    
    with read_gz_file(args.rnammerGff2) as fileIn, GzCapableWriter(args.outputFileName) as fileOut:
        for line in fileIn:
            # Parse line and skip if irrelevant
            l = line.strip()
            if l.startswith("#") or l == "":
                continue
            
            # Extract details out of line
            contig, source, feature, start, end, score, strand, \
                frame, attribute = sl = l.split("\t") # attribute is '8s_rRNA' or '18s_rRNA' etc
            ongoingCounts.setdefault(attribute, 0) # used for ID formatting
            ongoingCounts[attribute] += 1
            
            # Format unique IDs for this feature
            parentID = f"RNAmmer.{attribute}.{ongoingCounts[attribute]}"
            rrnaID = f"{parentID}.1"
            exonID = f"{rrnaID}.exon1"
            
            # Format each feature line
            parentLine = f"{contig}\t{source}\tncRNA_gene\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={parentID}\n"
            rrnaLine = f"{contig}\t{source}\trRNA\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={rrnaID};Parent={parentID}\n"
            exonLine = f"{contig}\t{source}\texon\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={exonID};Parent={parentID}\n"
            
            # Write to file
            fileOut.write(parentLine)
            fileOut.write(rrnaLine)
            fileOut.write(exonLine)
