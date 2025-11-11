#! python3

import os, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, GzCapableWriter

def irf_to_gff3(args):
    ongoingCount = 0
    with read_gz_file(args.irfDatFile) as fileIn, GzCapableWriter(args.outputFileName) as fileOut:
        for line in fileIn:
            # Parse line
            l = line.strip()
            sl = l.split()
            
            # Handle sequence ID line
            if l.startswith("Sequence: "):
                contig = l.split(": ")[-1].strip()
                continue
            
            # Skip if irrelevant
            if l == "" or not sl[0].isdigit():
                continue
            
            # Extract details out of line
            lstart, lstop, llen, rstart, rstop, rlen, gaplen, \
                identity, indels, score, atpct, cgpct, atcomppct, \
                gccomppct, gtpct, centrex2, avgcentrex2, left, right = sl
            llen, rlen, gaplen = int(llen), int(rlen), int(gaplen)
            identity = float(identity)
            
            # Filter if needed
            if args.minimumLength > 0:
                irLength = min(llen, rlen)
                if irLength < args.minimumLength:
                    continue
            if args.maximumLength > 0:
                irLength = max(llen, rlen)
                if irLength > args.maximumLength:
                    continue
            if args.minimumGap > 0:
                if gaplen < args.minimumGap:
                    continue
            if args.maximumGap > 0:
                if gaplen > args.maximumGap:
                    continue
            if args.identityCutoff > 0:
                if identity < args.identityCutoff:
                    continue
            
            # Format GFF3 attributes
            ongoingCount += 1
            thisID = f"{contig}.TIR.{ongoingCount}"
            thisIdentity = "{:.2f}".format(identity)
            attributes = f"ID={thisID};identity={thisIdentity}"
            
            # Format each feature line
            parentLine = f"{contig}\tIRF\tinverted_repeat\t{lstart}\t{rstop}\t{score}\t.\t.\tID={thisID};identity={thisIdentity}\n"
            leftLine = f"{contig}\tIRF\tleft_repeat\t{lstart}\t{lstop}\t.\t+\t.\tID={thisID}.left;Parent={thisID}\n"
            rightLine = f"{contig}\tIRF\tright_repeat\t{rstart}\t{rstop}\t.\t-\t.\tID={thisID}.right;Parent={thisID}\n"
            
            # Write to file
            fileOut.write(parentLine)
            fileOut.write(leftLine)
            fileOut.write(rightLine)
