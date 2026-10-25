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
