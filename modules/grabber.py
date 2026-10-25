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
from pyfaidx import Fasta

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import GzCapableWriter

def index_names(records):
    '''
    pyfaidx keys records by their ID (everything up to the first whitespace). This
    function additionally indexes each record by its full description, so that a
    lookup can be performed using either value.

    Note that pyfaidx's .long_name retains a trailing '\r' for files with Windows-style
    newlines, which is removed here. Also note that a description cannot collide with
    a different record's ID, since a description lacking whitespace is equal to its
    own ID; and pyfaidx will raise an error if duplicate IDs exist.

    Parameters:
        records -- a pyfaidx.Fasta object
    Returns:
        nameIndex -- a dictionary where keys are record IDs or descriptions, and values
                     are the record ID as used by pyfaidx
    '''
    nameIndex = {}
    for record in records:
        description = record.long_name.rstrip("\r")
        nameIndex[record.name] = record.name
        nameIndex[description] = record.name
    return nameIndex

# fasta > grabber
def fasta_grabber(args):
    '''
    Parameters:
        args.sequences -- a list of sequence IDs and/or descriptions to retrieve or remove
        args.behaviour -- a string of 'retrieve' or 'remove'
    '''
    records = Fasta(args.fastaFile)
    nameIndex = index_names(records)

    # Locate the requested sequences
    foundIDs = {} # dict acts as an ordered set
    missingValues = []
    for value in args.sequences:
        if value in nameIndex:
            foundIDs[nameIndex[value]] = None
        else:
            missingValues.append(value)

    # Report on any values that were not found
    if missingValues != []:
        exampleValues = ", ".join([ f"'{x}'" for x in missingValues[:10] ]) + (", ..." if len(missingValues) > 10 else "")
        if args.behaviour == "retrieve":
            raise ValueError(f"{len(missingValues)} value(s) given to -s were not found in the FASTA file " +
                             f"as a sequence ID or description; these include {exampleValues}")
        else:
            print(f"# WARNING: {len(missingValues)} value(s) given to -s were not found in the FASTA file " +
                  f"as a sequence ID or description; these include {exampleValues}")

    # Write sequences to file
    if args.behaviour == "retrieve":
        outputIDs = list(foundIDs) # in the order they were given
    else:
        outputIDs = [ record.name for record in records if not record.name in foundIDs ] # in file order

    with GzCapableWriter(args.outputFileName) as fileOut:
        for seqID in outputIDs:
            record = records[seqID]
            description = record.long_name.rstrip("\r")
            fileOut.write(f">{description}\n{str(record)}\n")
    print(f"# {len(outputIDs)} sequence(s) written to output")
