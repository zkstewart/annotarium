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

def describe_target(target):
    'Formats a (fastaBasename, value) target for messages'
    fastaName, value = target
    return f"'{value}'" if fastaName is None else f"'{value}' (in {fastaName})"

def locate_targets(targets, fastaNames, nameIndexes):
    '''
    Parameters:
        targets -- a list of (fastaBasename, value) tuples where fastaBasename is None
                   if the value may come from any FASTA file
        fastaNames -- a list of FASTA file basenames, in -i order
        nameIndexes -- a list of dictionaries as returned by index_names(), in -i order
    Returns:
        foundKeys -- a dictionary (acting as an ordered set) of (fileIndex, seqID) keys
        missingTargets -- a list of targets which were not found
    '''
    fileIndices = { fastaName: fileIndex for fileIndex, fastaName in enumerate(fastaNames) }
    foundKeys = {}
    missingTargets = []
    ambiguousTargets = []
    for target in targets:
        fastaName, value = target
        if fastaName is None:
            candidates = [ fileIndex for fileIndex, nameIndex in enumerate(nameIndexes) if value in nameIndex ]
        else:
            fileIndex = fileIndices[fastaName]
            candidates = [ fileIndex ] if value in nameIndexes[fileIndex] else []
        
        if len(candidates) == 0:
            missingTargets.append(target)
        elif len(candidates) > 1:
            ambiguousTargets.append(target)
        else:
            fileIndex = candidates[0]
            foundKeys[(fileIndex, nameIndexes[fileIndex][value])] = None
    
    if ambiguousTargets != []:
        exampleValues = ", ".join([ describe_target(x) for x in ambiguousTargets[:10] ]) + (", ..." if len(ambiguousTargets) > 10 else "")
        raise ValueError(f"{len(ambiguousTargets)} value(s) given to -s were found in more than one -i FASTA " +
                         f"file; these include {exampleValues}. Use a two column TSV of 'fasta_file' and " +
                         "'contig_id' values to specify which file each should come from")
    return foundKeys, missingTargets

# fasta > grabber
def fasta_grabber(args):
    '''
    Parameters:
        args.fastaFiles -- a list of FASTA file locations
        args.targets -- a list of (fastaBasename, value) tuples where values are sequence IDs
                        or descriptions to retrieve or remove, and fastaBasename is None if the
                        value may come from any FASTA file
        args.behaviour -- a string of 'retrieve' or 'remove'
    '''
    fastaNames = [ os.path.basename(fastaFile) for fastaFile in args.fastaFiles ]
    records = [ Fasta(fastaFile) for fastaFile in args.fastaFiles ]
    nameIndexes = [ index_names(fileRecords) for fileRecords in records ]

    # Locate the requested sequences
    foundKeys, missingTargets = locate_targets(args.targets, fastaNames, nameIndexes)

    # Report on any values that were not found
    if missingTargets != []:
        exampleValues = ", ".join([ describe_target(x) for x in missingTargets[:10] ]) + (", ..." if len(missingTargets) > 10 else "")
        message = (f"{len(missingTargets)} value(s) given to -s were not found in the FASTA file(s) " +
                   f"as a sequence ID or description; these include {exampleValues}")
        if args.behaviour == "retrieve":
            raise ValueError(message)
        else:
            print(f"# WARNING: {message}")

    # Determine which sequences to write
    if args.behaviour == "retrieve":
        outputKeys = list(foundKeys) # in the order they were given
    else:
        outputKeys = [ (fileIndex, record.name) # in -i order then file order
                       for fileIndex, fileRecords in enumerate(records)
                       for record in fileRecords
                       if not (fileIndex, record.name) in foundKeys ]

    # Make sure the output won't contain duplicate sequence IDs
    seenIDs = {}
    for fileIndex, seqID in outputKeys:
        if seqID in seenIDs:
            raise ValueError(f"'{seqID}' would be written more than once (from '{fastaNames[seenIDs[seqID]]}' " +
                             f"and '{fastaNames[fileIndex]}'); output sequence IDs must be unique")
        seenIDs[seqID] = fileIndex

    # Write sequences to file
    with GzCapableWriter(args.outputFileName) as fileOut:
        for fileIndex, seqID in outputKeys:
            record = records[fileIndex][seqID]
            description = record.long_name.rstrip("\r")
            fileOut.write(f">{description}\n{str(record)}\n")
    print(f"# {len(outputKeys)} sequence(s) written to output")
