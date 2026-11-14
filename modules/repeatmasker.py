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

import os, re, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, GzCapableWriter

TOKEN_REGEX = re.compile(r"\S+")
REPEAT_INDEX = 9 # 0-based index of the 'matching repeat' column
CLASS_INDEX = 10 # 0-based index of the 'repeat class/family' column
BEGIN_INDEX = 11 # 0-based index of the 'position in repeat begin' column
MIN_COLUMNS = 15 # an optional 16th column holds an asterisk for overlapping matches
HEADER_LINES = 3 # two lines of column labels then a blank line
NO_REPEATS_MESSAGE = "There were no repetitive sequences detected"

def parse_class_family_tsv(fileName):
    '''
    Parameters:
        fileName -- a string indicating the location of a two column TSV pairing
                    repeat IDs to class/family labels
    Returns:
        classFamilies -- a dictionary with structure like:
                         {
                             repeatID1: "LTR/Copia",
                             repeatID2: "DNA/hAT",
                             ...
                         }
    '''
    classFamilies = {}
    with read_gz_file(fileName) as fileIn:
        for lineNum, line in enumerate(fileIn, start=1):
            if line.startswith("#") or line.strip() == "":
                continue
            sl = line.rstrip("\r\n").split("\t")
            if len(sl) != 2:
                raise ValueError(f"Line {lineNum} of '{fileName}' has {len(sl)} columns; expected 2 " +
                                 "(repeat ID and class/family)")
            repeatID, classFamily = sl[0].strip(), sl[1].strip()
            if repeatID == "" or classFamily == "":
                raise ValueError(f"Line {lineNum} of '{fileName}' has an empty value")
            if TOKEN_REGEX.fullmatch(repeatID) is None or TOKEN_REGEX.fullmatch(classFamily) is None:
                raise ValueError(f"Line {lineNum} of '{fileName}' has a value containing whitespace, " +
                                 "which RepeatMasker .out files cannot hold")
            if repeatID in classFamilies and classFamilies[repeatID] != classFamily:
                raise ValueError(f"'{repeatID}' is paired to more than one class/family in '{fileName}' " +
                                 f"('{classFamilies[repeatID]}' and '{classFamily}')")
            classFamilies[repeatID] = classFamily
    return classFamilies

def iterate_out_lines(outFile):
    '''
    Parameters:
        outFile -- a string indicating the location of a RepeatMasker .out file
    Yields:
        line -- the line as a string, including its newline
        tokens -- None for header or blank lines; otherwise, a list of re.Match objects
                  for each whitespace-delimited value on the line
    '''
    with read_gz_file(outFile) as fileIn:
        for lineNum, line in enumerate(fileIn, start=1):
            tokens = list(TOKEN_REGEX.finditer(line))

            # Header lines start with 'SW' and 'score'; data lines start with an integer score
            if tokens == [] or not tokens[0].group().isdigit():
                if tokens != [] and lineNum > HEADER_LINES and not line.startswith(NO_REPEATS_MESSAGE):
                    raise ValueError(f"Line {lineNum} of '{outFile}' is not a RepeatMasker match line " +
                                     f"('{line.rstrip()}'). Is this a RepeatMasker .out file?")
                yield line, None
                continue

            if len(tokens) < MIN_COLUMNS:
                raise ValueError(f"Line {lineNum} of '{outFile}' has {len(tokens)} columns; " +
                                 f"expected at least {MIN_COLUMNS}. Is this a RepeatMasker .out file?")
            yield line, tokens

def shift_header_line(line, anchor, shift):
    '''
    Inserts spaces before an anchor word in a header line so that header labels stay
    above their columns after the class/family column is widened.

    Parameters:
        line -- a string of the header line
        anchor -- a string of the word to insert spaces before; its last occurrence is used
        shift -- an integer number of spaces to insert
    Returns:
        line -- the header line with spaces inserted, or unchanged if the anchor is absent
    '''
    index = line.rfind(anchor)
    if shift == 0 or index == -1:
        return line
    return line[:index] + " "*shift + line[index:]

def repeatmasker_classify(args):
    classFamilies = parse_class_family_tsv(args.classFamilyFile)
    
    # First pass: find missing IDs and the column width needed for new labels
    missing = {}
    numMatches = 0
    originalWidth = None # width from the start of class/family to the end of 'begin'
    regionWidth = 0
    for line, tokens in iterate_out_lines(args.outFile):
        if tokens is None:
            continue
        numMatches += 1
        repeatID = tokens[REPEAT_INDEX].group()
        oldClass, begin = tokens[CLASS_INDEX], tokens[BEGIN_INDEX]
        if not repeatID in classFamilies:
            missing[repeatID] = missing.get(repeatID, 0) + 1
        newClass = classFamilies.get(repeatID, oldClass.group())
        
        if originalWidth is None:
            originalWidth = begin.end() - oldClass.start()
        regionWidth = max(regionWidth,
                          begin.end() - oldClass.start(), # keep existing alignment
                          len(newClass) + 1 + len(begin.group())) # fit the new label
    headerShift = 0 if originalWidth is None else regionWidth - originalWidth
    
    # Stop if missing IDs aren't tolerated
    if missing != {} and not args.allowMissing:
        examples = ", ".join(f"'{x}'" for x in list(missing)[:5])
        raise KeyError(f"{len(missing)} repeat ID(s) in '{args.outFile}' are not in " +
                       f"'{args.classFamilyFile}' (e.g., {examples}); add them to the TSV, " +
                       "or use --allowMissing to leave their class/family unchanged")
    
    # Second pass: write relabelled output
    numChanged = 0
    with GzCapableWriter(args.outputFileName) as fileOut:
        for line, tokens in iterate_out_lines(args.outFile):
            # Keep header labels above their columns if class/family was widened
            if tokens is None:
                if "position in repeat" in line:
                    line = shift_header_line(line, "position in repeat", headerShift)
                elif "class/family" in line:
                    line = shift_header_line(line, "begin", headerShift)
                fileOut.write(line)
                continue
            
            repeatID = tokens[REPEAT_INDEX].group()
            oldClass, begin = tokens[CLASS_INDEX], tokens[BEGIN_INDEX]
            newClass = classFamilies.get(repeatID, oldClass.group())
            if newClass != oldClass.group():
                numChanged += 1
            
            # Left-align class/family and right-align 'begin' within the region
            region = newClass.ljust(regionWidth - len(begin.group())) + begin.group()
            fileOut.write(line[:oldClass.start()] + region + line[begin.end():])
    
    # Report on results
    print(f"Relabelled {numChanged} of {numMatches} repeat matches")
    if missing != {}:
        print(f"{len(missing)} repeat ID(s) ({sum(missing.values())} matches) were not in the TSV " +
              "and kept their original class/family")
