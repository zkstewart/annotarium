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

import os, sys, zipfile
from glob import glob

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file

ZIP_SUFFIX = "_fastqc.zip"
DATA_NAME = "fastqc_data.txt"
FLAGS = { "pass": "PASS", "warn": "WARNING", "fail": "FAIL" }
STAT_COLUMNS = ["total_reads", "approx_bp", "low_quality_reads"]
NULL_VALUE = "."

# Cell colours for flag values
COLOURS = {
    "PASS": "#00CC66",
    "WARNING": "#FFCC00",
    "FAIL": "#FF0000",
    "SAME": "#00CC66",
    "DIFFERS": "#FFCC00"
}

def null_if_none(value):
    '''
    Parameters:
        value -- any value, or None where a statistic is unreported (e.g., low_quality_reads
                 from FastQC 0.13 onwards)
    Returns:
        value -- the value unchanged, or NULL_VALUE ('.') if it was None
    '''
    return NULL_VALUE if value is None else value

def find_fastqc_samples(directories):
    '''
    Locates FastQC results, where each sample's '_fastqc.zip' files (one for single-end,
    two for paired-end reads) are in a directory named after the sample.

    Parameters:
        directories -- a list of strings indicating the locations of sample directories
                       and/or directories containing them
    Returns:
        samples -- a dictionary where keys are sample names (i.e., directory names) and
                   values are lists of one or two '_fastqc.zip' files, forward first
    '''
    samples = {}
    for directory in directories:
        if glob(os.path.join(directory, f"*{ZIP_SUFFIX}")) != []:
            sampleDirectories = [ directory ]
        else:
            sampleDirectories = sorted({ os.path.dirname(f) for f in glob(os.path.join(directory, "*", f"*{ZIP_SUFFIX}")) })

        for sampleDirectory in sampleDirectories:
            sampleName = os.path.basename(sampleDirectory)
            zipFiles = sorted(glob(os.path.join(sampleDirectory, f"*{ZIP_SUFFIX}"))) # e.g., _1 before _2
            if sampleName in samples:
                raise ValueError(f"More than one directory is named '{sampleName}'; sample names must be unique")
            if len(zipFiles) > 2:
                raise ValueError(f"'{sampleDirectory}' has {len(zipFiles)} '{ZIP_SUFFIX}' files; expected one " +
                                 "(single-end) or two (paired-end) per sample directory")
            samples[sampleName] = zipFiles
    return samples

def parse_sample_names_tsv(fileName):
    '''
    Parameters:
        fileName -- a string indicating the location of a two column TSV pairing read
                    prefixes (column 1) to sample names (column 2)
    Returns:
        sampleNames -- a dictionary where keys are read prefixes and values are sample names
    '''
    sampleNames = {}
    with read_gz_file(fileName) as fileIn:
        for lineNum, line in enumerate(fileIn, start=1):
            if line.startswith("#") or line.strip() == "":
                continue
            sl = line.rstrip("\r\n").split("\t")
            if len(sl) != 2:
                raise ValueError(f"Line {lineNum} of '{fileName}' has {len(sl)} columns; expected 2 " +
                                 "(read prefix and sample name)")
            prefix, sampleName = sl[0].strip(), sl[1].strip()
            if prefix in sampleNames and sampleNames[prefix] != sampleName:
                raise ValueError(f"'{prefix}' is paired to more than one sample name in '{fileName}'")
            sampleNames[prefix] = sampleName
    return sampleNames

def approximate_bp(lengthLine):
    '''
    Parameters:
        lengthLine -- a string of a 'Sequence Length Distribution' line, like '35-39\\t1500.0'
    Returns:
        bp -- a float of the line's read count multiplied by its (midpoint) read length
    '''
    lengthRange, count = lengthLine.split("\t")
    if "-" in lengthRange:
        start, end = map(int, lengthRange.split("-"))
        length = (start + end) / 2
    else:
        length = int(lengthRange)
    return length * float(count)

def parse_fastqc_zip(zipFile):
    '''
    Parameters:
        zipFile -- a string indicating the location of a '_fastqc.zip' file
    Returns:
        stats -- a dictionary with structure like:
                 {
                     "total_reads": 84552130,
                     "approx_bp": 11202942050,
                     "low_quality_reads": 0,
                     "basic_statistics": "PASS",
                     "per_base_sequence_quality": "WARNING",
                     ...
                 }
                 where flag keys are in the order FastQC reports them; low_quality_reads
                 is None for FastQC 0.13 onwards, which no longer reports it
    '''
    with zipfile.ZipFile(zipFile, "r") as archive:
        dataNames = [ name for name in archive.namelist() if os.path.basename(name) == DATA_NAME ]
        if len(dataNames) != 1:
            raise ValueError(f"'{zipFile}' should contain one '{DATA_NAME}' file; is it a FastQC result?")
        contents = archive.read(dataNames[0]).decode("utf-8")

    stats = { "total_reads": None, "approx_bp": 0, "low_quality_reads": None }
    flags = {}
    module = None
    for line in contents.splitlines():
        if line.startswith(">>END_MODULE"):
            module = None
        elif line.startswith(">>"):
            module, flag = line[2:].rsplit("\t", maxsplit=1)
            flags[module.lower().replace(" ", "_")] = FLAGS.get(flag, flag.upper())
        elif line.startswith("#") or line.strip() == "":
            continue
        elif module == "Basic Statistics":
            measure, value = line.split("\t", maxsplit=1)
            if measure == "Total Sequences":
                stats["total_reads"] = int(value)
            elif measure == "Sequences flagged as poor quality":
                stats["low_quality_reads"] = int(value)
        elif module == "Sequence Length Distribution":
            stats["approx_bp"] += approximate_bp(line)

    if stats["total_reads"] is None:
        raise ValueError(f"'{zipFile}' does not report 'Total Sequences'; is it a FastQC result?")
    stats["approx_bp"] = int(stats["approx_bp"])
    stats.update(flags)
    return stats

def compare_reads(forward, reverse, flagColumns):
    '''
    Combines a sample's forward and reverse read statistics.

    Parameters:
        forward -- a dictionary as returned by parse_fastqc_zip()
        reverse -- a dictionary as returned by parse_fastqc_zip(), or None for single-end reads
        flagColumns -- a list of the flag keys to compare
    Returns:
        combined -- a dictionary of summed read statistics, and flag keys with values of
                    'SAME' or 'DIFFERS'; values are '.' where reads are single-end, a
                    flag is absent from both files, or a statistic is unreported (None)
    '''
    if reverse is None:
        combined = { column: null_if_none(forward[column]) for column in STAT_COLUMNS }
        combined.update({ column: NULL_VALUE for column in flagColumns })
        return combined

    combined = { column: NULL_VALUE if forward[column] is None or reverse[column] is None else forward[column] + reverse[column]
                 for column in STAT_COLUMNS }
    for column in flagColumns:
        forwardFlag, reverseFlag = forward.get(column), reverse.get(column)
        if forwardFlag is None and reverseFlag is None:
            combined[column] = NULL_VALUE
        else:
            combined[column] = "SAME" if forwardFlag == reverseFlag else "DIFFERS"
    return combined

def write_fastqc_workbook(sheets, outputFileName):
    '''
    Parameters:
        sheets -- a list of (sheetName, header, rows) tuples, where header is a list of
                  column names and rows is a list of lists of values; columns whose values
                  are flags (e.g., 'PASS' or 'SAME') are coloured accordingly
        outputFileName -- a string indicating the location to write the .xlsx file to
    '''
    import xlsxwriter # only needed by this mode, so not a requirement of annotarium as a whole

    workbook = xlsxwriter.Workbook(outputFileName)
    headerFormat = workbook.add_format({ "bold": True, "border": 1 })
    colourFormats = { value: workbook.add_format({ "bg_color": colour }) for value, colour in COLOURS.items() }

    for sheetName, header, rows in sheets:
        worksheet = workbook.add_worksheet(sheetName)
        worksheet.write_row(0, 0, header, headerFormat)
        for rowIndex, row in enumerate(rows, start=1):
            worksheet.write_row(rowIndex, 0, row)

        # Size columns to fit their contents
        for columnIndex, column in enumerate(header):
            width = max([ len(str(column)) ] + [ len(str(row[columnIndex])) for row in rows ])
            worksheet.set_column(columnIndex, columnIndex, width + 1)
        worksheet.freeze_panes(1, 1)

        # Colour flag values in all columns other than the first (i.e., the row labels)
        if rows != []:
            for value, colourFormat in colourFormats.items():
                worksheet.conditional_format(1, 1, len(rows), len(header) - 1,
                                             { "type": "cell", "criteria": "==",
                                               "value": f'"{value}"', "format": colourFormat })
    workbook.close()
