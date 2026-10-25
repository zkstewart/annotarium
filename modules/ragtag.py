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

import os, subprocess, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, GzCapableWriter
from fasta import FASTATarium, Sequence

class RagTag:
    '''
    A wrapper for performing reference-guided scaffolding with RagTag. Code helpfully
    adapted from Various_scripts' phase_contigs.py with modification.
    '''
    SUFFIX = "_RagTag"

    def __init__(self, ragtagPath, minimap2Path=None, threads=1):
        self.exe = ragtagPath
        self.minimap2 = minimap2Path
        self.threads = threads

    @property
    def exe(self):
        return self._exe

    @exe.setter
    def exe(self, value):
        if not os.path.isfile(value):
            raise FileNotFoundError(f"ragtag.py not found at '{value}'")
        self._exe = value

    @property
    def minimap2(self):
        return self._minimap2

    @minimap2.setter
    def minimap2(self, value):
        if value is not None and not os.path.isfile(value):
            raise FileNotFoundError(f"minimap2 executable not found at '{value}'")
        self._minimap2 = value

    @property
    def threads(self):
        return self._threads

    @threads.setter
    def threads(self, value):
        if not isinstance(value, int) or value < 1:
            raise ValueError("threads must be an integer >= 1")
        self._threads = value

    def scaffold(self, referenceFile, queryFile, outputDirectory):
        '''
        Parameters:
            referenceFile -- a string indicating the location of the reference FASTA file
            queryFile -- a string indicating the location of the FASTA file to scaffold
            outputDirectory -- a string indicating the directory to write RagTag outputs to
        Returns:
            fastaFile -- a string indicating the location of the 'ragtag.scaffold.fasta' file
            agpFile -- a string indicating the location of the 'ragtag.scaffold.agp' file
        '''
        cmd = [self.exe, "scaffold", "-t", str(self.threads), "-o", outputDirectory]
        if self.minimap2 is not None:
            cmd += ["--aligner", self.minimap2]
        cmd += [referenceFile, queryFile]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0 or not "Finished running" in result.stderr:
            raise RuntimeError(("RagTag encountered an error; have a look at the stderr " +
                                f"({result.stderr}) to make sense of this."))

        return (os.path.join(outputDirectory, "ragtag.scaffold.fasta"),
                os.path.join(outputDirectory, "ragtag.scaffold.agp"))

    @staticmethod
    def rename_scaffolds(ragtagFastaFile, outputFileName, lineLength=None):
        '''
        Rewrites a RagTag scaffold FASTA such that scaffolds are named after their reference
        contig (i.e., without the '_RagTag' suffix). Unplaced contigs retain their names.

        Parameters:
            ragtagFastaFile -- a string indicating the location of a 'ragtag.scaffold.fasta' file
            outputFileName -- a string indicating the location to write the renamed FASTA to
            lineLength -- OPTIONAL; an integer of characters per line, or None for no multilining
        '''
        def renamed_records():
            with read_gz_file(ragtagFastaFile) as fileIn:
                for record in FASTATarium.fasta_parser(fileIn):
                    if record.id.endswith(RagTag.SUFFIX):
                        record = Sequence(record.id[:-len(RagTag.SUFFIX)], record.sequence)
                    yield record

        # Make sure renaming doesn't cause ID collisions before writing anything
        seenIDs = set()
        for record in renamed_records():
            if record.id in seenIDs:
                raise ValueError(f"Removing the '{RagTag.SUFFIX}' suffix from scaffolds in '{ragtagFastaFile}' " +
                                 f"would cause '{record.id}' to occur more than once")
            seenIDs.add(record.id)

        with GzCapableWriter(outputFileName) as fileOut:
            for record in renamed_records():
                fileOut.write(record.format(lineLength))

    @staticmethod
    def parse_agp(agpFile):
        '''
        Parameters:
            agpFile -- a string indicating the location of a 'ragtag.scaffold.agp' file
        Returns:
            scaffolds -- a dictionary where keys are scaffold IDs (without the '_RagTag'
                         suffix) and values are lists of the contig IDs they contain, in order;
                         unplaced contigs are keyed by their own ID
            placed -- a set of the scaffold IDs which were placed onto a reference contig
        '''
        scaffolds = {}
        placed = set()
        with open(agpFile, "r") as fileIn:
            for line in fileIn:
                if line.startswith("#") or line.strip() == "":
                    continue
                sl = line.rstrip("\r\n").split("\t")
                scaffoldID, componentType, componentID = sl[0], sl[4], sl[5]

                if scaffoldID.endswith(RagTag.SUFFIX):
                    scaffoldID = scaffoldID[:-len(RagTag.SUFFIX)]
                    placed.add(scaffoldID)

                # Skip gaps
                if componentType in ["N", "U"]:
                    continue
                scaffolds.setdefault(scaffoldID, []).append(componentID)
        return scaffolds, placed
