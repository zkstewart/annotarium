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

import io, os, subprocess, sys, tempfile

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, GzCapableWriter
from fasta import FASTATarium, Sequence

class MAFFT:
    '''
    A wrapper for performing plain MAFFT alignment of a FASTA file. Code helpfully
    copied from Various_scripts' ZS_AlignIO with modification.
    '''
    ALG_CONVERT = {
            "auto": ["--auto"],
            "einsi": ["--genafpair"],
            "linsi": ["--localpair"],
            "ginsi": ["--globalpair"],
            "fftns1": ["--retree", "1"],
            "fftns2": ["--retree", "2"],
            "fftnsi": []
        }
    ACCEPTED_ALGS = ["auto", "einsi", "linsi", "ginsi", "fftns1", "fftns2", "fftnsi"]
    ACCEPTED_MOLECULES = ["auto", "nucleotide", "protein"]

    def __init__(self, mafftPath, algorithm="auto", thread=1, maxiterate=0):
        self.exe = mafftPath
        self.algorithm = algorithm
        self.thread = thread
        self.maxiterate = maxiterate

    @property
    def exe(self):
        return self._exe

    @exe.setter
    def exe(self, value):
        if not os.path.isfile(value):
            raise FileNotFoundError(f"MAFFT executable not found at '{value}'")
        self._exe = value

    @property
    def algorithm(self):
        return self._algorithm

    @algorithm.setter
    def algorithm(self, value):
        if not value in MAFFT.ACCEPTED_ALGS:
            raise ValueError(f"algorithm '{value}' not recognised; please choose from {MAFFT.ACCEPTED_ALGS}")
        self._algorithm = value

    @property
    def thread(self):
        return self._thread

    @thread.setter
    def thread(self, value):
        if not isinstance(value, int) or value < 1:
            raise ValueError("thread must be an integer >= 1")
        self._thread = value

    @property
    def maxiterate(self):
        return self._maxiterate

    @maxiterate.setter
    def maxiterate(self, value):
        if not isinstance(value, int) or value < 0:
            raise ValueError("maxiterate must be an integer >= 0")
        self._maxiterate = value

    def _align(self, fastaFile, molecule="auto"):
        '''
        Runs MAFFT on a FASTA file and parses its stdout.

        Parameters:
            fastaFile -- a string indicating the location of an uncompressed FASTA file to align
            molecule -- OPTIONAL; a string of 'auto', 'nucleotide', or 'protein' to indicate
                        whether MAFFT should detect or be told the molecule type
        Returns:
            alignedSeqs -- a list of Sequence objects containing aligned sequences in input order
        '''
        if not molecule in MAFFT.ACCEPTED_MOLECULES:
            raise ValueError(f"molecule '{molecule}' not recognised; please choose from {MAFFT.ACCEPTED_MOLECULES}")
        if self.algorithm == "fftnsi" and self.maxiterate < 2:
            raise ValueError("maxiterate must be >= 2 for the 'fftnsi' algorithm")

        # Construct the cmd for subprocess
        cmd = [self.exe] + MAFFT.ALG_CONVERT[self.algorithm]
        if self.algorithm != "auto":
            cmd += ["--maxiterate", str(self.maxiterate)]
        if molecule == "nucleotide":
            cmd.append("--nuc")
        elif molecule == "protein":
            cmd.append("--amino")
        cmd += [
            "--preservecase", # MAFFT otherwise lowercases nucleotide output
            "--quiet", # necessary to prevent unusual MAFFT errors; see https://mafft.cbrc.jp/alignment/software/stderr.html
            "--thread", str(self.thread),
            fastaFile
        ]

        # Run the command
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0 or not result.stdout.startswith(">"):
            raise RuntimeError(("MAFFT encountered an error; have a look at the stdout " +
                                f"({result.stdout}) and stderr ({result.stderr}) to make sense of this."))

        # Parse out the MSA result
        return list(FASTATarium.fasta_parser(io.StringIO(result.stdout.replace("\r", ""))))

    def align(self, records, asCodons=False):
        '''
        Parameters:
            records -- a list of Sequence objects to align
            asCodons -- OPTIONAL; a boolean indicating whether to align nucleotides as codons,
                        by translating them in the +ve strand and frame 0 before alignment and
                        then mapping the aligned proteins back to their codons
        Returns:
            alignedSeqs -- a list of Sequence objects containing aligned sequences in input order
        '''
        if asCodons:
            toAlign = MAFFT.codons_to_proteins(records)
            molecule = "protein"
        else:
            toAlign = records
            molecule = "auto"

        # Write sequences to a temporary file for MAFFT
        tmpHandle, tmpFile = tempfile.mkstemp(suffix=".fasta")
        try:
            with os.fdopen(tmpHandle, "w") as fileOut:
                for record in toAlign:
                    fileOut.write(record.format())
            alignedSeqs = self._align(tmpFile, molecule=molecule)
        finally:
            os.unlink(tmpFile)

        if asCodons:
            alignedSeqs = MAFFT.proteins_to_codons(records, alignedSeqs)
        return alignedSeqs

    @staticmethod
    def codons_to_proteins(records):
        '''
        Translates nucleotide sequences in the +ve strand and frame 0 for codon alignment.
        Only complete codons are translated, so that each residue corresponds to exactly
        3 nucleotides when mapping back with proteins_to_codons().

        Parameters:
            records -- a list of Sequence objects containing nucleotides
        Returns:
            proteins -- a list of Sequence objects containing translated proteins
        '''
        proteins = []
        for record in records:
            completeLength = len(record) - (len(record) % 3)
            protein = record[0:completeLength].translate()

            '''
            MAFFT will NOT report stop codons in the output sequence, which poses problems
            when mapping codon positions back. Since proteins_to_codons() works by consuming the
            next 3 nucleotides per aligned residue regardless of what the residue is, we can
            just substitute *'s for X's.
            '''
            protein.sequence = protein.sequence.replace("*", "X")
            proteins.append(protein)
        return proteins

    @staticmethod
    def proteins_to_codons(records, alignedProteins):
        '''
        Maps aligned proteins back to their original nucleotide sequences, such that
        each aligned residue becomes its codon and each gap becomes '---'. Any leftover
        nucleotides from an incomplete final codon are placed after the last codon.

        Parameters:
            records -- a list of Sequence objects containing the original nucleotides
            alignedProteins -- a list of Sequence objects containing the aligned proteins,
                               in the same order as records
        Returns:
            alignedSeqs -- a list of Sequence objects containing aligned nucleotides
        '''
        alignedSeqs = []
        for record, alignedProtein in zip(records, alignedProteins):
            if record.id != alignedProtein.id:
                raise ValueError(f"MAFFT output order differs from the input; expected '{record.id}' " +
                                 f"but found '{alignedProtein.id}'")
            nuc = str(record)

            # Perform mapping procedure
            alignedNuc = []
            codonIndex = 0
            for residue in str(alignedProtein):
                if residue != "-":
                    alignedNuc.append(nuc[codonIndex:codonIndex+3])
                    codonIndex += 3
                else:
                    alignedNuc.append("---")

            # Add in any leftover bits as close to the sequence end as possible
            alignedNuc = "".join(alignedNuc).rstrip("-") + nuc[codonIndex:]
            alignedSeqs.append(Sequence(record.description, alignedNuc))

        # Pad out the end of sequences for anything we rstripped
        maxLen = max([ len(seq) for seq in alignedSeqs ])
        for seq in alignedSeqs:
            seq.sequence += "-" * (maxLen - len(seq))

        return alignedSeqs

    def __repr__(self):
        return "<MAFFT(exe='{0}',algorithm='{1}',thread={2},maxiterate={3})>".format(
            self.exe, self.algorithm, self.thread, self.maxiterate
        )

# fasta > to > msa
def fasta_to_msa(args):
    with read_gz_file(args.fastaFile) as fileIn:
        records = list(FASTATarium.fasta_parser(fileIn))

    aligner = MAFFT(args.mafft, algorithm=args.algorithm,
                    thread=args.threads, maxiterate=args.maxiterate)
    alignedSeqs = aligner.align(records, asCodons=args.asCodons)

    with GzCapableWriter(args.outputFileName) as fileOut:
        for seq in alignedSeqs:
            fileOut.write(seq.format())
