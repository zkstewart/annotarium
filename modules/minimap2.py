#! python3

import os, subprocess, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file

class Minimap2:
    '''
    A wrapper for performing assembly-to-assembly alignment with minimap2. Code helpfully
    adapted from Various_scripts' ZS_AlignIO with modification.
    '''
    ACCEPTED_PRESETS = ["asm5", "asm10", "asm20"]

    def __init__(self, minimap2Path, preset="asm20", threads=1):
        self.exe = minimap2Path
        self.preset = preset
        self.threads = threads

    @property
    def exe(self):
        return self._exe

    @exe.setter
    def exe(self, value):
        if not os.path.isfile(value):
            raise FileNotFoundError(f"minimap2 executable not found at '{value}'")
        self._exe = value

    @property
    def preset(self):
        return self._preset

    @preset.setter
    def preset(self, value):
        if not value in Minimap2.ACCEPTED_PRESETS:
            raise ValueError(f"preset '{value}' not recognised; please choose from {Minimap2.ACCEPTED_PRESETS}")
        self._preset = value

    @property
    def threads(self):
        return self._threads

    @threads.setter
    def threads(self, value):
        if not isinstance(value, int) or value < 1:
            raise ValueError("threads must be an integer >= 1")
        self._threads = value

    def align(self, targetFile, queryFiles, outputFileName):
        '''
        Parameters:
            targetFile -- a string indicating the location of the reference FASTA file
            queryFiles -- a list of strings indicating the locations of query FASTA files
            outputFileName -- a string indicating the location to write the PAF file to
        '''
        cmd = [
            self.exe, "-cx", self.preset, "-t", str(self.threads),
            "-o", outputFileName, targetFile
        ] + queryFiles

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(("minimap2 encountered an error; have a look at the stderr " +
                                f"({result.stderr}) to make sense of this."))

    def __repr__(self):
        return "<Minimap2(exe='{0}',preset='{1}',threads={2})>".format(
            self.exe, self.preset, self.threads
        )

def parse_paf(pafFile, minQ=0, minAlignLen=0):
    '''
    Parses a minimap2 PAF file, retaining primary (or supplementary) alignments which meet
    the MAPQ and alignment length minimums. Secondary alignments are always skipped.

    Parameters:
        pafFile -- a string indicating the location of a PAF file
        minQ -- OPTIONAL; an integer of the minimum MAPQ for an alignment to be kept
        minAlignLen -- OPTIONAL; an integer of the minimum alignment block length (column 11)
                       for an alignment to be kept
    Returns:
        alignments -- a dictionary where keys are query IDs and values are lists of tuples
                      with structure like (qstart, qend, tid, tstart, tend)
    '''
    alignments = {}
    with read_gz_file(pafFile) as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            if len(sl) < 12:
                continue
            qid, qlen, qstart, qend, strand, tid, tlen, tstart, tend, \
                numMatches, alignLen, mapq = sl[0:12]

            # Skip alignments that fail our criteria
            if "tp:A:S" in sl[12:]:
                continue
            if int(mapq) < minQ:
                continue
            if int(alignLen) < minAlignLen:
                continue

            alignments.setdefault(qid, []).append((int(qstart), int(qend), tid, int(tstart), int(tend)))
    return alignments
