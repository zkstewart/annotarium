#!/usr/bin/env python3
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

import os, random, re, sys, unittest, shutil, subprocess, tempfile, zipfile

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.fasta import FASTATarium, Sequence
from modules.msa import MAFFT
from modules.annotable import Annotable
from modules.minimap2 import parse_paf
from modules.parsing import parse_phasing_tsv, parse_phasing_matches
from modules.ragtag import RagTag
from modules.phasing import HaplotypePartitioner, deduplicate_query, merge_intervals, \
    build_footprints, infer_orientations, unique_coverage, group_coverage, pair_homologous_contigs
from modules.repeatmasker import parse_class_family_tsv
from modules.compleasm import CompleasmTable, SetOptimiser, parse_compleasm_summary
from modules.parsing import parse_file_contig_tsv, parse_fasta_stats_file, find_files_with_suffix
from modules.fastqc import parse_fastqc_zip, compare_reads

# Specify data locations
baseDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
dataDir = os.path.join(baseDir, "tests", "data")
codonIndelsFile = os.path.join(dataDir, "codon_indels.fasta")
annotableFile = os.path.join(dataDir, "annotable.tsv")
pafFile = os.path.join(dataDir, "phasing.paf")
agpFile = os.path.join(dataDir, "ragtag.agp")
repeatmaskerFile = os.path.join(dataDir, "repeatmasker.out")

MAFFT_EXE = shutil.which("mafft")
MINIMAP2_EXE = shutil.which("minimap2")
RAGTAG_EXE = shutil.which("ragtag.py")

try:
    import xlsxwriter, openpyxl
    EXCEL_AVAILABLE = True
except ImportError:
    EXCEL_AVAILABLE = False

# Define utility functions
def run_subprocess(command):
    '''
    Parameters:
        command -- a list of strings representing the command to run
    '''
    process = subprocess.Popen(command,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE)
    stdout, stderr = process.communicate()
    return process.returncode, stdout.decode("utf-8").rstrip("\r\n "), stderr.decode("utf-8").rstrip("\r\n ")

def parse_fasta(fastaFile):
    '''
    Returns:
        records -- a list of Sequence objects in file order
    '''
    with open(fastaFile, "r") as fileIn:
        return list(FASTATarium.fasta_parser(fileIn))

def run_annotable_filter(outputFile, *filterArgs):
    '''
    Parameters:
        outputFile -- a string indicating the location to write the filtered annotable
        filterArgs -- strings to pass as arguments to 'annotable filter'
    '''
    cmd = [
        sys.executable, os.path.join(baseDir, "annotarium.py"), "annotable", "filter",
        "-i", annotableFile,
        "-o", outputFile,
        *filterArgs
    ]
    return run_subprocess(cmd)

def write_phasing_inputs(tmpDir):
    '''
    Writes a random reference with two chromosomes, and contigs from two haplotypes that
    each carry SNPs relative to the reference. Haplotype 'A' has one contig per chromosome,
    whereas haplotype 'B' has chromosome 1 broken into two contigs.

    Returns:
        referenceFile -- a string indicating the location of the reference FASTA
        contigsFile -- a string indicating the location of the contigs FASTA
    '''
    rng = random.Random(1)
    def mutate(sequence, rate):
        return "".join([ rng.choice("ACGT".replace(x, "")) if rng.random() < rate else x for x in sequence ])

    reference = { chrom: "".join(rng.choices("ACGT", k=60000)) for chrom in ["chr1", "chr2"] }
    hapA = { chrom: mutate(sequence, 0.005) for chrom, sequence in reference.items() }
    hapB = { chrom: mutate(sequence, 0.005) for chrom, sequence in reference.items() }
    contigs = {
        "A_chr1": hapA["chr1"],
        "A_chr2": hapA["chr2"],
        "B_chr1_left": hapB["chr1"][:35000],
        "B_chr1_right": hapB["chr1"][35000:],
        "B_chr2": hapB["chr2"]
    }

    referenceFile = os.path.join(tmpDir, "reference.fasta")
    with open(referenceFile, "w") as fileOut:
        for chrom, sequence in reference.items():
            fileOut.write(f">{chrom}\n{sequence}\n")
    contigsFile = os.path.join(tmpDir, "contigs.fasta")
    with open(contigsFile, "w") as fileOut:
        for contig, sequence in contigs.items():
            fileOut.write(f">{contig}\n{sequence}\n")
    return referenceFile, contigsFile

def parse_phasing_rows(phasingFile):
    '''
    Returns:
        phasing -- a dictionary where keys are contig IDs and values are dictionaries
                   of that contig's column values
    '''
    with open(phasingFile, "r") as fileIn:
        header = fileIn.readline().lstrip("#").rstrip("\n").split("\t")
        return { sl[0]: dict(zip(header, sl)) for sl in [ line.rstrip("\n").split("\t") for line in fileIn ] }

def gaps_are_codon_aligned(sequence):
    '''
    Returns:
        isCodonAligned -- a boolean indicating whether every run of gaps in the sequence
                          starts at a codon boundary and has a length divisible by 3
    '''
    return all([
        match.start() % 3 == 0 and len(match.group()) % 3 == 0
        for match in re.finditer("-+", sequence)
    ])

# Define unit tests
class TestCodonTranslation(unittest.TestCase):
    def test_codons_to_proteins_basic(self):
        "Test that nucleotides are translated in the +ve strand and frame 0"
        # Arrange
        records = [Sequence("seq1", "ATGGCTAAACTG")]

        # Act
        proteins = MAFFT.codons_to_proteins(records)

        # Assert
        self.assertEqual(str(proteins[0]), "MAKL")
        self.assertEqual(proteins[0].description, "seq1", "Translation should retain the description")

    def test_codons_to_proteins_stop_codons(self):
        "Test that stop codons are translated as X so that MAFFT retains them"
        # Arrange
        records = [Sequence("seq1", "ATGTAAGCTTGATAG")]

        # Act
        proteins = MAFFT.codons_to_proteins(records)

        # Assert
        self.assertEqual(str(proteins[0]), "MXAXX")

    def test_codons_to_proteins_incomplete_codon(self):
        "Test that an incomplete final codon is not translated"
        # Arrange
        "'GC' would be translated as 'A' since truncated codons are inferrable"
        records = [Sequence("seq1", "ATGGCTGC"), Sequence("seq2", "ATGGCTG")]

        # Act
        proteins = MAFFT.codons_to_proteins(records)

        # Assert
        self.assertEqual(str(proteins[0]), "MA")
        self.assertEqual(str(proteins[1]), "MA")

    def test_codons_to_proteins_lowercase(self):
        "Test that softmasked (lowercase) nucleotides are translated"
        # Arrange
        records = [Sequence("seq1", "atgGCTaaa")]

        # Act
        proteins = MAFFT.codons_to_proteins(records)

        # Assert
        self.assertEqual(str(proteins[0]), "MAK")

class TestCodonMapping(unittest.TestCase):
    def test_proteins_to_codons_no_gaps(self):
        "Test that an ungapped protein alignment maps back to the original nucleotides"
        # Arrange
        records = [Sequence("seq1", "ATGGCTAAA"), Sequence("seq2", "ATGGCCAAG")]
        alignedProteins = [Sequence("seq1", "MAK"), Sequence("seq2", "MAK")]

        # Act
        alignedSeqs = MAFFT.proteins_to_codons(records, alignedProteins)

        # Assert
        self.assertEqual([ str(x) for x in alignedSeqs ], ["ATGGCTAAA", "ATGGCCAAG"])

    def test_proteins_to_codons_gaps(self):
        "Test that each protein gap becomes a codon gap at the same position"
        # Arrange
        records = [Sequence("seq1", "ATGGCTAAACTG"), Sequence("seq2", "ATGAAACTG")]
        alignedProteins = [Sequence("seq1", "MAKL"), Sequence("seq2", "M-KL")]

        # Act
        alignedSeqs = MAFFT.proteins_to_codons(records, alignedProteins)

        # Assert
        self.assertEqual(str(alignedSeqs[0]), "ATGGCTAAACTG")
        self.assertEqual(str(alignedSeqs[1]), "ATG---AAACTG")

    def test_proteins_to_codons_stop_codons(self):
        "Test that residues translated as X from stop codons map back to their codons"
        # Arrange
        records = [Sequence("seq1", "ATGTGAGCT"), Sequence("seq2", "ATGGCT")]
        alignedProteins = [Sequence("seq1", "MXA"), Sequence("seq2", "M-A")]

        # Act
        alignedSeqs = MAFFT.proteins_to_codons(records, alignedProteins)

        # Assert
        self.assertEqual(str(alignedSeqs[0]), "ATGTGAGCT")
        self.assertEqual(str(alignedSeqs[1]), "ATG---GCT")

    def test_proteins_to_codons_leftover_nucleotides(self):
        "Test that leftover nucleotides follow the last codon and other rows are padded"
        # Arrange
        records = [Sequence("seq1", "ATGGCTAAAGC"), Sequence("seq2", "ATGAAA")]
        alignedProteins = [Sequence("seq1", "MAK"), Sequence("seq2", "M-K")]

        # Act
        alignedSeqs = MAFFT.proteins_to_codons(records, alignedProteins)

        # Assert
        self.assertEqual(str(alignedSeqs[0]), "ATGGCTAAAGC")
        self.assertEqual(str(alignedSeqs[1]), "ATG---AAA--")

    def test_proteins_to_codons_leftover_after_trailing_gap(self):
        "Test that leftover nucleotides are placed before trailing gaps rather than after"
        # Arrange
        records = [Sequence("seq1", "ATGGCTAAA"), Sequence("seq2", "ATGGCTG")]
        alignedProteins = [Sequence("seq1", "MAK"), Sequence("seq2", "MA-")]

        # Act
        alignedSeqs = MAFFT.proteins_to_codons(records, alignedProteins)

        # Assert
        self.assertEqual(str(alignedSeqs[0]), "ATGGCTAAA")
        self.assertEqual(str(alignedSeqs[1]), "ATGGCTG--")

    def test_proteins_to_codons_preserves_case_and_description(self):
        "Test that mapping retains the original nucleotide case and full description"
        # Arrange
        records = [Sequence("seq1 a description", "atgGCTaaa")]
        alignedProteins = [Sequence("seq1 a description", "MAK")]

        # Act
        alignedSeqs = MAFFT.proteins_to_codons(records, alignedProteins)

        # Assert
        self.assertEqual(str(alignedSeqs[0]), "atgGCTaaa")
        self.assertEqual(alignedSeqs[0].description, "seq1 a description")

    def test_proteins_to_codons_order_mismatch(self):
        "Test that an error is raised if aligned proteins are not in the input order"
        # Arrange
        records = [Sequence("seq1", "ATGGCT"), Sequence("seq2", "ATGGCT")]
        alignedProteins = [Sequence("seq2", "MA"), Sequence("seq1", "MA")]

        # Act & Assert
        with self.assertRaises(ValueError):
            MAFFT.proteins_to_codons(records, alignedProteins)

    def test_translation_round_trip(self):
        "Test that translating then mapping an ungapped alignment returns the input"
        # Arrange
        records = parse_fasta(codonIndelsFile)

        # Act
        proteins = MAFFT.codons_to_proteins(records)
        alignedSeqs = MAFFT.proteins_to_codons(records, proteins)

        # Assert
        "Sequences are unaligned, so shorter ones will be padded with gaps at the end"
        for record, alignedSeq in zip(records, alignedSeqs):
            self.assertEqual(str(alignedSeq).rstrip("-"), str(record))

class TestMAFFT(unittest.TestCase):
    def test_fftnsi_requires_maxiterate(self):
        "Test that the fftnsi algorithm refuses a maxiterate < 2"
        # Arrange
        "Any existing file will do as the executable since MAFFT is not run"
        aligner = MAFFT(codonIndelsFile, algorithm="fftnsi", maxiterate=1)

        # Act & Assert
        with self.assertRaises(ValueError):
            aligner._align(codonIndelsFile)

    def test_invalid_settings(self):
        "Test that invalid MAFFT settings are refused"
        # Act & Assert
        with self.assertRaises(FileNotFoundError):
            MAFFT(os.path.join(dataDir, "does_not_exist"))
        with self.assertRaises(ValueError):
            MAFFT(codonIndelsFile, algorithm="not_an_algorithm")
        with self.assertRaises(ValueError):
            MAFFT(codonIndelsFile, thread=0)
        with self.assertRaises(ValueError):
            MAFFT(codonIndelsFile, maxiterate=-1)

    @unittest.skipIf(MAFFT_EXE is None, "mafft is not discoverable in the system PATH")
    def test_align_nucleotides(self):
        "Test plain alignment retains sequences, their order, and their case"
        # Arrange
        records = parse_fasta(codonIndelsFile)
        records[0].sequence = records[0].sequence.lower()
        aligner = MAFFT(MAFFT_EXE)

        # Act
        alignedSeqs = aligner.align(records)

        # Assert
        self.assertEqual([ x.description for x in alignedSeqs ], [ x.description for x in records ])
        self.assertEqual(len(set([ len(x) for x in alignedSeqs ])), 1, "All aligned sequences should be the same length")
        for record, alignedSeq in zip(records, alignedSeqs):
            self.assertEqual(str(alignedSeq).replace("-", ""), str(record))

    @unittest.skipIf(MAFFT_EXE is None, "mafft is not discoverable in the system PATH")
    def test_align_as_codons(self):
        "Test codon alignment places all gaps at codon boundaries in multiples of 3"
        # Arrange
        records = parse_fasta(codonIndelsFile)
        aligner = MAFFT(MAFFT_EXE)

        # Act
        alignedSeqs = aligner.align(records, asCodons=True)

        # Assert
        self.assertEqual([ x.description for x in alignedSeqs ], [ x.description for x in records ])
        self.assertEqual(len(set([ len(x) for x in alignedSeqs ])), 1, "All aligned sequences should be the same length")
        for record, alignedSeq in zip(records, alignedSeqs):
            self.assertEqual(str(alignedSeq).replace("-", ""), str(record))
            self.assertTrue(gaps_are_codon_aligned(str(alignedSeq)),
                            f"Gaps in '{alignedSeq.id}' are not codon-aligned: {str(alignedSeq)}")

    @unittest.skipIf(MAFFT_EXE is None, "mafft is not discoverable in the system PATH")
    def test_align_as_codons_with_internal_stop(self):
        "Test codon alignment of a sequence with an internal stop codon and incomplete final codon"
        # Arrange
        records = [
            Sequence("seq1", "ATGGCTAAACTGGTTGAATAA"),
            Sequence("seq2", "ATGGCTCTGGTTGAATAA"),
            Sequence("seq3", "ATGGCTAAATGAGTTGAATAAGC")
        ]
        aligner = MAFFT(MAFFT_EXE)

        # Act
        alignedSeqs = aligner.align(records, asCodons=True)

        # Assert
        self.assertEqual(str(alignedSeqs[0]), "ATGGCTAAACTGGTTGAATAA--")
        self.assertEqual(str(alignedSeqs[1]), "ATGGCT---CTGGTTGAATAA--")
        self.assertEqual(str(alignedSeqs[2]), "ATGGCTAAATGAGTTGAATAAGC")

class TestFastaToMSA(unittest.TestCase):
    @unittest.skipIf(MAFFT_EXE is None, "mafft is not discoverable in the system PATH")
    def test_fasta_to_msa_codons(self):
        "Test the 'fasta to msa' command produces a codon-aligned MSA file"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            outputFile = os.path.join(tmpDir, "msa.fasta")
            cmd = [
                sys.executable, os.path.join(baseDir, "annotarium.py"), "fasta", "to", "msa",
                "-i", codonIndelsFile,
                "-o", outputFile,
                "--codons"
            ]

            # Act
            returncode, stdout, stderr = run_subprocess(cmd)

            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            alignedSeqs = parse_fasta(outputFile)
            records = parse_fasta(codonIndelsFile)
            self.assertEqual(len(alignedSeqs), len(records))
            for alignedSeq in alignedSeqs:
                self.assertTrue(gaps_are_codon_aligned(str(alignedSeq)))

    def test_fasta_to_msa_fftnsi_validation(self):
        "Test the 'fasta to msa' command refuses fftnsi without sufficient --maxiterate"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            outputFile = os.path.join(tmpDir, "msa.fasta")
            expectedError = "--maxiterate must be >= 2"
            cmd = [
                sys.executable, os.path.join(baseDir, "annotarium.py"), "fasta", "to", "msa",
                "-i", codonIndelsFile,
                "-o", outputFile,
                "--alg", "fftnsi",
                "--mafft", MAFFT_EXE if MAFFT_EXE is not None else codonIndelsFile
            ]

            # Act
            returncode, stdout, stderr = run_subprocess(cmd)

            # Assert
            self.assertNotEqual(returncode, 0)
            self.assertTrue(expectedError in stderr, f"Expected '{expectedError}' in stderr but got: {stderr}")
            self.assertFalse(os.path.exists(outputFile), "No output file should be written")

class TestAnnotableHits(unittest.TestCase):
    def test_split_hits(self):
        "Test splitting of multi-hit values into individual hits"
        # Act & Assert
        self.assertEqual(Annotable.split_hits(".", 0), [])
        self.assertEqual(Annotable.split_hits("A1", 1), ["A1"])
        self.assertEqual(Annotable.split_hits("A1 [A2][A3]", 3), ["A1", "A2", "A3"])
        self.assertEqual(Annotable.split_hits("Zea mays [Oryza sativa]", 2), ["Zea mays", "Oryza sativa"])
    
    def test_split_hits_brackets_in_first_hit(self):
        "Test that brackets within the first hit's value do not disrupt splitting"
        # Act
        hits = Annotable.split_hits("protein kinase [partial] [kinase-like][hypothetical protein]", 3)
        
        # Assert
        self.assertEqual(hits, ["protein kinase [partial]", "kinase-like", "hypothetical protein"])
    
    def test_split_hits_count_mismatch(self):
        "Test that values which cannot be split into the expected number of hits give None"
        # Act & Assert
        self.assertIsNone(Annotable.split_hits("only one name", 2))
        self.assertIsNone(Annotable.split_hits("A1 [A2]", 3))
        self.assertIsNone(Annotable.split_hits("A1", 0))
    
    def test_join_hits(self):
        "Test joining of individual hits into the multi-hit format"
        # Act & Assert
        self.assertEqual(Annotable.join_hits([]), ".")
        self.assertEqual(Annotable.join_hits(["A1"]), "A1")
        self.assertEqual(Annotable.join_hits(["A1", "A2", "A3"]), "A1 [A2][A3]")
    
    def test_split_join_round_trip(self):
        "Test that splitting then joining every hit column returns the original values"
        # Arrange
        annotable = Annotable(annotableFile)
        
        # Act & Assert
        for rowID, row in annotable.df.iterrows():
            accessions = row["Target_accession"]
            numHits = 0 if accessions == "." else len(accessions.replace(" ", "").replace("]", "").split("["))
            for col in Annotable.HIT_COLUMNS:
                if col in annotable.df.columns:
                    hits = Annotable.split_hits(row[col], numHits)
                    self.assertEqual(Annotable.join_hits(hits), row[col], f"'{rowID}' column '{col}'")
    
    def test_load_retains_values(self):
        "Test that loading then writing an annotable retains values exactly as written"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            outputFile = os.path.join(tmpDir, "annotable.tsv")
            
            # Act
            Annotable(annotableFile).write(outputFile)
            
            # Assert
            with open(annotableFile) as originalIn, open(outputFile) as writtenIn:
                self.assertEqual(originalIn.read(), writtenIn.read())
    
    def test_filter_hits_maxhits(self):
        "Test that --maxhits retains the best hits across all hit columns"
        # Arrange
        annotable = Annotable(annotableFile)
        
        # Act
        numModified = annotable.filter_hits(maxhits=1)
        
        # Assert
        self.assertEqual(numModified, 2, "gene1 and gene2 have more than 1 hit")
        self.assertEqual(annotable["gene1"]["Target_accession"], "A1")
        self.assertEqual(annotable["gene1"]["Gene_names"], "protein kinase [partial]")
        self.assertEqual(annotable["gene1"]["Expect_value"], "1.00e-100")
        self.assertEqual(annotable["gene4"]["Target_accession"], "C1", "Single hit rows should be unchanged")
    
    def test_filter_hits_idmapping_reset(self):
        "Test that ID mapping columns are reset if their best mapped hit is removed"
        # Arrange
        annotable = Annotable(annotableFile)
        
        # Act
        annotable.filter_hits(maxhits=1)
        
        # Assert
        "gene1's best mapped hit (A1) is retained, but gene2's (B2) is removed"
        self.assertEqual(annotable["gene1"]["Best_mapped_GOs"], "GO:1")
        for col in Annotable.IDMAPPING_COLUMNS:
            self.assertEqual(annotable["gene2"][col], ".")
    
    def test_filter_hits_evalue(self):
        "Test that --evalue removes worse hits and empties rows without any passing hits"
        # Arrange
        annotable = Annotable(annotableFile)
        
        # Act
        numModified = annotable.filter_hits(evalue=1e-20)
        
        # Assert
        self.assertEqual(numModified, 2)
        self.assertEqual(annotable["gene1"]["Target_accession"], "A1 [A2]")
        self.assertEqual(annotable["gene1"]["Expect_value"], "1.00e-100 [2e-30]")
        self.assertEqual(len(annotable), 4, "Rows should not be removed by --evalue")
        for col in ["Source", *Annotable.HIT_COLUMNS, *Annotable.IDMAPPING_COLUMNS]:
            if col in annotable.df.columns: # test file is a slim table lacking some hit columns
                self.assertEqual(annotable["gene2"][col], ".", f"gene2 column '{col}' should be emptied")
    
    def test_filter_hits_evalue_then_maxhits(self):
        "Test that --maxhits is applied to the hits remaining after --evalue"
        # Arrange
        annotable = Annotable(annotableFile)
        
        # Act
        annotable.filter_hits(evalue=1e-50, maxhits=2)
        
        # Assert
        self.assertEqual(annotable["gene1"]["Target_accession"], "A1")
    
    def test_filter_hits_inconsistent_columns(self):
        "Test that an error is raised if a hit column disagrees with Target_accession"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            badFile = os.path.join(tmpDir, "bad.tsv")
            with open(badFile, "w") as fileOut:
                fileOut.write("#Query\tTarget_accession\tGene_names\nq1\tA1 [A2]\tonly one name\n")
            annotable = Annotable(badFile)
            
            # Act & Assert
            with self.assertRaises(ValueError):
                annotable.filter_hits(maxhits=1)

class TestAnnotableFilter(unittest.TestCase):
    def test_include_columns(self):
        "Test that --include only matches within the specified --columns"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            outputFile = os.path.join(tmpDir, "filtered.tsv")
            
            # Act
            returncode, stdout, stderr = run_annotable_filter(outputFile, "--include", "Oryza",
                                                              "--columns", "Taxa_names")
            
            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            self.assertEqual(Annotable(outputFile).rowNames, ["gene1", "gene4"])
    
    def test_include_all_columns(self):
        "Test that --include matches within any column by default"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            outputFile = os.path.join(tmpDir, "filtered.tsv")
            
            # Act
            returncode, stdout, stderr = run_annotable_filter(outputFile, "--include", "GO:3", "transporter")
            
            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            self.assertEqual(Annotable(outputFile).rowNames, ["gene2", "gene4"])
    
    def test_exclude_query_column(self):
        "Test that --exclude can match the #Query column"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            outputFile = os.path.join(tmpDir, "filtered.tsv")
            
            # Act
            returncode, stdout, stderr = run_annotable_filter(outputFile, "--exclude", "gene1", "gene3",
                                                              "--columns", "#Query")
            
            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            self.assertEqual(Annotable(outputFile).rowNames, ["gene2", "gene4"])
    
    def test_include_after_maxhits(self):
        "Test that --include only matches against hits retained by --maxhits"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            outputFile = os.path.join(tmpDir, "filtered.tsv")
            
            # Act
            "'transporter' is only gene2's second hit"
            returncode, stdout, stderr = run_annotable_filter(outputFile, "--include", "transporter",
                                                              "--maxhits", "1")
            
            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            self.assertEqual(Annotable(outputFile).rowNames, [])
    
    def test_invalid_column(self):
        "Test that --columns values not in the annotable header are refused"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            outputFile = os.path.join(tmpDir, "filtered.tsv")
            expectedError = "--columns value 'Nope' is not a header"
            
            # Act
            returncode, stdout, stderr = run_annotable_filter(outputFile, "--include", "x",
                                                              "--columns", "Nope")
            
            # Assert
            self.assertNotEqual(returncode, 0)
            self.assertTrue(expectedError in stderr, f"Expected '{expectedError}' in stderr but got: {stderr}")
            self.assertFalse(os.path.exists(outputFile), "No output file should be written")
    
    def test_nothing_to_do(self):
        "Test that the command is refused if no filtering arguments are given"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            outputFile = os.path.join(tmpDir, "filtered.tsv")
            expectedError = "must be provided for 'annotable filter' to do anything"
            
            # Act
            returncode, stdout, stderr = run_annotable_filter(outputFile)
            
            # Assert
            self.assertNotEqual(returncode, 0)
            self.assertTrue(expectedError in stderr, f"Expected '{expectedError}' in stderr but got: {stderr}")

class TestAnnotableIds(unittest.TestCase):
    def test_ids(self):
        "Test that 'annotable ids' writes the first column without its header"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            outputFile = os.path.join(tmpDir, "ids.txt")
            cmd = [
                sys.executable, os.path.join(baseDir, "annotarium.py"), "annotable", "ids",
                "-i", annotableFile,
                "-o", outputFile
            ]
            
            # Act
            returncode, stdout, stderr = run_subprocess(cmd)
            
            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            with open(outputFile) as fileIn:
                self.assertEqual(fileIn.read(), "gene1\ngene2\ngene3\ngene4\n")

class TestPhasingIntervals(unittest.TestCase):
    def test_merge_intervals(self):
        "Test that overlapping and abutting intervals are merged"
        # Arrange
        intervals = [(50, 60), (0, 10), (5, 20), (20, 30)]

        # Act
        merged = merge_intervals(intervals)

        # Assert
        self.assertEqual(merged, [(0, 30), (50, 60)])

    def test_parse_paf_filters(self):
        "Test that PAF parsing skips secondary, low MAPQ, and short alignments"
        # Act
        alignments = parse_paf(pafFile, minQ=1, minAlignLen=10000)

        # Assert
        self.assertEqual(alignments, {
            "q1": [(0, 20000, "chr1", 1000, 21000, "+")],
            "q2": [(5000, 30000, "chr1", 55000, 80000, "+")]
        })

    def test_parse_paf_secondary(self):
        "Test that secondary alignments are skipped even without other filters"
        # Act
        alignments = parse_paf(pafFile)

        # Assert
        self.assertEqual(len(alignments["q1"]), 2)
        self.assertFalse(any([ aln[2] == "chr2" and aln[3] == 5000 for aln in alignments["q1"] ]))

    def test_deduplicate_query(self):
        "Test that a query region aligning to multiple reference regions only keeps the longest"
        # Arrange
        alignments = { "q1": [
            (0, 10000, "chr1", 0, 10000, "+"),
            (0, 8000, "chr2", 0, 8000, "+"), # entirely within the first alignment's query region
            (9000, 20000, "chr3", 0, 11000, "+") # mostly outside the first alignment's query region
        ]}

        # Act
        dedupAlignments = deduplicate_query(alignments)

        # Assert
        self.assertEqual(sorted([ aln[2] for aln in dedupAlignments["q1"] ]), ["chr1", "chr3"])

    def test_build_footprints(self):
        "Test that footprints are merged per reference contig"
        # Arrange
        alignments = { "q1": [
            (0, 1000, "chr1", 0, 1000, "+"),
            (1000, 2000, "chr1", 900, 1900, "+"),
            (2000, 3000, "chr2", 500, 1500, "-")
        ]}

        # Act
        footprints = build_footprints(alignments)

        # Assert
        self.assertEqual(footprints, { "q1": { "chr1": [(0, 1900)], "chr2": [(500, 1500)] } })

    def test_infer_orientations(self):
        "Test that orientation per reference contig follows the strand with the most aligned query bases"
        # Arrange
        alignments = { "q1": [
            (0, 1000, "chr1", 0, 1000, "-"),
            (1000, 1500, "chr1", 1000, 1500, "+"), # minority inversion
            (1500, 3000, "chr2", 0, 1500, "+"),
            (3000, 3500, "chr3", 0, 500, "+"),
            (3500, 4000, "chr3", 500, 1000, "-") # tie
        ]}

        # Act
        orientations = infer_orientations(alignments)

        # Assert
        self.assertEqual(orientations, { "q1": { "chr1": "-", "chr2": "+", "chr3": "+" } })

class TestHaplotypePartitioner(unittest.TestCase):
    def test_diploid_partition(self):
        "Test that overlapping contigs are split into separate partitions"
        # Arrange
        footprints = {
            "A": { "chr1": [(0, 1000)] },
            "B1": { "chr1": [(0, 600)] },
            "B2": { "chr1": [(600, 1000)] }
        }

        # Act
        states = HaplotypePartitioner(2).partition(footprints)

        # Assert
        self.assertEqual(states["B1"], states["B2"])
        self.assertNotEqual(states["A"], states["B1"])
        self.assertEqual(states["A"], 0, "The partition with most coverage should be the first")

    def test_redundant_unassigned(self):
        "Test that a contig which only adds redundancy is left unassigned"
        # Arrange
        footprints = {
            "A": { "chr1": [(0, 1000)] },
            "B": { "chr1": [(0, 1000)] },
            "dup": { "chr1": [(100, 300)] }
        }

        # Act
        states = HaplotypePartitioner(2).partition(footprints)

        # Assert
        self.assertEqual(sorted([states["A"], states["B"]]), [0, 1])
        self.assertIsNone(states["dup"])

    def test_haploid_partition(self):
        "Test that ploidy 1 keeps the contigs giving the most non-redundant coverage"
        # Arrange
        footprints = {
            "A": { "chr1": [(0, 1000)] },
            "B1": { "chr1": [(0, 300)] },
            "B2": { "chr1": [(900, 1500)] }
        }

        # Act
        states = HaplotypePartitioner(1).partition(footprints)

        # Assert
        self.assertEqual(states, { "A": 0, "B1": None, "B2": 0 })

    def test_penalty(self):
        "Test that the penalty controls whether a partially redundant contig is assigned"
        # Arrange
        footprints = {
            "A": { "chr1": [(0, 1000)] },
            "B": { "chr1": [(800, 1100)] } # 100 bp of new coverage and 200 bp of redundancy
        }

        # Act
        lenientStates = HaplotypePartitioner(1, penalty=0.25).partition(footprints)
        strictStates = HaplotypePartitioner(1, penalty=1.0).partition(footprints)

        # Assert
        self.assertEqual(lenientStates["B"], 0)
        self.assertIsNone(strictStates["B"])

    def test_translocated_contig(self):
        "Test that a contig contributing coverage to two reference contigs is phased with both"
        # Arrange
        footprints = {
            "A_chr1": { "chr1": [(0, 1000)] },
            "A_chr2": { "chr2": [(0, 1000)] },
            "B_chr1": { "chr1": [(0, 700)] },
            "B_trans": { "chr1": [(700, 1000)], "chr2": [(0, 400)] },
            "B_chr2": { "chr2": [(400, 1000)] }
        }

        # Act
        states = HaplotypePartitioner(2).partition(footprints)

        # Assert
        self.assertEqual(len(set([states["B_chr1"], states["B_trans"], states["B_chr2"]])), 1)
        self.assertNotEqual(states["B_trans"], states["A_chr1"])
        self.assertNotEqual(states["B_trans"], states["A_chr2"])

    def test_heuristic_matches_exact(self):
        "Test that the heuristic search finds the optimal score for random small components"
        # Arrange
        rng = random.Random(0)
        for _ in range(50):
            footprints = {}
            for i in range(rng.randint(3, 7)):
                start = rng.randint(0, 900)
                footprints[f"c{i}"] = { rng.choice(["chr1", "chr2"]): [(start, start + rng.randint(20, 400))] }
            partitioner = HaplotypePartitioner(rng.randint(1, 3), penalty=rng.choice([0.5, 1.0, 3.0]))
            neighbours = partitioner.find_neighbours(footprints)

            for component in partitioner.find_components(neighbours):
                # Act
                exactScore = partitioner.score(partitioner._exact(component, footprints), footprints)
                heuristicScore = partitioner.score(partitioner._heuristic(component, footprints, neighbours), footprints)

                # Assert
                self.assertAlmostEqual(heuristicScore, exactScore)

    def test_unique_and_group_coverage(self):
        "Test the calculation of unique and per-partition coverage"
        # Arrange
        footprints = {
            "A1": { "chr1": [(0, 600)] },
            "A2": { "chr1": [(500, 1000)] },
            "B": { "chr1": [(0, 1000)], "chr2": [(0, 200)] }
        }
        states = { "A1": 0, "A2": 0, "B": 1 }

        # Act
        uniqueBp = unique_coverage(states, footprints)
        coverage = group_coverage(states, footprints, 2)

        # Assert
        self.assertEqual(uniqueBp, { "A1": 500, "A2": 400, "B": 1200 })
        self.assertEqual(coverage, { "chr1": [1000, 1000], "chr2": [0, 200] })

class TestAssemblyPhase(unittest.TestCase):
    @unittest.skipIf(MINIMAP2_EXE is None, "minimap2 is not discoverable in the system PATH")
    def test_assembly_phase(self):
        "Test the 'assembly phase' command separates haplotypes and reports coverage"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            referenceFile, contigsFile = write_phasing_inputs(tmpDir)
            outputDir = os.path.join(tmpDir, "phase")
            cmd = [
                sys.executable, os.path.join(baseDir, "annotarium.py"), "assembly", "phase",
                "-i", contigsFile,
                "-r", referenceFile,
                "-o", outputDir,
                "--ploidy", "2",
                "--minAlignLen", "1000"
            ]

            # Act
            returncode, stdout, stderr = run_subprocess(cmd)

            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            phasing = parse_phasing_rows(os.path.join(outputDir, "phasing.tsv"))
            self.assertTrue(all([ x["status"] == "assigned" for x in phasing.values() ]))
            self.assertEqual(phasing["B_chr1_left"]["partition"], phasing["B_chr1_right"]["partition"])
            self.assertNotEqual(phasing["A_chr1"]["partition"], phasing["B_chr1_left"]["partition"])
            self.assertNotEqual(phasing["A_chr2"]["partition"], phasing["B_chr2"]["partition"])

            with open(os.path.join(outputDir, "coverage.tsv"), "r") as fileIn:
                fileIn.readline()
                for line in fileIn:
                    sl = line.rstrip("\n").split("\t")
                    self.assertTrue(float(sl[3]) > 95 and float(sl[5]) > 95, f"Expected near-complete coverage but got: {line}")

    @unittest.skipIf(MINIMAP2_EXE is None, "minimap2 is not discoverable in the system PATH")
    def test_assembly_phase_no_overwrite(self):
        "Test the 'assembly phase' command refuses to overwrite a previous result"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            referenceFile, contigsFile = write_phasing_inputs(tmpDir)
            outputDir = os.path.join(tmpDir, "phase")
            expectedError = "already exists"
            cmd = [
                sys.executable, os.path.join(baseDir, "annotarium.py"), "assembly", "phase",
                "-i", contigsFile,
                "-r", referenceFile,
                "-o", outputDir,
                "--ploidy", "2",
                "--minAlignLen", "1000"
            ]
            run_subprocess(cmd)

            # Act
            returncode, stdout, stderr = run_subprocess(cmd)

            # Assert
            self.assertNotEqual(returncode, 0)
            self.assertTrue(expectedError in stderr, f"Expected '{expectedError}' in stderr but got: {stderr}")

class TestRagTag(unittest.TestCase):
    def test_parse_phasing_tsv(self):
        "Test that only assigned contigs are parsed, grouped by partition in numeric order"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            phasingFile = os.path.join(tmpDir, "phasing.tsv")
            with open(phasingFile, "w") as fileOut:
                fileOut.write("#contig\tlength\tpartition\tstatus\taligned_bp\tunique_bp\treference_regions\treference_bp\n")
                fileOut.write("c1\t100\t10\tassigned\t100\t100\tchr1:1-100\t100\n")
                fileOut.write("c2\t100\t2\tassigned\t100\t100\tchr1:1-100\t100\n")
                fileOut.write("c3\t100\t.\tredundant\t100\t0\tchr1:1-100\t100\n")
                fileOut.write("c4\t100\t.\tunaligned\t0\t0\t.\t.\n")
                fileOut.write("c5\t100\t2\tassigned\t100\t100\tchr2:1-100\t100\n")

            # Act
            partitions = parse_phasing_tsv(phasingFile)

            # Assert
            self.assertEqual(partitions, { "2": ["c2", "c5"], "10": ["c1"] })
            self.assertEqual(list(partitions), ["2", "10"])

    def test_parse_phasing_tsv_header(self):
        "Test that a file without the 'assembly phase' header is refused"
        # Act & Assert
        with self.assertRaises(ValueError):
            parse_phasing_tsv(annotableFile)

    def test_parse_agp(self):
        "Test that AGP parsing groups contigs by scaffold, skipping gaps and recognising unplaced contigs"
        # Act
        scaffolds, placed = RagTag.parse_agp(agpFile)

        # Assert
        self.assertEqual(scaffolds, { "chr1": ["c1", "c2"], "chr2": ["c3"], "c4": ["c4"] })
        self.assertEqual(placed, {"chr1", "chr2"})

    def test_rename_scaffolds(self):
        "Test that the RagTag suffix is removed from scaffolds but not unplaced contigs"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            ragtagFile = os.path.join(tmpDir, "ragtag.scaffold.fasta")
            outputFile = os.path.join(tmpDir, "renamed.fasta")
            with open(ragtagFile, "w") as fileOut:
                fileOut.write(">chr1_RagTag\nACGTACGT\n>c4\nTTTT\n")

            # Act
            RagTag.rename_scaffolds(ragtagFile, outputFile, lineLength=3)

            # Assert
            with open(outputFile, "r") as fileIn:
                self.assertEqual(fileIn.read(), ">chr1\nACG\nTAC\nGT\n>c4\nTTT\nT\n")

    def test_rename_scaffolds_collision(self):
        "Test that renaming is refused if it would cause duplicate IDs, without writing output"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            ragtagFile = os.path.join(tmpDir, "ragtag.scaffold.fasta")
            outputFile = os.path.join(tmpDir, "renamed.fasta")
            with open(ragtagFile, "w") as fileOut:
                fileOut.write(">chr1_RagTag\nACGT\n>chr1\nTTTT\n")

            # Act & Assert
            with self.assertRaises(ValueError):
                RagTag.rename_scaffolds(ragtagFile, outputFile)
            self.assertFalse(os.path.exists(outputFile), "No output file should be written")

class TestAssemblyScaffold(unittest.TestCase):
    @unittest.skipIf(MINIMAP2_EXE is None or RAGTAG_EXE is None,
                     "minimap2 and/or ragtag.py are not discoverable in the system PATH")
    def test_assembly_scaffold(self):
        "Test the 'assembly scaffold' command produces one scaffolded FASTA per partition"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            referenceFile, contigsFile = write_phasing_inputs(tmpDir)
            phaseDir = os.path.join(tmpDir, "phase")
            scaffoldDir = os.path.join(tmpDir, "scaffold")
            phaseCmd = [
                sys.executable, os.path.join(baseDir, "annotarium.py"), "assembly", "phase",
                "-i", contigsFile,
                "-r", referenceFile,
                "-o", phaseDir,
                "--ploidy", "2",
                "--minAlignLen", "1000"
            ]
            scaffoldCmd = [
                sys.executable, os.path.join(baseDir, "annotarium.py"), "assembly", "scaffold",
                "-i", os.path.join(phaseDir, "phasing.tsv"),
                "-f", contigsFile,
                "-r", referenceFile,
                "-o", scaffoldDir
            ]
            returncode, stdout, stderr = run_subprocess(phaseCmd)
            self.assertEqual(returncode, 0, f"Expected no error from 'assembly phase' but got: {stderr}")

            # Act
            returncode, stdout, stderr = run_subprocess(scaffoldCmd)

            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            for partition in ["1", "2"]:
                records = parse_fasta(os.path.join(scaffoldDir, f"partition{partition}.fasta"))
                self.assertEqual(sorted([ record.id for record in records ]), ["chr1", "chr2"])
                self.assertTrue(all([ len(record) >= 59000 for record in records ]))

class TestHomologyPair(unittest.TestCase):
    def write_phasing_tsv(self, tmpDir):
        '''
        Returns:
            phasingFile -- a string indicating the location of a phasing TSV written to tmpDir
        '''
        phasingFile = os.path.join(tmpDir, "phasing.tsv")
        with open(phasingFile, "w") as fileOut:
            fileOut.write("#contig\tlength\tpartition\tstatus\taligned_bp\tunique_bp\treference_regions\treference_bp\treference_strand\n")
            fileOut.write("scafA\t1000\t1\tassigned\t1000\t1000\tchr1:1-900,chr2:5-104\t900,100\t+,-\n")
            fileOut.write("scafB\t500\t1\tassigned\t500\t500\tchr2:1-500\t500\t-\n")
            fileOut.write("frag\t300\t2\tassigned\t300\t300\tchr1:1-300\t300\t+\n")
            fileOut.write("dup\t300\t.\tredundant\t300\t0\tchr3:1-300\t300\t+\n")
            fileOut.write("junk\t100\t.\tunaligned\t0\t0\t.\t.\t.\n")
        return phasingFile

    def test_parse_phasing_matches(self):
        "Test that reference regions of assigned contigs are parsed alongside their aligned bases and strand"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            phasingFile = self.write_phasing_tsv(tmpDir)

            # Act
            matches = parse_phasing_matches(phasingFile)
            partitionMatches = parse_phasing_matches(phasingFile, partition="2")

            # Assert
            self.assertEqual(list(matches), ["scafA", "scafB", "frag"])
            self.assertEqual(matches["scafA"], [
                { "reference": "chr1", "start": 1, "end": 900, "bp": 900, "strand": "+" },
                { "reference": "chr2", "start": 5, "end": 104, "bp": 100, "strand": "-" }
            ])
            self.assertEqual(list(partitionMatches), ["frag"])

    def test_parse_phasing_matches_old_format(self):
        "Test that a phasing TSV lacking the reference_strand column is refused"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            phasingFile = os.path.join(tmpDir, "phasing.tsv")
            with open(phasingFile, "w") as fileOut:
                fileOut.write("#contig\tlength\tpartition\tstatus\taligned_bp\tunique_bp\treference_regions\treference_bp\n")
                fileOut.write("scafA\t1000\t1\tassigned\t1000\t1000\tchr1:1-1000\t1000\n")

            # Act & Assert
            with self.assertRaises(ValueError):
                parse_phasing_matches(phasingFile)

    def test_pair_homologous_contigs(self):
        "Test that each reference contig is paired with the query contig best matching it"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            matches = parse_phasing_matches(self.write_phasing_tsv(tmpDir))

            # Act
            pairs, unpaired = pair_homologous_contigs(matches)

            # Assert
            self.assertEqual({ reference: (pair["contig"], pair["strand"]) for reference, pair in pairs.items() },
                             { "chr1": ("scafA", "+"), "chr2": ("scafB", "-") })
            self.assertEqual(unpaired, ["frag"])

    @unittest.skipIf(MINIMAP2_EXE is None, "minimap2 is not discoverable in the system PATH")
    def test_assembly_homology_pair(self):
        "Test the 'assembly homology_pair' command renames and reverse complements contigs to match the reference"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            referenceFile, contigsFile = write_phasing_inputs(tmpDir)
            records = { record.id: record for record in parse_fasta(contigsFile) }
            queryFile = os.path.join(tmpDir, "query.fasta")
            with open(queryFile, "w") as fileOut:
                fileOut.write(records["A_chr2"].reverse_complement().format().replace(">A_chr2", ">scaffold_2"))
                fileOut.write(Sequence("scaffold_1", str(records["A_chr1"])).format())
                fileOut.write(records["B_chr1_left"].format())

            phaseDir = os.path.join(tmpDir, "phase")
            outputFile = os.path.join(tmpDir, "paired.fasta")
            orderedFile = os.path.join(tmpDir, "paired.ordered.fasta")
            phaseCmd = [
                sys.executable, os.path.join(baseDir, "annotarium.py"), "assembly", "phase",
                "-i", queryFile,
                "-r", referenceFile,
                "-o", phaseDir,
                "--ploidy", "2",
                "--minAlignLen", "1000"
            ]
            pairCmd = [
                sys.executable, os.path.join(baseDir, "annotarium.py"), "assembly", "homology_pair",
                "-i", os.path.join(phaseDir, "phasing.tsv"),
                "-f", queryFile,
                "-o", outputFile
            ]
            orderedCmd = pairCmd[:-1] + [orderedFile, "-r", referenceFile]
            returncode, stdout, stderr = run_subprocess(phaseCmd)
            self.assertEqual(returncode, 0, f"Expected no error from 'assembly phase' but got: {stderr}")

            # Act
            returncode, stdout, stderr = run_subprocess(pairCmd)
            orderedReturncode, orderedStdout, orderedStderr = run_subprocess(orderedCmd)

            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            self.assertEqual(orderedReturncode, 0, f"Expected no error with -r but got: {orderedStderr}")
            paired = [ (record.id, str(record)) for record in parse_fasta(outputFile) ]
            self.assertEqual(paired, [ ("chr2", str(records["A_chr2"])), ("chr1", str(records["A_chr1"])) ],
                             "Output should be in -f file order without -r")
            ordered = [ (record.id, str(record)) for record in parse_fasta(orderedFile) ]
            self.assertEqual(ordered, [ ("chr1", str(records["A_chr1"])), ("chr2", str(records["A_chr2"])) ],
                             "Output should be in -r file order with -r")

class TestRepeatMaskerClassify(unittest.TestCase):
    LABELS = { "CITME_RE696_RIX": "LINE/Unknown", "MURKO_RE46_RIX": "LINE/Unknown",
               "CITME_RE4_RLC_sire": "LTR/Copia" }

    def write_tsv(self, tmpDir, labels):
        '''
        Returns:
            tsvFile -- a string indicating the location of a class/family TSV written to tmpDir
        '''
        tsvFile = os.path.join(tmpDir, "classes.tsv")
        with open(tsvFile, "w") as fileOut:
            fileOut.write("#repeat_id\tclass_family\n")
            for repeatID, classFamily in labels.items():
                fileOut.write(f"{repeatID}\t{classFamily}\n")
        return tsvFile

    def run_classify(self, tmpDir, labels, extraArgs=[]):
        '''
        Returns:
            returncode, stdout, stderr -- as per run_subprocess
            outputFile -- a string indicating the location of the output .out file
        '''
        outputFile = os.path.join(tmpDir, "classified.out")
        cmd = [
            sys.executable, os.path.join(baseDir, "annotarium.py"), "repeatmasker", "classify",
            "-i", repeatmaskerFile,
            "-t", self.write_tsv(tmpDir, labels),
            "-o", outputFile
        ] + extraArgs
        return *run_subprocess(cmd), outputFile

    def test_parse_class_family_tsv(self):
        "Test that a TSV is parsed with comments skipped and conflicting labels rejected"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            tsvFile = self.write_tsv(tmpDir, self.LABELS)
            conflictFile = os.path.join(tmpDir, "conflict.tsv")
            with open(conflictFile, "w") as fileOut:
                fileOut.write("rep1\tLTR/Copia\nrep1\tLTR/Gypsy\n")

            # Act
            classFamilies = parse_class_family_tsv(tsvFile)

            # Assert
            self.assertEqual(classFamilies, self.LABELS)
            with self.assertRaises(ValueError):
                parse_class_family_tsv(conflictFile)

    def test_classify(self):
        "Test that class/family labels are replaced, other values are kept, and columns stay aligned"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            with open(repeatmaskerFile, "r") as fileIn:
                original = fileIn.read().splitlines()

            # Act
            returncode, stdout, stderr, outputFile = self.run_classify(tmpDir, self.LABELS)

            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            with open(outputFile, "r") as fileIn:
                classified = fileIn.read().splitlines()
            self.assertEqual(classified[:3], original[:3], "Header should be unchanged when labels fit")
            for before, after in zip(original[3:], classified[3:]):
                beforeFields, afterFields = before.split(), after.split()
                self.assertEqual(afterFields[10], self.LABELS[beforeFields[9]])
                self.assertEqual(afterFields[:10] + afterFields[11:], beforeFields[:10] + beforeFields[11:])
                self.assertEqual(len(after), len(before), "Column alignment should be kept")

    def test_classify_widens_column(self):
        "Test that a label longer than the class/family column widens it, header included"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            labels = { **self.LABELS, "CITME_RE4_RLC_sire": "Retrotransposon/Unknown_long_label" }

            # Act
            returncode, stdout, stderr, outputFile = self.run_classify(tmpDir, labels)

            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            with open(repeatmaskerFile, "r") as fileIn:
                original = fileIn.read().splitlines()
            with open(outputFile, "r") as fileIn:
                classified = fileIn.read().splitlines()
            def begin_offsets(lines): # distance between the end of the 'begin' header and data values
                headerEnd = lines[1].rindex("begin") + len("begin")
                return { list(re.finditer(r"\S+", line))[11].end() - headerEnd for line in lines[3:] }
            self.assertGreater(len(classified[3]), len(original[3]), "The class/family column should widen")
            self.assertEqual(len(begin_offsets(original)), 1, "Test data should have an aligned 'begin' column")
            self.assertEqual(begin_offsets(classified), begin_offsets(original),
                             "The 'begin' column should stay aligned under its header")

    def test_classify_missing(self):
        "Test that repeat IDs absent from the TSV cause an error unless --allowMissing is given"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            labels = { k:v for k,v in self.LABELS.items() if k != "MURKO_RE46_RIX" }

            # Act
            returncode, stdout, stderr, outputFile = self.run_classify(tmpDir, labels)
            self.assertNotEqual(returncode, 0, "Expected an error for a missing repeat ID")
            self.assertFalse(os.path.exists(outputFile), "No output should be written on error")
            allowReturncode, allowStdout, allowStderr, outputFile = self.run_classify(tmpDir, labels, ["--allowMissing"])

            # Assert
            self.assertIn("MURKO_RE46_RIX", stderr)
            self.assertEqual(allowReturncode, 0, f"Expected no error with --allowMissing but got: {allowStderr}")
            with open(outputFile, "r") as fileIn:
                classes = { line.split()[9]: line.split()[10] for line in fileIn.read().splitlines()[3:] }
            self.assertEqual(classes["MURKO_RE46_RIX"], "Unspecified", "Missing IDs should keep their label")
            self.assertEqual(classes["CITME_RE4_RLC_sire"], "LTR/Copia")

def write_compleasm_run(runDir, rows, lineage="test_odb12"):
    '''
    Writes a minimal compleasm output directory.
    
    Parameters:
        runDir -- a string indicating the compleasm directory to create
        rows -- a list of (geneID, status, contigID) tuples; contigID is None for missing genes
    '''
    lineageDir = os.path.join(runDir, lineage)
    os.makedirs(lineageDir)
    with open(os.path.join(lineageDir, "full_table_busco_format.tsv"), "w") as fileOut:
        fileOut.write("# compleasm version is: 0.2.10\n")
        fileOut.write(f"# The lineage dataset is: {lineage} (Creation date: 2026-05-22, number of genomes: 1, number of BUSCOs: 3)\n")
        fileOut.write("# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n")
        for geneID, status, contigID in rows:
            if contigID is None:
                fileOut.write(f"{geneID}\t{status}\n")
            else:
                fileOut.write(f"{geneID}\t{status}\t{contigID}\t1\t100\t+\t100\t100\n")

class TestCompleasmSetOptimise(unittest.TestCase):
    # A translocation leaves gene A twice in hap1 and gene B twice in hap2; swapping
    # either chromosome between sets makes both genes single-copy in both sets
    HAP1 = [ ("A", "Duplicated", "chr1"), ("A", "Duplicated", "chr2"), ("B", "Missing", None), ("C", "Complete", "ptg1") ]
    HAP2 = [ ("A", "Missing", None), ("B", "Duplicated", "chr1"), ("B", "Duplicated", "chr2"), ("C", "Missing", None) ]

    def write_inputs(self, tmpDir):
        '''
        Returns:
            runDirs -- a list of compleasm directory locations
            fastaFiles -- a list of FASTA file locations, with sequences named after their file and contig
        '''
        runDirs, fastaFiles = [], []
        for hap, rows, contigs in [ ("hap1", self.HAP1, ["chr1", "chr2", "ptg1"]), ("hap2", self.HAP2, ["chr1", "chr2"]) ]:
            runDir = os.path.join(tmpDir, f"{hap}_compleasm")
            write_compleasm_run(runDir, rows)
            fastaFile = os.path.join(tmpDir, f"{hap}.fasta")
            with open(fastaFile, "w") as fileOut:
                for contig in contigs:
                    fileOut.write(f">{contig}\n{'ACGT' if hap == 'hap1' else 'TTTT'}{contig[-1]}\n")
            runDirs.append(runDir)
            fastaFiles.append(fastaFile)
        return runDirs, fastaFiles

    def test_optimiser(self):
        "Test that swapping a translocated chromosome makes duplicated genes single-copy in each set"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            runDirs, fastaFiles = self.write_inputs(tmpDir)
            tables = [ CompleasmTable(os.path.join(runDir, "test_odb12", "full_table_busco_format.tsv")) for runDir in runDirs ]
            contigGenes = { (fileIndex, contigID): genes for fileIndex, table in enumerate(tables)
                            for contigID, genes in table.completeHits.items() }
            groups = { "chr1": [0, 1], "chr2": [0, 1], "ptg1": [0] }
            optimiser = SetOptimiser(contigGenes, groups, 2)

            # Act
            assignment = optimiser.optimise()

            # Assert
            self.assertEqual(tables[0].lineage, "test_odb12")
            self.assertEqual(tables[0].genes, {"A", "B", "C"})
            self.assertEqual(optimiser.score(optimiser.identity()), (3, 1), "Only C is single-copy before optimising")
            self.assertEqual(optimiser.score(assignment), (5, 5), "A and B should be single-copy in both sets")
            self.assertNotIn("ptg1", optimiser.groups, "Contigs in one file only should not move")
            contigSets = optimiser.contig_sets(assignment, [["chr1", "chr2", "ptg1"], ["chr1", "chr2"]])
            self.assertEqual(contigSets[(0, "ptg1")], 0)
            self.assertNotEqual(contigSets[(0, "chr1")], contigSets[(0, "chr2")])

    def test_set_optimise_to_grabber(self):
        "Test that 'qc compleasm set_optimise' output TSVs can be given to 'fasta grabber'"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            runDirs, fastaFiles = self.write_inputs(tmpDir)
            outputDir = os.path.join(tmpDir, "sets")
            optimiseCmd = [ sys.executable, os.path.join(baseDir, "annotarium.py"), "qc", "compleasm", "set_optimise",
                            "-i", *runDirs, "-f", *fastaFiles, "-o", outputDir ]
            set1Fasta = os.path.join(tmpDir, "set1.fasta")
            grabCmd = [ sys.executable, os.path.join(baseDir, "annotarium.py"), "fasta", "grabber",
                        "-i", *fastaFiles, "-s", os.path.join(outputDir, "set1.tsv"), "-o", set1Fasta ]

            # Act
            returncode, stdout, stderr = run_subprocess(optimiseCmd)
            self.assertEqual(returncode, 0, f"Expected no error from 'set_optimise' but got: {stderr}")
            grabReturncode, grabStdout, grabStderr = run_subprocess(grabCmd)

            # Assert
            pairs = parse_file_contig_tsv(os.path.join(outputDir, "set1.tsv"))
            self.assertEqual(len(pairs), 3, "set1 should hold chr1, chr2 and ptg1")
            self.assertEqual({ fastaName for fastaName, _ in pairs }, {"hap1.fasta", "hap2.fasta"},
                             "set1 should draw contigs from both files")
            self.assertEqual(grabReturncode, 0, f"Expected no error from 'fasta grabber' but got: {grabStderr}")
            grabbed = { record.id: str(record) for record in parse_fasta(set1Fasta) }
            expected = { contigID: ("ACGT" if fastaName == "hap1.fasta" else "TTTT") + contigID[-1] for fastaName, contigID in pairs }
            self.assertEqual(grabbed, expected)

class TestCompleasmTabulate(unittest.TestCase):
    GENOME = "## lineage: test_odb12\nS:90.00%, 90\nD:5.00%, 5\nF:2.00%, 2\nI:1.00%, 1\nM:2.00%, 2\nN:100\n"
    PROTEIN = "## lineage: test_odb12\nS:80.00%, 80\nD:10.00%, 10\nF:5.00%, 5\nM:5.00%, 5\nN:100\n"

    def write_summary(self, runDir, content):
        os.makedirs(runDir)
        with open(os.path.join(runDir, "summary.txt"), "w") as fileOut:
            fileOut.write("# compleasm version is: 0.2.10\n" + content)

    def run_tabulate(self, inputDirs, outputFile):
        cmd = [ sys.executable, os.path.join(baseDir, "annotarium.py"), "qc", "compleasm", "tabulate",
                "-i", *inputDirs, "-o", outputFile ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
        with open(outputFile, "r") as fileIn:
            return [ line.rstrip("\n").split("\t") for line in fileIn ]

    def test_parse_summary(self):
        "Test that each lineage of a summary file is parsed separately"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            runDir = os.path.join(tmpDir, "run")
            self.write_summary(runDir, self.GENOME + self.PROTEIN.replace("test_odb12", "other_odb12"))

            # Act
            lineageStats = parse_compleasm_summary(os.path.join(runDir, "summary.txt"))

            # Assert
            self.assertEqual(list(lineageStats), ["test_odb12", "other_odb12"])
            self.assertEqual(lineageStats["test_odb12"]["I"], "1.00")
            self.assertEqual(lineageStats["other_odb12"]["N"], "100")

    def test_tabulate(self):
        "Test that protein runs get '.' for 'I' alongside genome runs, and lack 'I' on their own"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            mixedDir, proteinDir = os.path.join(tmpDir, "mixed"), os.path.join(tmpDir, "proteins")
            self.write_summary(os.path.join(mixedDir, "genome_compleasm"), self.GENOME)
            self.write_summary(os.path.join(mixedDir, "protein.compleasm"), self.PROTEIN)
            self.write_summary(os.path.join(proteinDir, "protein"), self.PROTEIN)

            # Act
            mixed = self.run_tabulate([mixedDir], os.path.join(tmpDir, "mixed.tsv"))
            proteins = self.run_tabulate([proteinDir], os.path.join(tmpDir, "proteins.tsv"))

            # Assert
            self.assertEqual(mixed, [
                ["sample", "lineage", "S", "D", "F", "M", "I", "N"],
                ["genome", "test_odb12", "90.00", "5.00", "2.00", "2.00", "1.00", "100"],
                ["protein", "test_odb12", "80.00", "10.00", "5.00", "5.00", ".", "100"]
            ])
            self.assertEqual(proteins, [
                ["sample", "lineage", "S", "D", "F", "M", "N"],
                ["protein", "test_odb12", "80.00", "10.00", "5.00", "5.00", "100"]
            ])

def write_fastqc_zip(outputDir, readName, flags, totalReads, lengthCounts, version="0.12.1"):
    '''
    Writes a minimal FastQC '_fastqc.zip' file in the FastQC 0.12.1 format, or the 0.13.0
    format (which no longer reports sequences flagged as poor quality).
    
    Parameters:
        flags -- a list of 'pass', 'warn' or 'fail' values, one per module in MODULES
        lengthCounts -- a list of (length, count) tuples where length can be a range like '35-39'
        version -- OPTIONAL; a string of '0.12.1' or '0.13.0'
    '''
    MODULES = ["Basic Statistics", "Per base sequence quality", "Per sequence quality scores",
               "Sequence Length Distribution", "Adapter Content"]
    lines = [f"##FastQC\t{version}"]
    for module, flag in zip(MODULES, flags):
        lines.append(f">>{module}\t{flag}")
        if module == "Basic Statistics":
            lines += ["#Measure\tValue", f"Filename\t{readName}.fq.gz", f"Total Sequences\t{totalReads}",
                      "Total Bases\t1 Mbp"]
            if version == "0.12.1":
                lines.append("Sequences flagged as poor quality\t3")
            lines.append("%GC\t38")
        elif module == "Sequence Length Distribution":
            lines += ["#Length\tCount"] + [ f"{length}\t{count}" for length, count in lengthCounts ]
        lines.append(">>END_MODULE")
    
    zipFile = os.path.join(outputDir, f"{readName}_fastqc.zip")
    with zipfile.ZipFile(zipFile, "w") as archive:
        archive.writestr(f"{readName}_fastqc/fastqc_data.txt", "\n".join(lines) + "\n")
    return zipFile

class TestFastQCTabulate(unittest.TestCase):
    def write_runs(self, tmpDir):
        '''
        Returns:
            runsDir -- a directory holding a paired-end 'pairedA' and single-end 'singleB' sample
        '''
        runsDir = os.path.join(tmpDir, "fastqc")
        pairedDir, singleDir = os.path.join(runsDir, "pairedA"), os.path.join(runsDir, "singleB")
        os.makedirs(pairedDir)
        os.makedirs(singleDir)
        write_fastqc_zip(pairedDir, "pairedA_1", ["pass", "pass", "warn", "pass", "fail"], 100, [("10-20", 50), ("100", 50)])
        write_fastqc_zip(pairedDir, "pairedA_2", ["pass", "warn", "warn", "pass", "fail"], 100, [("100", 100)])
        write_fastqc_zip(singleDir, "singleB", ["pass", "pass", "pass", "pass", "pass"], 10, [("50", 10)])
        return runsDir

    def test_parse_fastqc_zip(self):
        "Test that read counts, approximate base pairs and flags are parsed from a FastQC zip"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            zipFile = write_fastqc_zip(tmpDir, "reads_1", ["pass", "warn", "fail", "pass", "pass"], 100, [("10-20", 50), ("100", 50)])

            # Act
            stats = parse_fastqc_zip(zipFile)

            # Assert
            self.assertEqual(stats["total_reads"], 100)
            self.assertEqual(stats["low_quality_reads"], 3)
            self.assertEqual(stats["approx_bp"], 15 * 50 + 100 * 50, "Length ranges should use their midpoint")
            self.assertEqual(list(stats)[3:], ["basic_statistics", "per_base_sequence_quality",
                                               "per_sequence_quality_scores", "sequence_length_distribution",
                                               "adapter_content"])
            self.assertEqual(stats["per_base_sequence_quality"], "WARNING")
            self.assertEqual(stats["per_sequence_quality_scores"], "FAIL")

    def test_compare_reads(self):
        "Test that read pairs are summed and their flags compared"
        # Arrange
        forward = { "total_reads": 10, "approx_bp": 1000, "low_quality_reads": 1, "a": "PASS", "b": "WARNING" }
        reverse = { "total_reads": 10, "approx_bp": 900, "low_quality_reads": 0, "a": "PASS", "b": "FAIL" }

        # Act
        paired = compare_reads(forward, reverse, ["a", "b", "c"])
        single = compare_reads(forward, None, ["a", "b", "c"])

        # Assert
        self.assertEqual(paired, { "total_reads": 20, "approx_bp": 1900, "low_quality_reads": 1,
                                   "a": "SAME", "b": "DIFFERS", "c": "." })
        self.assertEqual(single, { "total_reads": 10, "approx_bp": 1000, "low_quality_reads": 1,
                                   "a": ".", "b": ".", "c": "." })

    def test_compare_reads_unreported(self):
        "Test that statistics FastQC 0.13 doesn't report (None) are '.' rather than summed"
        # Arrange
        forward = { "total_reads": 10, "approx_bp": 1000, "low_quality_reads": None, "a": "PASS" }
        reverse = { "total_reads": 10, "approx_bp": 900, "low_quality_reads": None, "a": "PASS" }

        # Act
        paired = compare_reads(forward, reverse, ["a"])
        single = compare_reads(forward, None, ["a"])

        # Assert
        self.assertEqual(paired, { "total_reads": 20, "approx_bp": 1900, "low_quality_reads": ".", "a": "SAME" })
        self.assertEqual(single, { "total_reads": 10, "approx_bp": 1000, "low_quality_reads": ".", "a": "." })

    @unittest.skipIf(not EXCEL_AVAILABLE, "xlsxwriter and openpyxl are needed to write and read workbooks")
    def test_tabulate_fastqc_0_13(self):
        "Test that FastQC 0.13 results, which lack sequences flagged as poor quality, are tabulated"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            pairedDir, singleDir = os.path.join(tmpDir, "pairedA"), os.path.join(tmpDir, "singleB")
            os.makedirs(pairedDir)
            os.makedirs(singleDir)
            write_fastqc_zip(pairedDir, "pairedA_1", ["pass"] * 5, 100, [("100", 100)], version="0.13.0")
            write_fastqc_zip(pairedDir, "pairedA_2", ["pass"] * 5, 100, [("100", 100)], version="0.13.0")
            write_fastqc_zip(singleDir, "singleB", ["pass"] * 5, 10, [("50", 10)], version="0.13.0")
            outputFile = os.path.join(tmpDir, "fastqc.xlsx")
            cmd = [ sys.executable, os.path.join(baseDir, "annotarium.py"), "qc", "fastqc", "tabulate",
                    "-i", pairedDir, singleDir, "-o", outputFile ]

            # Act
            returncode, stdout, stderr = run_subprocess(cmd)

            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            workbook = openpyxl.load_workbook(outputFile)
            combined = list(workbook["Combined Reads"].iter_rows(values_only=True))
            self.assertEqual(combined[1][:4], ("pairedA", 200, 20000, "."))
            self.assertEqual(combined[2][:4], ("singleB", 10, 500, "."))
            forward = list(workbook["Forward Reads"].iter_rows(values_only=True))
            self.assertEqual(forward[1][:4], ("pairedA", 100, 10000, "."))

    @unittest.skipIf(not EXCEL_AVAILABLE, "xlsxwriter and openpyxl are needed to write and read workbooks")
    def test_tabulate(self):
        "Test that a combined sheet with sample names comes first when any reads are paired"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            runsDir = self.write_runs(tmpDir)
            namesFile = os.path.join(tmpDir, "names.tsv")
            with open(namesFile, "w") as fileOut:
                fileOut.write("pairedA\tSample A\nsingleB\tSample B\n")
            outputFile = os.path.join(tmpDir, "fastqc") # .xlsx should be appended
            cmd = [ sys.executable, os.path.join(baseDir, "annotarium.py"), "qc", "fastqc", "tabulate",
                    "-i", runsDir, "-o", outputFile, "-s", namesFile ]

            # Act
            returncode, stdout, stderr = run_subprocess(cmd)

            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            workbook = openpyxl.load_workbook(outputFile + ".xlsx")
            self.assertEqual(workbook.sheetnames, ["Combined Reads", "Forward Reads", "Reverse Reads"])
            combined = list(workbook["Combined Reads"].iter_rows(values_only=True))
            self.assertEqual(combined[0][:5], ("read_prefix", "sample_name", "total_reads", "approx_bp", "low_quality_reads"))
            self.assertEqual(combined[1], ("pairedA", "Sample A", 200, 15750, 6, "SAME", "DIFFERS", "SAME", "SAME", "SAME"))
            self.assertEqual(combined[2], ("singleB", "Sample B", 10, 500, 3, ".", ".", ".", ".", "."))
            forward = list(workbook["Forward Reads"].iter_rows(values_only=True))
            self.assertNotIn("sample_name", forward[0], "Only the first sheet should have sample names")
            self.assertEqual(len(list(workbook["Reverse Reads"].iter_rows())), 2, "Only paired samples have reverse reads")

    @unittest.skipIf(not EXCEL_AVAILABLE, "xlsxwriter and openpyxl are needed to write and read workbooks")
    def test_tabulate_single_end(self):
        "Test that single-end reads only get a forward sheet, and unnamed read prefixes cause an error"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            runsDir = self.write_runs(tmpDir)
            namesFile = os.path.join(tmpDir, "names.tsv")
            with open(namesFile, "w") as fileOut:
                fileOut.write("singleB\tSample B\n")
            singleFile, pairedFile = os.path.join(tmpDir, "single.xlsx"), os.path.join(tmpDir, "paired.xlsx")
            baseCmd = [ sys.executable, os.path.join(baseDir, "annotarium.py"), "qc", "fastqc", "tabulate", "-s", namesFile ]

            # Act
            returncode, stdout, stderr = run_subprocess(baseCmd + [ "-i", os.path.join(runsDir, "singleB"), "-o", singleFile ])
            unnamedReturncode, unnamedStdout, unnamedStderr = run_subprocess(baseCmd + [ "-i", runsDir, "-o", pairedFile ])

            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            workbook = openpyxl.load_workbook(singleFile)
            self.assertEqual(workbook.sheetnames, ["Forward Reads"])
            self.assertEqual(list(workbook["Forward Reads"].iter_rows(values_only=True))[1][:3], ("singleB", "Sample B", 10))
            self.assertNotEqual(unnamedReturncode, 0, "Expected an error for a read prefix without a sample name")
            self.assertIn("pairedA", unnamedStderr)
            self.assertFalse(os.path.exists(pairedFile))

class TestFastaTabulate(unittest.TestCase):
    def write_stats(self, fastaContent, statsFile):
        '''
        Writes a FASTA file then runs 'fasta stats' on it to produce statsFile.
        '''
        fastaFile = statsFile + ".fasta"
        with open(fastaFile, "w") as fileOut:
            fileOut.write(fastaContent)
        cmd = [ sys.executable, os.path.join(baseDir, "annotarium.py"), "fasta", "stats", "-i", fastaFile, "-o", statsFile ]
        returncode, stdout, stderr = run_subprocess(cmd)
        self.assertEqual(returncode, 0, f"Expected no error from 'fasta stats' but got: {stderr}")

    def test_parse_fasta_stats_file(self):
        "Test that 'fasta stats' output is parsed with N base pairs tallied"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            statsFile = os.path.join(tmpDir, "genome.stats")
            self.write_stats(">a\nACGTNNNNac\n>b\nNNACG\n>c\nACGTACGTACGTAAAA\n", statsFile)

            # Act
            stats = parse_fasta_stats_file(statsFile)

            # Assert
            self.assertEqual(stats, { "total_bp": 31, "num_contigs": 3, "shortest_contig": 5, "longest_contig": 16,
                                      "n50_contig": 16, "median_contig": 10, "mean_contig": 10, "n_bp": 6 })

    def test_tabulate(self):
        "Test that statistics files are found recursively, limited by depth, and tabulated by sample"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            nestedDir = os.path.join(tmpDir, "runs", "nested")
            os.makedirs(nestedDir)
            self.write_stats(">x\nACGTNNNN\n", os.path.join(tmpDir, "runs", "b.stats"))
            self.write_stats(">y\nACGT\n>z\nNN\n", os.path.join(nestedDir, "a.stats"))
            oldFile = os.path.join(tmpDir, "old.stats")
            with open(oldFile, "w") as fileOut:
                fileOut.write("Genome size (bp): 63\nNumber of contigs: 1\n")
            outputFile = os.path.join(tmpDir, "table.tsv")
            cmd = [ sys.executable, os.path.join(baseDir, "annotarium.py"), "qc", "fasta", "tabulate",
                    "-i", os.path.join(tmpDir, "runs"), "-o", outputFile ]

            # Act
            returncode, stdout, stderr = run_subprocess(cmd)
            shallowFiles = find_files_with_suffix(os.path.join(tmpDir, "runs"), ".stats", maxDepth=1)

            # Assert
            self.assertEqual(returncode, 0, f"Expected no error but got: {stderr}")
            with open(outputFile, "r") as fileIn:
                table = [ line.rstrip("\n").split("\t") for line in fileIn ]
            self.assertEqual(table, [
                ["sample", "total_bp", "num_contigs", "shortest_contig", "longest_contig",
                 "n50_contig", "median_contig", "mean_contig", "n_bp"],
                ["a", "6", "2", "2", "4", "4", "3", "3", "2"],
                ["b", "8", "1", "8", "8", "8", "8", "8", "4"]
            ])
            self.assertEqual([ os.path.basename(f) for f in shallowFiles ], ["b.stats"])
            with self.assertRaises(ValueError):
                parse_fasta_stats_file(oldFile)

class TestFastaGrabber(unittest.TestCase):
    def write_fastas(self, tmpDir):
        fastaFiles = []
        for name, seq in [ ("a.fasta", "AAAA"), ("b.fasta", "CCCC") ]:
            fastaFile = os.path.join(tmpDir, name)
            with open(fastaFile, "w") as fileOut:
                fileOut.write(f">chr1 shared\n{seq}\n>only_{name[0]}\n{seq}\n")
            fastaFiles.append(fastaFile)
        return fastaFiles

    def run_grabber(self, fastaFiles, values, outputFile, extraArgs=[]):
        cmd = [ sys.executable, os.path.join(baseDir, "annotarium.py"), "fasta", "grabber",
                "-i", *fastaFiles, "-s", *values, "-o", outputFile ] + extraArgs
        return run_subprocess(cmd)

    def test_parse_file_contig_tsv(self):
        "Test that two column TSVs are parsed and that list files are not mistaken for them"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            tsvFile = os.path.join(tmpDir, "pairs.tsv")
            listFile = os.path.join(tmpDir, "ids.txt")
            with open(tsvFile, "w") as fileOut:
                fileOut.write("#fasta_file\tcontig_id\na.fasta\tchr1\nb.fasta\tonly_b\n")
            with open(listFile, "w") as fileOut:
                fileOut.write("chr1\nonly_b\n")

            # Act
            pairs = parse_file_contig_tsv(tsvFile)
            notPairs = parse_file_contig_tsv(listFile)

            # Assert
            self.assertEqual(pairs, [ ("a.fasta", "chr1"), ("b.fasta", "only_b") ])
            self.assertIsNone(notPairs)

    def test_multiple_files(self):
        "Test that IDs unique across files are grabbed, while shared IDs need a TSV naming their file"
        with tempfile.TemporaryDirectory() as tmpDir:
            # Arrange
            fastaFiles = self.write_fastas(tmpDir)
            tsvFile = os.path.join(tmpDir, "pairs.tsv")
            with open(tsvFile, "w") as fileOut:
                fileOut.write("b.fasta\tchr1\n")
            badTsvFile = os.path.join(tmpDir, "bad.tsv")
            with open(badTsvFile, "w") as fileOut:
                fileOut.write("c.fasta\tchr1\n")
            outputs = [ os.path.join(tmpDir, f"out{i}.fasta") for i in range(4) ]

            # Act
            uniqueResult = self.run_grabber(fastaFiles, ["only_a", tsvFile], outputs[0])
            ambiguousResult = self.run_grabber(fastaFiles, ["chr1"], outputs[1])
            badNameResult = self.run_grabber(fastaFiles, [badTsvFile], outputs[2])
            removeResult = self.run_grabber(fastaFiles, ["only_a"], outputs[3], ["-b", "remove"])

            # Assert
            self.assertEqual(uniqueResult[0], 0, f"Expected no error but got: {uniqueResult[2]}")
            self.assertEqual([ (record.id, str(record)) for record in parse_fasta(outputs[0]) ],
                             [ ("only_a", "AAAA"), ("chr1", "CCCC") ])
            self.assertNotEqual(ambiguousResult[0], 0, "Expected an error for an ID in multiple files")
            self.assertIn("more than one -i FASTA file", ambiguousResult[2])
            self.assertNotEqual(badNameResult[0], 0, "Expected an error for a TSV file name not given to -i")
            self.assertIn("do not match any -i file name", badNameResult[2])
            self.assertNotEqual(removeResult[0], 0, "Expected an error for duplicate output IDs")
            self.assertIn("would be written more than once", removeResult[2])

if __name__ == '__main__':
    unittest.main()
