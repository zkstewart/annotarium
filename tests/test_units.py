#!/usr/bin/env python3

import os, re, sys, unittest, shutil, subprocess, tempfile

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.fasta import FASTATarium, Sequence
from modules.msa import MAFFT
from modules.annotable import Annotable

# Specify data locations
baseDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
dataDir = os.path.join(baseDir, "tests", "data")
codonIndelsFile = os.path.join(dataDir, "codon_indels.fasta")
annotableFile = os.path.join(dataDir, "annotable.tsv")

MAFFT_EXE = shutil.which("mafft")

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

if __name__ == '__main__':
    unittest.main()
