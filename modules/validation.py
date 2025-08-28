import os

def validate_f(args):
    '''
    Validation for arguments common to all "fasta" mode commands.
    '''
    # Validate FASTA file
    args.fastaFile = os.path.abspath(args.fastaFile)
    if not os.path.isfile(args.fastaFile):
        raise FileNotFoundError(f"FASTA file (-i {args.fastaFile}) does not exist!")
    
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_f_stats(args):
    '''
    Validation for arguments used in "fasta stats" mode.
    '''
    pass # no specific validation needed for this mode

def validate_g(args):
    '''
    Validation for arguments common to all "gff3" mode commands.
    '''
    # Validate GFF3 file
    args.gff3File = os.path.abspath(args.gff3File)
    if not os.path.isfile(args.gff3File):
        raise FileNotFoundError(f"GFF3 file (-i {args.gff3File}) does not exist!")

def validate_g_stats(args):
    '''
    Validation for arguments used in "gff3 stats" mode.
    '''
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_g_to(args):
    '''
    Validation for arguments common to all "gff3 to" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_g_to_fasta(args):
    '''
    Validation for arguments used in "gff3 to fasta" mode.
    '''
    # Validate FASTA file
    args.fastaFile = os.path.abspath(args.fastaFile)
    if not os.path.isfile(args.fastaFile):
        raise FileNotFoundError(f"FASTA file (-f {args.fastaFile}) does not exist!")

    # Expand args.types shortcut
    if "all" in args.types:
        args.types = ["exon", "cds", "protein"]
    args.types = list(set(args.types)) # remove any duplicate arguments; order doesn't matter
    
    # Format output file names based on args.types
    args.outputFilePrefix = os.path.abspath(args.outputFilePrefix)
    args.outputFileNames = { "exon": None, "cds": None, "protein": None }
    if "exon" in args.types:
        args.outputFileNames["exon"] = args.outputFilePrefix + ".exon.fasta"
    if "cds" in args.types:
        args.outputFileNames["cds"] = args.outputFilePrefix + ".cds.fasta"
    if "protein" in args.types:
        args.outputFileNames["protein"] = args.outputFilePrefix + ".protein.fasta"
    
    # Validate that output files do not already exist
    for seqType, outputFileName in args.outputFileNames.items():
        if os.path.exists(outputFileName):
            raise FileExistsError(f"--types {seqType} output '{outputFileName}' already exists!")
    
    # Ensure that translationTable value is sensible
    if args.translationTable < 0:
        raise ValueError("--translation should be given a value >= 0")
    if args.translationTable in [7, 8, 17, 18, 19, 20]:
        raise ValueError(f"--translation {args.translationTable} is not a valid NCBI table")
    if args.translationTable > 33:
        raise ValueError("--translation was given a value > 33; this is not a valid NCBI table")
    