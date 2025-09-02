import os, re

def parse_regions(regions):
    '''
    Returns:
        parsedRegions -- a list of dictionaries with structure like:
                         [
                             {
                                 "contig": contigID, # string
                                 "start": start, # int
                                 "end": end # int or None
                             }, { ... }, ...
                         ]
    '''
    # Parse regions
    parsedRegions = []
    regionsRegex = re.compile(r"^([^:]+):(\d+)-(\d+)$")
    for region in regions:
        reMatch = regionsRegex.match(region)
        
        # Handle chr:start-end format
        if reMatch != None:
            contigID, start, end = reMatch.groups()
            start = int(start)
            end = int(end)
            
            # Detect and fix reverse orientation
            if start > end:
                start, end = end, start
            
            # Store region
            parsedRegions.append({"contig": contigID, "start": start, "end": end})
        
        # Handle chr format
        else:
            parsedRegions.append({"contig": region, "start": 1, "end": None}) # None is interpreted as no end
    
    # Handle empty regions
    if parsedRegions == []:
        parsedRegions = None # None is interpreted as no selection"
    
    return parsedRegions

def validate_f(args):
    '''
    Validation for arguments common to all "fasta" mode commands.
    '''
    # Validate FASTA file
    args.fastaFile = os.path.abspath(args.fastaFile)
    if not os.path.isfile(args.fastaFile):
        raise FileNotFoundError(f"FASTA file (-i {args.fastaFile}) does not exist!")

def validate_f_stats(args):
    '''
    Validation for arguments used in "fasta stats" mode.
    '''
    # Validate output file name
    if args.outputFileName != None:
        args.outputFileName = os.path.abspath(args.outputFileName)
        if os.path.exists(args.outputFileName):
            raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

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

def validate_g_merge(args):
    '''
    Validation for arguments used in "gff3 merge" mode.
    '''
    # Validate second GFF3 file
    args.gff3File2 = os.path.abspath(args.gff3File2)
    if not os.path.isfile(args.gff3File2):
        raise FileNotFoundError(f"Second GFF3 file (-2 {args.gff3File2}) does not exist!")
    
    # Validate numeric arguments
    if not 0 <= args.isoformPercent <= 1:
        raise ValueError("--isoformPercent must be in the range of 0.0 to 1.0 inclusive")
    if not 0 <= args.duplicatePercent <= 1:
        raise ValueError("--duplicatePercent must be in the range of 0.0 to 1.0 inclusive")
    
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")
    
    # Validate optional details file
    if args.outputDetailsName != None:
        if os.path.exists(args.outputDetailsName):
            raise FileExistsError(f"Details output file (--outputDetails {args.outputDetailsName}) already exists!")

def validate_g_filter(args):
    '''
    Validation for arguments used in "gff3 filter" mode.
    '''
    # Validate argument logic
    if args.regions == [] and args.listFile == None and args.values == []:
        raise ValueError("'gff3 filter' needs at least one of --regions or --list or --values to be set to do anything")
    
    # Parse regions
    args.regions = parse_regions(args.regions)
    
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")
    
    # Validate optional list file
    if args.listFile != None:
        if not os.path.isfile(args.listFile):
            raise FileExistsError(f"Values list file (--list {args.listFile}) is not a file")

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

def validate_g_to_tsv(args):
    '''
    Validation for arguments used in "gff3 to tsv" mode.
    '''
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_g_to_gff3(args):
    '''
    Validation for arguments used in "gff3 to gff3" mode.
    '''
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")
