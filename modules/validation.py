import os, re

class DirectoryNotFoundError(Exception):
    pass

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

def parse_annotate(columnAttributeDelimiter):
    '''
    Parameters:
        columnAttributeDelimiter -- a list of strings in format 'column:attribute:delimiter'
    Returns:
        annotateValues -- a list of lists like [column, attribute, delimiter] where delimiter
                          can be None
    '''
    annotateValues = []
    for value in columnAttributeDelimiter:
        column, attributeDelimiter = value.split(":", maxsplit=1)
        attribute, delimiter = attributeDelimiter.split(":", maxsplit=1)
        annotateValues.append([column, attribute, delimiter if delimiter != "" else None])
    return annotateValues

def validate_b(args):
    '''
    Validation for arguments common to all "blast" mode commands.
    '''
    # Validate output file name
    if args.outputFileName != None:
        args.outputFileName = os.path.abspath(args.outputFileName)
        if os.path.exists(args.outputFileName):
            raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_d(args):
    '''
    Validation for arguments common to all "domains" mode commands.
    '''
    # Validate domains file
    args.domainsFile = os.path.abspath(args.domainsFile)
    if not os.path.isfile(args.domainsFile):
        raise FileNotFoundError(f"Domains file (-i {args.domainsFile}) does not exist!")
    
    # Validate output file name
    if args.outputFileName != None:
        args.outputFileName = os.path.abspath(args.outputFileName)
        if os.path.exists(args.outputFileName):
            raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_d_resolve(args):
    '''
    Validation for arguments used in "domains resolve" mode.
    '''
    # Validate numeric arguments
    if args.evalue != None:
        if args.evalue < 0:
            raise ValueError("--evalue must >= 0")
    
    if args.overlapCutoff < 0:
        raise ValueError("--overlapPercent must >= 0")
    if args.overlapCutoff > 1:
        raise ValueError("--overlapPercent must <= 1")

def validate_b_to(args):
    '''
    Validation for arguments common to all "blast to" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_b_to_homologs(args):
    '''
    Validation for arguments used in "blast to homologs" mode.
    '''
    # Validate outfmt6 files
    args.inputFile1 = os.path.abspath(args.inputFile1)
    if not os.path.isfile(args.inputFile1):
        raise FileNotFoundError(f"outfmt6 file (-i1 {args.inputFile1}) does not exist!")
    
    args.inputFile2 = os.path.abspath(args.inputFile2)
    if not os.path.isfile(args.inputFile2):
        raise FileNotFoundError(f"outfmt6 file (-i2 {args.inputFile2}) does not exist!")
    
    # Validate numeric arguments
    if args.evalue != None:
        if args.evalue < 0:
            raise ValueError("--evalue must >= 0")

def validate_f(args):
    '''
    Validation for arguments common to all "fasta" mode commands.
    '''
    # Validate FASTA file
    args.fastaFile = os.path.abspath(args.fastaFile)
    if not os.path.isfile(args.fastaFile):
        raise FileNotFoundError(f"FASTA file (-i {args.fastaFile}) does not exist!")

def validate_f_softmask(args):
    '''
    Validation for arguments used in "fasta softmask" mode.
    '''
    # Validate output file name
    if args.outputFileName != None:
        args.outputFileName = os.path.abspath(args.outputFileName)
        if os.path.exists(args.outputFileName):
            raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_f_stats(args):
    '''
    Validation for arguments used in "fasta stats" mode.
    '''
    # Validate output file name
    if args.outputFileName != None:
        args.outputFileName = os.path.abspath(args.outputFileName)
        if os.path.exists(args.outputFileName):
            raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_f_explode(args):
    '''
    Validation for arguments used in "fasta explode" mode.
    '''
    # Validate output file name
    args.outputDirectory = os.path.abspath(args.outputDirectory)
    if not os.path.exists(args.outputDirectory):
        parentDir = os.path.dirname(args.outputDirectory)
        if not os.path.isdir(parentDir):
            raise DirectoryNotFoundError((f"Output directory (-o {args.outputDirectory}) cannot be created since " + 
                                          f"its parent directory ({parentDir}) does not exist. Create this location " +
                                          "first or specify a different output location"))
        os.mkdir(args.outputDirectory)
        print(f"# Created output directory (-o {args.outputDirectory}) as part of argument validation")

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

def validate_g_pcr(args):
    '''
    Validation for arguments used in "gff3 pcr" mode.
    '''
    # Validate FASTA file
    args.fastaFile = os.path.abspath(args.fastaFile)
    if not os.path.isfile(args.fastaFile):
        raise FileNotFoundError(f"FASTA file (-f {args.fastaFile}) does not exist!")
    
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_g_relabel(args):
    '''
    Validation for arguments used in "gff3 relabel" mode.
    '''
    # Validate that an argument was set which gives us something to do
    if args.listFile == None and args.contigSuffix == None and args.geneSuffix == None:
        raise ValueError("At least one of --list and/or --csuffix and/or --gsuffix must be set")
    
    # Validate that suffix doesn't break GFF3 formatting
    if args.contigSuffix != None:
        if "\t" in args.contigSuffix:
            raise ValueError("Tab characters cannot be used with --csuffix")
    if args.geneSuffix != None:
        if "\t" in args.geneSuffix:
            raise ValueError("Tab characters cannot be used with --gsuffix")
        if ";" in args.geneSuffix:
            raise ValueError("Semicolon characters cannot be used with --gsuffix")
    
    # Validate list file
    if args.listFile != None:
        args.listFile = os.path.abspath(args.listFile)
        if not os.path.isfile(args.listFile):
            raise FileNotFoundError(f"List file (--list {args.listFile}) does not exist!")
    
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

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

def validate_g_annotate(args):
    '''
    Validation for arguments used in "gff3 annotate" mode.
    '''
    # Validate table file name
    if not os.path.isfile(args.tableFile):
        raise FileExistsError(f"Table file (-t {args.tableFile}) is not a file")
    
    # Parse columnAttributeDelimiter
    args.columnAttributeDelimiter = parse_annotate(args.columnAttributeDelimiter)
    
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
        if outputFileName != None and os.path.exists(outputFileName):
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

def validate_p(args):
    '''
    Validation for arguments common to all "homologs" mode commands.
    '''
    # Validate homologs file
    args.homologsFile = os.path.abspath(args.homologsFile)
    if not os.path.isfile(args.homologsFile):
        raise FileNotFoundError(f"Homologs file (-i {args.homologsFile}) does not exist!")
    
    # Validate output file name
    if args.outputFileName != None:
        args.outputFileName = os.path.abspath(args.outputFileName)
        if os.path.exists(args.outputFileName):
            raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_p_annotate(args):
    '''
    Validation for arguments for "homologs annotate" mode commands.
    '''
    # Validate GFF3 files
    args.gff3File1 = os.path.abspath(args.gff3File1)
    if not os.path.isfile(args.gff3File1):
        raise FileNotFoundError(f"GFF3 file (-g1 {args.gff3File1}) does not exist!")
    
    args.gff3File2 = os.path.abspath(args.gff3File2)
    if not os.path.isfile(args.gff3File2):
        raise FileNotFoundError(f"GFF3 file (-g2 {args.gff3File2}) does not exist!")

def validate_p_to(args):
    '''
    Validation for arguments common to all "homologs to" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_p_to_bedpe(args):
    '''
    Validation for arguments for "homologs to bedpe" mode commands.
    '''
    # Validate GFF3 files
    args.gff3File1 = os.path.abspath(args.gff3File1)
    if not os.path.isfile(args.gff3File1):
        raise FileNotFoundError(f"GFF3 file (-g1 {args.gff3File1}) does not exist!")
    
    args.gff3File2 = os.path.abspath(args.gff3File2)
    if not os.path.isfile(args.gff3File2):
        raise FileNotFoundError(f"GFF3 file (-g2 {args.gff3File2}) does not exist!")

def validate_rnammer(args):
    '''
    Validation for arguments used in "rnammer" mode.
    '''
    # Validate RNAmmer GFF2 file
    args.rnammerGff2 = os.path.abspath(args.rnammerGff2)
    if not os.path.isfile(args.rnammerGff2):
        raise FileNotFoundError(f"RNAmmer GFF2 file (-i {args.rnammerGff2}) does not exist!")
    
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")

def validate_irf(args):
    '''
    Validation for arguments used in "irf" mode.
    '''
    # Validate IRF .dat file
    args.irfDatFile = os.path.abspath(args.irfDatFile)
    if not os.path.isfile(args.irfDatFile):
        raise FileNotFoundError(f"IRF .dat file (-i {args.irfDatFile}) does not exist!")
    
    # Validate numeric arguments
    if args.minimumLength < -1:
        raise ValueError("--minLen must be a positive integer, or 0 or -1")
    if args.maximumLength < -1:
        raise ValueError("--maxLen must be a positive integer, or 0 or -1")
    if args.minimumGap < -1:
        raise ValueError("--minGap must be a positive integer, or 0 or -1")
    if args.maximumGap < -1:
        raise ValueError("--maxGap must be a positive integer, or 0 or -1")
    if args.identityCutoff < 0:
        raise ValueError("--identity must be no less than 0")
    if args.identityCutoff > 100:
        raise ValueError("--identity must be no more than 100")
    
    # Validate output file name
    args.outputFileName = os.path.abspath(args.outputFileName)
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"Output file (-o {args.outputFileName}) already exists!")
