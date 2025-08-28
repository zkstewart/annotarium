#! python3
# annotarium.py
# Front-end interface for GFF3-handling tools
# that have been developed over several years in the
# Genome_analysis_scripts and Various_scripts Z.K.S
# repository, but which need to be consolidated and
# made more user-friendly.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.validation import validate_f, validate_f_stats, \
    validate_g, validate_g_stats, \
    validate_g_to, validate_g_to_fasta
from modules.fasta import fasta_stats
from modules.gff3 import gff3_stats, gff3_to_fasta
from _version import __version__

def main():
    usage = """%(prog)s encapsulates a variety of annotation analysis tools for
    working with annotations in GFF3 format.
    """
    # Establish main parser
    p = argparse.ArgumentParser()
    
    # Set arguments shared by subparsers
    p.add_argument("-v", "--version",
                   action="version",
                   version="annotarium.py {version}".format(version=__version__))
    
    # Establish subparsers
    subParentParser = argparse.ArgumentParser(description=usage)
    subParentParser.add_argument("-v", "--version",
                                 action="version",
                                 version="annotarium.py {version}".format(version=__version__))
    
    subparsers = subParentParser.add_subparsers(dest="mode",
                                                required=True)
    
    # FASTA subparser
    fparser = subparsers.add_parser("fasta",
                                    parents=[p],
                                    add_help=False,
                                    help="FASTA file handling")
    fparser.set_defaults(func=fmain)
    
    subFASTAParsers = fparser.add_subparsers(dest="fastaMode",
                                             required=True)
    
    # FASTA > stats mode
    fstatsparser = subFASTAParsers.add_parser("stats",
                                             parents=[p],
                                             add_help=False,
                                             help="Report FASTA statistics")
    fstatsparser.add_argument("-i", dest="fastaFile",
                              required=True,
                              help="Location of FASTA file")
    fstatsparser.add_argument("-o", dest="outputFileName",
                              required=True,
                              help="Location to write statistics output")
    
    # GFF3 subparser
    gparser = subparsers.add_parser("gff3",
                                    parents=[p],
                                    add_help=False,
                                    help="GFF3 file handling")
    gparser.set_defaults(func=gmain)
    
    subGFF3Parsers = gparser.add_subparsers(dest="gff3Mode",
                                            required=True)
    
    # GFF3 > stats mode
    gstatsparser = subGFF3Parsers.add_parser("stats",
                                             parents=[p],
                                             add_help=False,
                                             help="Report GFF3 statistics")
    gstatsparser.add_argument("-i", dest="gff3File",
                              required=True,
                              help="Location of GFF3 file")
    gstatsparser.add_argument("-o", dest="outputFileName",
                              required=True,
                              help="Location to write statistics output")
    
    # GFF3 > to subparser
    gff3toparser = subGFF3Parsers.add_parser("to",
                                             parents=[p],
                                             add_help=False,
                                             help="GFF3 conversion")
    gff3toparser.set_defaults(func=fmain)
    
    subGFF3ToParsers = gff3toparser.add_subparsers(dest="gff3ToMode",
                                                   required=True)
    
    # GFF3 > to > fasta mode
    gtofparser = subGFF3ToParsers.add_parser("fasta",
                                             parents=[p],
                                             add_help=False,
                                             help="GFF3 to FASTA conversion")
    gtofparser.add_argument("-i", dest="gff3File",
                            required=True,
                            help="Location of GFF3 file")
    gtofparser.add_argument("-f", dest="fastaFile",
                            required=True,
                            help="Location of FASTA (e.g., genome) file")
    gtofparser.add_argument("-o", dest="outputFilePrefix",
                            required=True,
                            help="Prefix to output files (suffixes will be set by output type)")
    gtofparser.add_argument("--features", dest="features",
                            required=False,
                            nargs="+",
                            help="""Specify one or more top level feature types to output sequences for;
                            child features with exons will be found, and if multiple exist only the
                            longest will be output. Hence, specifying 'gene' will output representative
                            mRNAs, whereas specifying 'mRNA' will output each isoform. Default ==
                            'mRNA'""",
                            default=["mRNA"])
    gtofparser.add_argument("--types", dest="types",
                            required=False,
                            nargs="+",
                            choices=["exon", "cds", "protein", "all"],
                            help="""Specify one or more types of sequence to output; each type will
                            be output into separate files, and only features which have all
                            available types will be emitted e.g., exon-only features will be
                            skipped if any other types are requested. Default == 'all' which 
                            is a shortcut to set the three other types.""",
                            default=["all"])
    gtofparser.add_argument("--translation", dest="translationTable",
                            required=False,
                            type=int,
                            help="""Specify the NCBI numeric genetic code to utilise for CDS
                            translation (if relevant); this should be an integer from 1 to 31
                            (default == 1 i.e., Standard Code)""",
                            default=1)
    
    args = subParentParser.parse_args()
    
    # Split into mode-specific functions
    if args.mode == "fasta":
        print("## annotarium.py - FASTA handling ##")
        validate_f(args)
        fmain(args)
    elif args.mode == "gff3":
        print("## annotarium.py - GFF3 handling ##")
        validate_g(args)
        gmain(args)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def fmain(args):
    # Split into sub-mode-specific functions
    if args.fastaMode == "stats":
        print("## FASTA statistics ##")
        validate_f_stats(args)
        fasta_stats(args)
    
    print("FASTA handling complete!")

def gmain(args):
    # Split into sub-mode-specific functions
    if args.gff3Mode == "stats":
        print("## GFF3 statistics ##")
        validate_g_stats(args)
        gff3_stats(args)
    elif args.gff3Mode == "to":
        validate_g_to(args)
        if args.gff3ToMode == "fasta":
            print("## GFF3 to FASTA conversion ##")
            validate_g_to_fasta(args)
            gff3_to_fasta(args)
    
    print("GFF3 handling complete!")

if __name__ == "__main__":
    main()
