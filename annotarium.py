#! python3
# annotarium.py
# Front-end interface for GFF3-handling tools
# that have been developed over several years in the
# Genome_analysis_scripts and Various_scripts Z.K.S
# repository, but which need to be consolidated and
# made more user-friendly.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.validation import validate_f, validate_f_stats, validate_f_explode, \
    validate_g, validate_g_stats, validate_g_merge, validate_g_filter, validate_g_annotate, \
    validate_g_to, validate_g_to_tsv, validate_g_to_fasta, validate_g_to_gff3, \
    validate_rnammer
from modules.fasta import fasta_stats, fasta_explode
from modules.gff3 import gff3_stats, gff3_merge, gff3_filter, gff3_annotate, \
    gff3_to_fasta, gff3_to_tsv, gff3_to_gff3
from modules.rnammer import rnammer_reformat
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
    fstatsparser.add_argument("--out", "-o", dest="outputFileName",
                              required=False,
                              help="Optionally, write statistics output to file")
    
    # FASTA > explode mode
    fexplodeparser = subFASTAParsers.add_parser("explode",
                                                parents=[p],
                                                add_help=False,
                                                help="Explode FASTA into contigs")
    fexplodeparser.add_argument("-i", dest="fastaFile",
                                 required=True,
                                 help="Location of FASTA file")
    fexplodeparser.add_argument("-o", dest="outputDirectory",
                                 required=False,
                                 help="Directory to write contig files to")
    
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
    
    # GFF3 > merge mode
    gmergeparser = subGFF3Parsers.add_parser("merge",
                                             parents=[p],
                                             add_help=False,
                                             help="Merge GFF3 files")
    gmergeparser.add_argument("-i", dest="gff3File",
                              required=True,
                              help="Location of first GFF3 file")
    gmergeparser.add_argument("-2", dest="gff3File2",
                              required=True,
                              help="Location of second GFF3 file to merge into the first file")
    gmergeparser.add_argument("-o", dest="outputFileName",
                              required=True,
                              help="Location to write merged GFF3 output")
    gmergeparser.add_argument("--outputDetails", dest="outputDetailsName",
                              required=False,
                              help="Optional location to write merging details",
                              default=None)
    gmergeparser.add_argument("--isoformPercent", dest="isoformPercent",
                              required=False,
                              type=float,
                              help="""Specify the percentage overlap of two models before they are clustered
                              as isoforms; default == 0.3, equivalent to 30 percent.""",
                              default=0.3)
    gmergeparser.add_argument("--duplicatePercent", dest="duplicatePercent",
                              required=False,
                              type=float,
                              help="""Specify the percentage overlap of two models before they are considered
                              as duplicates and hence rejected or replaced; default == 0.6; equivalent to 60 percent.""",
                              default=0.6)
    
    # GFF3 > filter mode
    gfilterparser = subGFF3Parsers.add_parser("filter",
                                              parents=[p],
                                              add_help=False,
                                              help="Filter GFF3 file")
    gfilterparser.add_argument("-i", dest="gff3File",
                               required=True,
                               help="Location of GFF3 file")
    gfilterparser.add_argument("-r", dest="retrieveOrRemove",
                               required=True,
                               choices=["retrieve", "remove"],
                               help="Specify whether selected regions are retrieved or removed")
    gfilterparser.add_argument("-o", dest="outputFileName",
                               required=True,
                               help="Location to write filtered GFF3 output")
    gfilterparser.add_argument("--regions", dest="regions",
                               required=False,
                               nargs="+",
                               help="""Optionally, specify one or more regions to select with
                               contig:start-end format (e.g., 'chr1:1000000-2000000') or just
                               with the contig alone
                               """,
                               default=[])
    gfilterparser.add_argument("--list", dest="listFile",
                               required=False,
                               help="Optional location of a text file to parse for listed values",
                               default=None)
    gfilterparser.add_argument("--values", dest="values",
                               required=False,
                               nargs="+",
                               help="Optionally, specify one or more values to select",
                               default=[])
    
    # GFF3 > annotate mode
    gannotateparser = subGFF3Parsers.add_parser("annotate",
                                                parents=[p],
                                                add_help=False,
                                                help="Annotate attributes into a GFF3 file")
    gannotateparser.add_argument("-i", dest="gff3File",
                                 required=True,
                                 help="Location of GFF3 file")
    gannotateparser.add_argument("-t", dest="tableFile",
                                 required=True,
                                 help="Location of table file")
    gannotateparser.add_argument("-c", dest="columnAttributeDelimiter",
                                 required=True,
                                 nargs="+",
                                 help="""Specify one or more 'tableColumn:attributeKey:delimiter' trios, where
                                 the -t column will be mapped as a GFF3 attribute with 'attributeKey=value'
                                 format. The delimiter will determine whether multiple values exist in the table
                                 column and to only return the first, or if you leave it blank like
                                 'gene_name:Name:' then no splitting will occur and the entire tableColumn
                                 value will be used as-is.""")
    gannotateparser.add_argument("-o", dest="outputFileName",
                                 required=True,
                                 help="Location to write modified GFF3 output")
    
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
    
    # GFF3 > to > TSV mode
    gtotsvparser = subGFF3ToParsers.add_parser("tsv",
                                               parents=[p],
                                               add_help=False,
                                               help="GFF3 to TSV conversion")
    gtotsvparser.add_argument("-i", dest="gff3File",
                              required=True,
                              help="Location of GFF3 file")
    gtotsvparser.add_argument("-o", dest="outputFileName",
                              required=True,
                              help="Location to write TSV output")
    gtotsvparser.add_argument("-forEach", dest="forEach",
                              required=True,
                              help="""The value here should relate to any feature type in column 3 of
                              the GFF3; examples include 'gene' or 'mRNA'.""")
    gtotsvparser.add_argument("-map", dest="map",
                              required=True,
                              help="""A valid key that will be the left-most column and serve as a unique
                              value to anchor other details""")
    gtotsvparser.add_argument("-to", dest="to",
                              required=True,
                              nargs="+",
                              help="""One or more valid keys separated with a space; each key should correspond
                              to a GFF3 column or attribute value, and output columns will be based on the order
                              provided herein""")
    gtotsvparser.add_argument("--noHeader", dest="noHeader",
                              required=False,
                              action="store_true",
                              help="Set this flag to omit the header row",
                              default=False)
    gtotsvparser.add_argument("--null", dest="nullChar",
                              required=False,
                              help="""Optionally, specify the character(s) used to denote a lack of 'map' value""",
                              default="_")
    gtotsvparser.add_argument("--sep", dest="sepChar",
                              required=False,
                              help="""Optionally, specify the character(s) used to separate multiple values of
                              the same 'map' key""",
                              default=";")
    
    # GFF3 > to > GFF3 mode
    gtogff3parser = subGFF3ToParsers.add_parser("gff3",
                                                parents=[p],
                                                add_help=False,
                                                help="GFF3 to GFF3 reformatting")
    gtogff3parser.add_argument("-i", dest="gff3File",
                               required=True,
                               help="Location of GFF3 file")
    gtogff3parser.add_argument("-o", dest="outputFileName",
                               required=True,
                               help="Location to write GFF3 output")
    
    # RNAmmer mode
    rnammerparser = subparsers.add_parser("rnammer",
                                          parents=[p],
                                          add_help=False,
                                          help="Reformat RNAmmer to GFF3")
    rnammerparser.set_defaults(func=rnammermain)
    rnammerparser.add_argument("-i", dest="rnammerGff2",
                               required=True,
                               help="Location of RNAmmer GFF2 file")
    rnammerparser.add_argument("-o", dest="outputFileName",
                               required=True,
                               help="Location to write GFF3 output")
    
    args = subParentParser.parse_args()
    
    # Split into mode-specific functions
    if args.mode == "fasta":
        print("## annotarium.py - FASTA handling ##")
        validate_f(args)
        fmain(args)
    if args.mode == "gff3":
        print("## annotarium.py - GFF3 handling ##")
        validate_g(args)
        gmain(args)
    if args.mode == "rnammer":
        print("## annotarium.py - RNAmmer handling ##")
        validate_rnammer(args)
        rnammermain(args)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def fmain(args):
    # Split into sub-mode-specific functions
    if args.fastaMode == "stats":
        print("## FASTA statistics ##")
        validate_f_stats(args)
        fasta_stats(args)
    if args.fastaMode == "explode":
        print("## FASTA explosion ##")
        validate_f_explode(args)
        fasta_explode(args)
    
    print("FASTA handling complete!")

def gmain(args):
    # Split into sub-mode-specific functions
    if args.gff3Mode == "stats":
        print("## GFF3 statistics ##")
        validate_g_stats(args)
        gff3_stats(args)
    if args.gff3Mode == "merge":
        print("## GFF3 merge ##")
        validate_g_merge(args)
        gff3_merge(args)
    if args.gff3Mode == "filter":
        print("## GFF3 filtering ##")
        validate_g_filter(args) # sets args.regions
        gff3_filter(args)
    if args.gff3Mode == "annotate":
        print("## GFF3 attribute annotation ##")
        validate_g_annotate(args) # sets args.columnAttributeDelimiter
        gff3_annotate(args)
    if args.gff3Mode == "to":
        validate_g_to(args)
        if args.gff3ToMode == "fasta":
            print("## GFF3 to FASTA conversion ##")
            validate_g_to_fasta(args)
            gff3_to_fasta(args)
        if args.gff3ToMode == "tsv":
            print("## GFF3 to TSV conversion ##")
            validate_g_to_tsv(args)
            gff3_to_tsv(args)
        if args.gff3ToMode == "gff3":
            print("## GFF3 to GFF3 re-formatter ##")
            validate_g_to_gff3(args)
            gff3_to_gff3(args)
    
    print("GFF3 handling complete!")

def rnammermain(args):
    rnammer_reformat(args)
    
    print("RNAmmer handling complete!")

if __name__ == "__main__":
    main()
