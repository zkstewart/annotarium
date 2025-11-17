#! python3
# annotarium.py
# Front-end interface for GFF3-handling tools
# that have been developed over several years in the
# Genome_analysis_scripts and Various_scripts Z.K.S
# repository, but which need to be consolidated and
# made more user-friendly.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.validation import validate_b, validate_b_to, validate_b_to_homologs, \
    validate_d, validate_d_resolve, \
    validate_f, validate_f_softmask, validate_f_stats, validate_f_explode, \
    validate_g, validate_g_stats, validate_g_merge, validate_g_filter, validate_g_annotate, validate_g_pcr, validate_g_relabel, \
    validate_g_to, validate_g_to_tsv, validate_g_to_fasta, validate_g_to_gff3, \
    validate_g_mp, validate_g_mp_reformat, validate_g_mp_resolve, \
    validate_p, validate_p_annotate, validate_p_to, validate_p_to_bedpe, \
    validate_rnammer, \
    validate_irf
from modules.blast import blast_to_homologs
from modules.domains import domains_resolve
from modules.fasta import fasta_softmask_to_bed, fasta_stats, fasta_explode
from modules.gff3 import gff3_stats, gff3_merge, gff3_filter, gff3_annotate, gff3_pcr, gff3_relabel, \
    gff3_to_fasta, gff3_to_tsv, gff3_to_gff3, \
    gff3_mp_reformat, gff3_mp_resolve
from modules.homologs import homologs_annotate, homologs_to_bedpe
from modules.irf import irf_to_gff3
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
    
    # Apollo subparser
    aparser = subparsers.add_parser("apollo",
                                    parents=[p],
                                    add_help=False,
                                    help="Apollo file handling")
    aparser.set_defaults(func=fmain)
    
    subApolloParsers = aparser.add_subparsers(dest="apolloMode",
                                              required=True)
    
    # Apollo > rename mode
    arenameparser = subApolloParsers.add_parser("rename",
                                                parents=[p],
                                                add_help=False,
                                                help="Rename Apollo annotations")
    
    # Apollo > update mode
    aupdateparser = subApolloParsers.add_parser("update",
                                                parents=[p],
                                                add_help=False,
                                                help="Update genome based on Apollo artifacts")
    
    # Apollo > pipeline mode
    apipelineparser = subApolloParsers.add_parser("pipeline",
                                                  parents=[p],
                                                  add_help=False,
                                                  help="Automatically integrate Apollo GFF3 into existing GFF3")
    
    # Blast subparser
    bparser = subparsers.add_parser("blast",
                                    parents=[p],
                                    add_help=False,
                                    help="BLAST (/MMseqs2) handling")
    bparser.set_defaults(func=bmain)
    
    subBlastParsers = bparser.add_subparsers(dest="blastMode",
                                             required=True)
    
    # Blast > to subparser
    btoparser = subBlastParsers.add_parser("to",
                                           parents=[p],
                                           add_help=False,
                                           help="Blast results file conversion")
    btoparser.set_defaults(func=bmain)
    
    subBlastToParsers = btoparser.add_subparsers(dest="blastToMode",
                                                 required=True)
    
    # Blast > to > homologs mode
    btohparser = subBlastToParsers.add_parser("homologs",
                                              parents=[p],
                                              add_help=False,
                                              help="Blast to homologs TSV (reciprocal best) conversion")
    btohparser.add_argument("-i1", dest="inputFile1",
                            required=True,
                            help="Location of first outfmt6 file")
    btohparser.add_argument("-i2", dest="inputFile2",
                            required=True,
                            help="Location of second outfmt6 file")
    btohparser.add_argument("-o", dest="outputFileName",
                            required=True,
                            help="Output TSV (2 columns; queryid targetid) file name")
    btohparser.add_argument("--evalue", dest="evalue",
                            required=False,
                            type=float,
                            help="Optionally ignore hits with worse than this evalue",
                            default=None)
    
    # Domains subparser
    dparser = subparsers.add_parser("domains",
                                    parents=[p],
                                    add_help=False,
                                    help="Protein domain prediction handling")
    dparser.set_defaults(func=dmain)
    
    subDomainsParsers = dparser.add_subparsers(dest="domainsMode",
                                               required=True)
    
    # Domains > resolve mode
    dresolveparser = subDomainsParsers.add_parser("resolve",
                                                  parents=[p],
                                                  add_help=False,
                                                  help="Resolve overlapping domain predictions")
    dresolveparser.add_argument("-i", dest="domainsFile",
                                required=True,
                                help="Location of .domtblout or .domains.tsv file")
    dresolveparser.add_argument("-o", dest="outputFileName",
                                required=False,
                                help="Write softmasked BED region output to file")
    dresolveparser.add_argument("--evalue", dest="evalue",
                                required=False,
                                type=float,
                                help="Optionally ignore hits with worse than this evalue",
                                default=None)
    dresolveparser.add_argument("--overlapPercent", dest="overlapCutoff",
                                required=False,
                                type=float,
                                help="""Specify the percentage value (as a ratio from 0 to 1) under which
                                overlaps are handled by trimming, and over which overlaps are handled
                                by culling the domain with worse E-value; default == 0.25, equivalent
                                to 25 percent""",
                                default=0.25)
    dresolveparser.add_argument("--tsv", dest="isDomainsTsvFormat",
                                required=False,
                                action="store_true",
                                help="""Optionally, provide this flag if the input domains file
                                has already been parsed into a .domains.tsv format""",
                                default=False)
    
    # FASTA subparser
    fparser = subparsers.add_parser("fasta",
                                    parents=[p],
                                    add_help=False,
                                    help="FASTA file handling")
    fparser.set_defaults(func=fmain)
    
    subFASTAParsers = fparser.add_subparsers(dest="fastaMode",
                                             required=True)
    
    # FASTA > softmask mode
    fsoftmaskparser = subFASTAParsers.add_parser("softmask",
                                                 parents=[p],
                                                 add_help=False,
                                                 help="Encode softmasked regions as BED 3-column format")
    fsoftmaskparser.add_argument("-i", dest="fastaFile",
                                 required=True,
                                 help="Location of FASTA file")
    fsoftmaskparser.add_argument("-o", dest="outputFileName",
                                 required=False,
                                 help="Write softmasked BED region output to file")
    
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
                              as duplicates; duplicates will not be merged from file 2 into file 1;
                              default == 0.6, equivalent to 60 percent.""",
                              default=0.6)
    
    # GFF3 > pcr mode
    gpcrparser = subGFF3Parsers.add_parser("pcr",
                                           parents=[p],
                                           add_help=False,
                                           help="Extract gene model fragments for e.g., PCR primer design")
    gpcrparser.add_argument("-i", dest="gff3File",
                            required=True,
                            help="Location of GFF3 file")
    gpcrparser.add_argument("-f", dest="fastaFile",
                            required=True,
                            help="Location of FASTA (e.g., genome) file")
    gpcrparser.add_argument("-m", dest="modelIdentifier",
                            required=True,
                            help="Identifier for GFF3 feature to model")
    gpcrparser.add_argument("-o", dest="outputFileName",
                            required=True,
                            help="Location to write model output")
    gpcrparser.add_argument("--buffer", dest="buffer",
                            required=False,
                            type=int,
                            help="Optionally obtain the provided length of surrounding sequence")
    
    # GFF3 > relabel mode
    grelabelparser = subGFF3Parsers.add_parser("relabel",
                                           parents=[p],
                                           add_help=False,
                                           help="Relabel parts of a GFF3 file")
    grelabelparser.add_argument("-i", dest="gff3File",
                                required=True,
                                help="Location of GFF3 file")
    grelabelparser.add_argument("-o", dest="outputFileName",
                                required=True,
                                help="Location to write modified GFF3")
    grelabelparser.add_argument("--list", dest="listFile",
                                required=False,
                                help="""Optionally, specify the location of a 2-column headerless
                                list with oldid:newid values for naive in-place line replacement""")
    grelabelparser.add_argument("--csuffix", dest="contigSuffix",
                                required=False,
                                help="""Optionally, specify a suffix to add to every contig e.g.,
                                if you have renamed the underlying genomic sequences""")
    grelabelparser.add_argument("--gsuffix", dest="geneSuffix",
                                required=False,
                                help="""Optionally, specify a suffix to add to every gene e.g.,
                                if you want to merge haplotype-resolved annotations together""")
    
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
    
    # GFF3 > miniprot subparser
    gff3mpparser = subGFF3Parsers.add_parser("miniprot",
                                             parents=[p],
                                             add_help=False,
                                             help="Miniprot GFF3 handling")
    gff3mpparser.set_defaults(func=gmain)
    
    subGFF3MpParsers = gff3mpparser.add_subparsers(dest="gff3MiniprotMode",
                                                   required=True)
    
    # GFF3 > miniprot > reformat mode
    gmpreformatparser = subGFF3MpParsers.add_parser("reformat",
                                                    parents=[p],
                                                    add_help=False,
                                                    help="Add gene and exon features")
    gmpreformatparser.add_argument("-i", dest="gff3File",
                                   required=True,
                                   help="Location of miniprot GFF3 file")
    gmpreformatparser.add_argument("-o", dest="outputFileName",
                                   required=True,
                                   help="Location to write reformatted miniprot GFF3 file")
    
    # GFF3 > miniprot > resolve mode
    gmpresolveparser = subGFF3MpParsers.add_parser("resolve",
                                             parents=[p],
                                             add_help=False,
                                             help="Resolve overlapping miniprot annotations; also reformats")
    gmpresolveparser.add_argument("-i", dest="gff3File",
                                  required=True,
                                  help="Location of miniprot GFF3 file")
    gmpresolveparser.add_argument("-o", dest="outputFileName",
                                  required=True,
                                  help="Location to write overlap-resolved miniprot GFF3 file")
    
    # GFF3 > to subparser
    gff3toparser = subGFF3Parsers.add_parser("to",
                                             parents=[p],
                                             add_help=False,
                                             help="GFF3 conversion")
    gff3toparser.set_defaults(func=gmain)
    
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
                              help="""Optionally, specify the character(s) used to denote a lack of 'map' value;
                              default = '_'""",
                              default="_")
    gtotsvparser.add_argument("--sep", dest="sepChar",
                              required=False,
                              help="""Optionally, specify the character(s) used to separate multiple values of
                              the same 'map' key; default = ';'""",
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
    
    # Homologs subparser
    hparser = subparsers.add_parser("homologs",
                                    parents=[p],
                                    add_help=False,
                                    help="Homologs TSV file handling")
    hparser.set_defaults(func=pmain)
    
    subHomologsParsers = hparser.add_subparsers(dest="homologsMode",
                                            required=True)
    
    # Homologs > annotate mode
    hannotateparser = subHomologsParsers.add_parser("annotate",
                                             parents=[p],
                                             add_help=False,
                                             help="Annotate a homologs file with GFF3 details")
    hannotateparser.add_argument("-i", dest="homologsFile",
                                 required=True,
                                 help="Location of homologs file")
    hannotateparser.add_argument("-g1", dest="gff3File1",
                                 required=True,
                                 help="Location of GFF3 file (first column of homologs)")
    hannotateparser.add_argument("-g2", dest="gff3File2",
                                 required=True,
                                 help="Location of GFF3 file (second column of homologs)")
    hannotateparser.add_argument("-o", dest="outputFileName",
                                 required=True,
                                 help="Location to write statistics output")
    
    # Homologs > to subparser
    homologstoparser = subHomologsParsers.add_parser("to",
                                                     parents=[p],
                                                     add_help=False,
                                                     help="Homologs conversion")
    homologstoparser.set_defaults(func=pmain)
    
    subHomologsToParsers = homologstoparser.add_subparsers(dest="homologsToMode",
                                                           required=True)
    
    # Homologs > to > BEDPE mode
    htobedpeparser = subHomologsToParsers.add_parser("bedpe",
                                             parents=[p],
                                             add_help=False,
                                             help="Homologs to BEDPE conversion")
    htobedpeparser.add_argument("-i", dest="homologsFile",
                                required=True,
                                help="Location of homologs file")
    htobedpeparser.add_argument("-g1", dest="gff3File1",
                                required=True,
                                help="Location of GFF3 file (first column of homologs)")
    htobedpeparser.add_argument("-g2", dest="gff3File2",
                                required=True,
                                help="Location of GFF3 file (second column of homologs)")
    htobedpeparser.add_argument("-o", dest="outputFileName",
                                required=True,
                                help="Location to write BEDPE output")
    
    # IRF mode
    irfparser = subparsers.add_parser("irf",
                                      parents=[p],
                                      add_help=False,
                                      help="Reformat IRF to GFF3")
    irfparser.set_defaults(func=irfmain)
    irfparser.add_argument("-i", dest="irfDatFile",
                           required=True,
                           help="Location of IRF .dat file")
    irfparser.add_argument("-o", dest="outputFileName",
                           required=True,
                           help="Location to write GFF3 output")
    irfparser.add_argument("--minLen", dest="minimumLength",
                           required=False,
                           type=int,
                           help="Optionally filter IRs < this length; set to 0 or -1 for no filtering",
                           default=0)
    irfparser.add_argument("--maxLen", dest="maximumLength",
                           required=False,
                           type=int,
                           help="Optionally filter IRs > this length; set to 0 or -1 for no filtering",
                           default=0)
    irfparser.add_argument("--minGap", dest="minimumGap",
                           required=False,
                           type=int,
                           help="Optionally filter IRs < this length; set to 0 or -1 for no filtering",
                           default=0)
    irfparser.add_argument("--maxGap", dest="maximumGap",
                           required=False,
                           type=int,
                           help="Optionally filter IRs > this length; set to 0 or -1 for no filtering",
                           default=0)
    irfparser.add_argument("--identity", dest="identityCutoff",
                           required=False,
                           type=float,
                           help="""Optionally filter IRs < this percentage identity (0-100); 
                           set to 0 for no filtering""",
                           default=0)
    
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
    if args.mode == "blast":
        print("## annotarium.py - BLAST handling ##")
        validate_b(args)
        bmain(args)
    if args.mode == "domains":
        print("## annotarium.py - Domain prediction handling ##")
        validate_d(args)
        dmain(args)
    if args.mode == "fasta":
        print("## annotarium.py - FASTA handling ##")
        validate_f(args)
        fmain(args)
    if args.mode == "gff3":
        print("## annotarium.py - GFF3 handling ##")
        validate_g(args)
        gmain(args)
    if args.mode == "irf":
        print("## annotarium.py - IRF results handling ##")
        validate_irf(args)
        irfmain(args)
    if args.mode == "homologs":
        print("## annotarium.py - Homologs (2 column IDs pair) handling ##")
        validate_p(args)
        pmain(args)
    if args.mode == "rnammer":
        print("## annotarium.py - RNAmmer handling ##")
        validate_rnammer(args)
        rnammermain(args)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def bmain(args):
    # Split into sub-mode-specific functions
    if args.blastMode == "to":
        validate_b_to(args)
        if args.blastToMode == "homologs":
            print("## Blast to homologs (reciprocal best 2 column) conversion ##")
            validate_b_to_homologs(args)
            blast_to_homologs(args)
    
    print("BLAST handling complete!")

def dmain(args):
    # Split into sub-mode-specific functions
    if args.domainsMode == "resolve":
        print("## Domain overlap resolution ##")
        validate_d_resolve(args)
        domains_resolve(args)
    
    print("Domain handling complete!")

def fmain(args):
    # Split into sub-mode-specific functions
    if args.fastaMode == "softmask":
        print("## FASTA softmask to BED tabulation ##")
        validate_f_softmask(args)
        fasta_softmask_to_bed(args)
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
    if args.gff3Mode == "pcr":
        print("## GFF3 PCR model ##")
        validate_g_pcr(args)
        gff3_pcr(args)
    if args.gff3Mode == "relabel":
        print("## GFF3 relabelling ##")
        validate_g_relabel(args)
        gff3_relabel(args)
    
    if args.gff3Mode == "miniprot":
        validate_g_mp(args)
        if args.gff3MiniprotMode == "reformat":
            print("## Reformat miniprot into a proper GFF3 ##")
            validate_g_mp_reformat(args)
            gff3_mp_reformat(args)
        if args.gff3MiniprotMode == "resolve":
            print("## Resolve overlapping miniprot annotations ##")
            validate_g_mp_resolve(args)
            gff3_mp_resolve(args)
    
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

def pmain(args):
    # Split into sub-mode-specific functions
    if args.homologsMode == "annotate":
        print("## Homologs annotation with GFF3 details ##")
        validate_p_annotate(args)
        homologs_annotate(args)
    if args.homologsMode == "to":
        validate_p_to(args)
        if args.homologsToMode == "bedpe":
            print("## Homologs to BEDPE conversion ##")
            validate_p_to_bedpe(args)
            homologs_to_bedpe(args)
    
    print("Homologs handling complete!")

def irfmain(args):
    irf_to_gff3(args)
    
    print("IRF handling complete!")

def rnammermain(args):
    rnammer_reformat(args)
    
    print("RNAmmer handling complete!")

if __name__ == "__main__":
    main()
