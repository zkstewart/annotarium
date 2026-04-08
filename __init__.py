import os, sys

# Make classes accessible to variantopia
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.domains import Domains
from modules.coordinates import OverlapResolver
from modules.gff3 import GFF3Feature, GFF3Tarium
from modules.fasta import TranslationTable, Sequence, Records, FASTATarium
