#! python3

import os, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, GzCapableWriter
from gff3 import GFF3Tarium

class HomologsParser:
    '''
    The HomologsParser Class provides parsing capability for homologs files
    (2 columns: seqID1 seqID2).
    
    Initialisation:
        fileLocation -- a string indicating the location of the outfmt6 file that was
                        parsed to generate this object instance
    '''
    def __init__(self, fileLocation):
        self.fileLocation = fileLocation
        self.left = {}
        self.right = {}
        self._parse()
        
        # Also set helper attribute
        self.isHomologsParser = True
    
    @property
    def fileLocation(self):
        return self._fileLocation
    
    @fileLocation.setter
    def fileLocation(self, value):
        if not isinstance(value, str):
            raise TypeError(f"fileLocation should be 'str', not '{type(value).__name__}'")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"fileLocation '{value}' does not point to an existing file location")
        
        self._fileLocation = value
    
    def _parse(self):
        '''
        Sets:
            self.left / self.right -- dictionary with keys equal to sequence ID (of its column) and value equal to
                                      corresponding paralog sequence ID
        '''
        self.left = {}
        self.right = {}
        
        with read_gz_file(self.fileLocation) as fileIn:
            for line in fileIn:
                leftID, rightID = line.rstrip("\r\n\t ").split("\t")
                
                if leftID in self.left:
                    raise ValueError(f"'{leftID}' already encountered in left column")
                else:
                    self.left[leftID] = rightID
                
                if rightID in self.right:
                    raise ValueError(f"'{rightID}' already encountered in right column")
                else:
                    self.right[rightID] = leftID
    
    def __iter__(self):
        for leftID, rightID in self.left.items():
            yield leftID, rightID
    
    def __len__(self):
        return len(self.left)
    
    def __contains__(self, value):
        return True if (value in self.left or value in self.right) else False
    
    def __str__(self):
        return (f"HomologsParser parsed '{self.fileLocation}' and found {len(self)} paralog pairs")
    
    def __repr__(self):
        return "<HomologsParser object; fileLocation='{0}'>".format(
            self.fileLocation
        )

def homologs_annotate(args):
    # Parse homologs file
    homologs = HomologsParser(args.homologsFile)
    
    # Parse GFF3 files 1 and 2
    gff3_1 = GFF3Tarium(args.gff3File1)
    gff3_2 = GFF3Tarium(args.gff3File2)
    
    # Format paralog details
    details = []
    for leftID, rightID in homologs:
        leftFeature = gff3_1[leftID]
        rightFeature = gff3_2[rightID]
        
        leftDetails = [leftID, leftFeature.contig, leftFeature.start, leftFeature.end]
        rightDetails = [rightID, rightFeature.contig, rightFeature.start, rightFeature.end]
        
        details.append([leftDetails, rightDetails])
    details.sort(key = lambda x: (x[0][1], x[0][2])) # sort by contig, start of left feature details
    
    # Format output file
    with GzCapableWriter(args.outputFileName) as fileOut:
        for leftDetails, rightDetails in details:
            leftDetails = "\t".join(map(str, leftDetails))
            rightDetails = "\t".join(map(str, rightDetails))
            fileOut.write(f"{leftDetails}\t{rightDetails}\n")

def homologs_to_bedpe(args):
    # Parse homologs file
    homologs = HomologsParser(args.homologsFile)
    
    # Parse GFF3 files 1 and 2
    gff3_1 = GFF3Tarium(args.gff3File1)
    gff3_2 = GFF3Tarium(args.gff3File2)
    
    # Format paralog details in BEDPE format
    details = []
    for leftID, rightID in homologs:
        leftFeature = gff3_1[leftID]
        rightFeature = gff3_2[rightID]
        
        leftDetails = [leftFeature.contig, leftFeature.start, leftFeature.end]
        rightDetails = [rightFeature.contig, rightFeature.start, rightFeature.end]
        
        details.append([leftDetails, rightDetails])
    details.sort(key = lambda x: (x[0][0], x[0][1])) # sort by contig, start of left feature details
    
    # Format output file
    with GzCapableWriter(args.outputFileName) as fileOut:
        for leftDetails, rightDetails in details:
            leftDetails = "\t".join(map(str, leftDetails))
            rightDetails = "\t".join(map(str, rightDetails))
            fileOut.write(f"{leftDetails}\t{rightDetails}\n")
