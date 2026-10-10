#! python3

import codecs
import gzip
import math
import re
from Bio import SeqIO
from contextlib import contextmanager

def get_codec(fileName):
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            break
        f.close()
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                break
            f.close()
            return "utf-16"
        except UnicodeDecodeError:
            print(f"'{fileName}' is neither utf-8 nor utf-16 encoded; please convert to one of these formats.")

@contextmanager
def read_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename, "r", encoding=get_codec(filename)) as f:
            yield f

class GzCapableWriter:
    def __init__(self, filename):
        self.filename = filename
        self.file = None
    
    def __enter__(self):
        if self.filename is None:
            return None
        else:
            if self.filename.endswith(".gz"):
                self.file = gzip.open(self.filename, "wt")
            else:
                self.file = open(self.filename, "w")
            return self.file
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file:
            self.file.close()
        if exc_type is not None:
            raise exc_type(exc_val).with_traceback(exc_tb)

class Emptyfile(object):
    def write(self, data):
        pass # ignore the data
    def __enter__(self): return self
    def __exit__(*x): pass

@contextmanager
def write_conditionally(fileName):
    if fileName == None:
        empty = Emptyfile()
        yield empty
    else:
        if fileName.endswith(".gz"):
            with gzip.open(fileName, "wt") as f:
                yield f
        else:
            with open(fileName, "w") as f:
                yield f

def parse_annotation_table(fileName, delimiter="\t"):
    '''
    Returns:
        dataDict -- a dictionary where keys are column headers and
                    values are the text contents of this row's value.
                    An additional "leftcolumn" key is added with the
                    first column's value which should correspond to
                    the ID.
    '''
    alreadyWarned = False
    header = None
    with open(fileName, "r") as fileIn:
        for line in fileIn:
            sl = line.strip().split(delimiter)
            if header == None:
                header = sl
                headerIndex = { h:i for i,h in enumerate(header) }
                continue
            
            dataDict = {"leftcolumn": sl[0]}
            try:
                dataDict.update({ h:sl[i] for h, i in headerIndex.items() })
            except:
                dataDict.update({ h:"." for h, i in headerIndex.items() })
                if not alreadyWarned:
                    print(f"WARNING: line '{sl}' does not have '{len(headerIndex)}' columns as expected;" + 
                          " this and similar gene lines will have '.' values imputed for all columns")
                    alreadyWarned = True
            yield dataDict

def parse_list_file(fileName):
    '''
    Returns:
        listedValues -- a list where each value is the string content of a line within
                        the provided file which has had .strip() applied to it.
    '''
    listedValues = []
    with read_gz_file(fileName) as fileIn:
        for line in fileIn:
            l = line.strip()
            if l != "":
                listedValues.append(l)
    return listedValues

def parse_fasta_to_lengths(fastaFile):
    '''
    Returns:
        lengthsDict -- a dictionary where keys are contig IDs and values are integers giving
                       the contig's length
    '''
    with read_gz_file(fastaFile) as fileIn:
        genomeRecords = SeqIO.parse(fileIn, "fasta")
        lengthsDict = { record.id:len(record) for record in genomeRecords }
    
    if lengthsDict == {}:
        raise ValueError(f"No contigs found in genome FASTA '{fastaFile}'; is it actually a FASTA file?")
    else:
        return lengthsDict

def parse_fasta_regions(regions, lengthsDict, argName="--regions"):
    '''
    Code helpfully copied from psQTL with minor modification.
    
    Parameters:
        regions -- a list of strings in 'contig:start-end' format, and/or
                   in 'contig' format; end can be > start to indicate e.g.,
                   reverse complementation
    Returns:
        parsedRegions -- a list of dictionaries with structure like:
                         [
                             {
                                 "contig": contigID, # string
                                 "start": start, # int
                                 "end": end, # int
                                 "reverse": reverse] # bool
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
            
            # Validate contig ID
            if not contigID in lengthsDict:
                raise ValueError(f"{argName} contig ID '{contigID}' not found in the -f FASTA!")
            
            # Validate start position
            if start < 0:
                raise ValueError(f"{argName} start position '{start}' is < 0!")
            if start == 0:
                print(f"# Note: {argName} start position '{start}' was set to 1, to use 1-based indexing.")
                start = 1
            if start == end:
                raise ValueError(f"{argName} start position '{start}' is equal to end position '{end}'!")
            
            # Detect reverse orientation and swap start/end if necessary
            reverse = False
            if start > end:
                start, end = end, start
                reverse = True
            
            # Validate end position
            if end > lengthsDict[contigID]:
                raise ValueError(f"{argName} '{contigID, start, end}' end position is > contig length '{lengthsDict[contigID]}'!")
            
            # Store region
            if start <= 1 and end == lengthsDict[contigID]:
                # If the region covers the full contig, set 'full' to True
                parsedRegions.append({"contig": contigID, "start": start, "end": end, "reverse": reverse, "full": True})
            else:
                # Otherwise, set 'full' to False
                parsedRegions.append({"contig": contigID, "start": start, "end": end, "reverse": reverse, "full": False})
        
        # Handle invalid format
        elif ":" in region:
            raise ValueError(f"Invalid region input '{region}'; you included a ':' but did " + 
                             "not format the region as 'chr:start-end'!")
        # Handle chr format
        else:
            if not region in lengthsDict:
                raise ValueError(f"{argName} contig ID '{region}' not found in the -f FASTA!")
            parsedRegions.append({"contig": region, "start": 1, "end": lengthsDict[region], "reverse": False, "full": True})
    
    # Handle empty regions
    "Empty is interpreted as all regions"
    if parsedRegions == []:
        parsedRegions = [
            {"contig": contigID, "start": 1, "end": lengthsDict[contigID], "reverse": False, "full": True}
            for contigID in lengthsDict
        ]
    
    return parsedRegions

def parse_gff3_regions(regions):
    '''
    Simplified regions parsing without needing to know reverse complementation
    or overall contig length
    
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

def get_chunking_points(numberToChunk, chunks, isNumOfChunks=True):
    '''
    Code helpfully copied from variantopia with minor modification.
    
    This is a general purpose function to take in a number of "things"
    that you want to chunk, and find out how to chunk them evenly.
    
    The resulting list should be interpreted as the 0-based indices where
    a new chunk should form. You should check for this index at the start
    of a loop, and form a new file if your index == the value in this list.
    
    Also, this uses "allocated chunking" such that it will try to keep
    the number of things per chunk approximately equal. Even if you specify
    X number of things per chunk, it might be more optimal to have X-1 in each
    chunk so as to make sure the last chunk doesn't contain a single thing.
    
    Parameters:
        numberToChunk -- an integer value, possibly derived from a list length as example.
        chunks -- an integer value for the desired number of chunks OR the number of
                  things to contain within each chunk, determined by
        isNumOfChunks -- a boolean indicating whether you want the number to be the number
                         of chunks (True), or the number of things within each chunk (False)
    Returns:
        chunkPoints -- a list of integers indicating the 0-based indices where a new chunk
                       should begin; the first chunk (at index 0) is implicit and not included
    '''
    assert isinstance(numberToChunk, int)
    assert isinstance(chunks, int)
    if isNumOfChunks and numberToChunk < chunks:
        raise ValueError(f"Chunking only valid if chunks <= numberToChunk i.e., {chunks} <= {numberToChunk}")
    
    # Derive how many chunks we want to split the things into
    if isNumOfChunks:
        numChunks = chunks
    else:
        numChunks = math.ceil(numberToChunk / chunks)
    
    # Find how many things go into each chunk, and how many chunks get one extra thing
    baseSize, numRoundedUp = divmod(numberToChunk, numChunks)
    
    # Store positions at which to start a new chunk
    chunkPoints = []
    ongoingCount = 0
    for i in range(numChunks - 1): # the last chunk runs until the end
        ongoingCount += baseSize + 1 if i < numRoundedUp else baseSize
        chunkPoints.append(ongoingCount)
    
    return chunkPoints
