#! python3

import os, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, GzCapableWriter

class SearchResult:
    '''
    The SearchResult class contains information of a BLAST/MMseqs2/equivalent result for
    the alignment of two sequences. It expects information akin to what would be found in
    outfmt6, but tolerates None values to allow for flexible usage.
    
    Minimal validation is done in the interest of fast parsing operation. Make sure to use this
    class properly as it will not hold your hand.
    
    Initialisation:
        queryID -- a string identifying the query sequence.
        targetID -- a string identifying the target sequence.
        identity -- a float indicating the percentage of similarity between the two sequences; typically
                    identity is calculated by counting the number of 'columns' that are identical
                    divided by the total length of the alignment.
        alignedLength -- an int counting the length of the alignment start to end.
        mismatch -- an int counting the number of substitutions in the alignment.
        gapopen -- an int counting the number of gaps in the alignment.
        queryStart -- an int (1-based) for the position where alignment begins on the query sequence.
        queryEnd -- an int (1-based, inclusive) for the position where alignment ends on the query sequence.
        targetStart -- an int (1-based) for the position where alignment begins on the target sequence.
        targetEnd -- an int (1-based, inclusive) for the position where alignment ends on the target sequence.
        evalue -- a float value for the 'expect value'.
        score -- a float or int numeric value; with BLAST this would be the bitscore.
        null -- any value to return in lieu of None if something in this class would return or utilise
                an attribute that is set to None; default == None.
    '''
    def __init__(self, queryID, targetID, identity=None, alignedLength=None, mismatch=None,
                 gapopen=None, queryStart=None, queryEnd=None, targetStart=None, targetEnd=None,
                 evalue=None, score=None, null=None):
        self.queryID = queryID
        self.targetID = targetID
        self.identity = identity
        self.alignedLength = alignedLength
        self.mismatch = mismatch
        self.gapopen = gapopen
        self.queryStart = queryStart
        self.queryEnd = queryEnd
        self.targetStart = targetStart
        self.targetEnd = targetEnd
        self.evalue = evalue
        self.score = score
        self.null = null
    
    @property
    def queryID(self):
        if self._queryID == None:
            return self.null
        return self._queryID
    
    @queryID.setter
    def queryID(self, value):
        self._queryID = value
    
    @property
    def targetID(self):
        if self._targetID == None:
            return self.null
        return self._targetID
    
    @targetID.setter
    def targetID(self, value):
        self._targetID = value
    
    @property
    def identity(self):
        if self._identity == None:
            return self.null
        return self._identity
    
    @identity.setter
    def identity(self, value):
        self._identity = value
    
    @property
    def alignedLength(self):
        if self._alignedLength == None:
            return self.null
        return self._alignedLength
    
    @alignedLength.setter
    def alignedLength(self, value):
        self._alignedLength = value
    
    @property
    def mismatch(self):
        if self._mismatch == None:
            return self.null
        return self._mismatch
    
    @mismatch.setter
    def mismatch(self, value):
        self._mismatch = value
    
    @property
    def gapopen(self):
        if self._gapopen == None:
            return self.null
        return self._gapopen
    
    @gapopen.setter
    def gapopen(self, value):
        self._gapopen = value
    
    @property
    def queryStart(self):
        if self._queryStart == None:
            return self.null
        return self._queryStart
    
    @queryStart.setter
    def queryStart(self, value):
        self._queryStart = value
    
    @property
    def queryEnd(self):
        if self._queryEnd == None:
            return self.null
        return self._queryEnd
    
    @queryEnd.setter
    def queryEnd(self, value):
        self._queryEnd = value
    
    @property
    def targetStart(self):
        if self._targetStart == None:
            return self.null
        return self._targetStart
    
    @targetStart.setter
    def targetStart(self, value):
        self._targetStart = value
    
    @property
    def targetEnd(self):
        if self._targetEnd == None:
            return self.null
        return self._targetEnd
    
    @targetEnd.setter
    def targetEnd(self, value):
        self._targetEnd = value
    
    @property
    def evalue(self):
        if self._evalue == None:
            return self.null
        return self._evalue
    
    @evalue.setter
    def evalue(self, value):
        self._evalue = value
    
    @property
    def score(self):
        if self._score == None:
            return self.null
        return self._score
    
    @score.setter
    def score(self, value):
        self._score = value
    
    def to_outfmt6(self):
        '''
        Returns:
            outfmt6Str -- a string (without newline termination) in outfmt6 format.
        '''
        return (f"{self.queryID}\t{self.targetID}\t{self.identity}\t{self.alignedLength}\t" + 
                f"{self.mismatch}\t{self.gapopen}{self.queryStart}\t{self.queryEnd}\t" + 
                f"{self.targetStart}\t{self.targetEnd}\t{self.evalue}\t{self.score}")
    
    def __str__(self):
        return self.to_outfmt6()
    
    def __repr__(self):
        return "<SearchResult object; queryID='{0}', targetID={1}>".format(
            self.queryID, self.targetID
        )

class NotParsedError(Exception):
    pass

class Outfmt6Parser:
    '''
    The BlastParser Class provides parsing capability for outfmt6 files, regardless
    of whether they were made by BLAST, MMSeqs2, or otherwise.
    
    Initialisation:
        fileLocation -- a string indicating the location of the outfmt6 file that was
                        parsed to generate this object instance
        evalue -- (OPTIONAL) an int or float value indicating the cutoff to enforce for
                  retrieving hits for a query OR None for no cutoff
        numHits -- (OPTIONAL) an int value to limit retention of only the best numHits
                   results per query OR None for no limit (i.e., retain all hits)
    '''
    def __init__(self, fileLocation, evalue=None, numHits=None, identity=None):
        self.fileLocation = fileLocation
        self.evalue = evalue
        self.numHits = numHits
        self.identity = identity
        self.results = None
        
        # Also set helper attribute
        self.isOutfmt6Parser = True
    
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
    
    @property
    def evalue(self):
        return self._evalue
    
    @evalue.setter
    def evalue(self, value):
        if value is None:
            self._evalue = None
            return
        
        if not isinstance(value, float) or isinstance(value, int):
            raise TypeError(f"evalue must be a 'float' or 'int', not '{type(value).__name__}'")
        
        if value < 0:
            raise ValueError("evalue must be a value >= 0")
        
        self._evalue = value
    
    @property
    def numHits(self):
        return self._numHits
    
    @numHits.setter
    def numHits(self, value):
        if value is None:
            self._numHits = None
            return
        
        if not isinstance(value, int):
            raise TypeError(f"numHits must be an 'int', not '{type(value).__name__}'")
        
        if value <= 0:
            raise ValueError("numHits cannot be 0 or negative (which would retrieve no hits)")
        
        self._numHits = value
    
    @property
    def identity(self):
        return self._identity
    
    @identity.setter
    def identity(self, value):
        if value is None:
            self._identity = None
            return
        
        if not isinstance(value, float) or isinstance(value, int):
            raise TypeError(f"identity must be a 'float' or 'int', not '{type(value).__name__}'")
        
        if value < 0:
            raise ValueError("identity must be a value >= 0")
        if value > 1:
            raise ValueError("identity must be a value <= 0")
        
        self._identity = value
    
    @property
    def results(self):
        if self._results == None:
            raise NotParsedError("Outfmt6Parser .results is None; call .parse() first!")
        return self._results
    
    @results.setter
    def results(self, value):
        if value is None:
            self._results = None
            return
        
        if not isinstance(value, dict):
            raise TypeError(f"results must be a 'dict', not '{type(value).__name__}'")
        self._results = value
    
    @identity.setter
    def identity(self, value):
        if value is None:
            self._identity = None
            return
        
        if not isinstance(value, float) or isinstance(value, int):
            raise TypeError(f"identity must be a 'float' or 'int', not '{type(value).__name__}'")
        
        if value < 0:
            raise ValueError("identity must be a value >= 0")
        if value > 1:
            raise ValueError("identity must be a value <= 0")
        
        self._identity = value
    
    def parse(self, indexQuery=True, indexTarget=False):
        '''
        Parses the given self.fileLocation outfmt6 file according to self.evalue
        and self.numHits. Requires one or both of indexQuery and/or indexTarget to be
        True. Results will be flawed if identical sequence IDs exist in both index and
        target files and both are indexed.
        
        Parameters:
            indexQuery -- (OPTIONAL) store results by query ID in the underlying dict
            indexTarget -- (OPTIONAL) store results by target ID in the underlying dict
        Sets:
            self.results -- a dict with structure like:
                                query_id1: [
                                    [target_id, identity_pct, query_start, query_end, target_start, target_end, evalue],
                                    ...
                                ],
                                query_id2: [ ... ],
                                ...
        '''
        if (not indexQuery) and (not indexTarget):
            raise ValueError("Cannot parse Outfmt6 unless one or both of indexQuery and indexTarget are set")
        
        blastDict = {}
        
        with read_gz_file(self.fileLocation) as fileIn:
            for line in fileIn:
                # Extract details
                qid, tid, identity, length, mismatch, gapopen, qstart, qend, \
                    tstart, tend, evalue, bitscore = line.rstrip("\r\n").split("\t")
                
                identity, evalue, bitscore = float(identity), float(evalue), float(bitscore)
                length, mismatch, gapopen, qstart, qend, tstart, tend = \
                    int(length), int(mismatch), int(gapopen), int(qstart), int(qend), int(tstart), int(tend)
                
                # Skip based on filtration values
                if self.evalue != None and evalue > self.evalue: # self.evalue might differ between BLAST run and parsing
                    continue
                if self.identity != None and identity < self.identity:
                    continue
                
                # Store result
                if indexQuery:
                    blastDict.setdefault(qid, [])
                    blastDict[qid].append(
                        SearchResult(
                            queryID=qid, targetID=tid,
                            identity=identity, alignedLength=length,
                            mismatch=mismatch, gapopen=gapopen,
                            queryStart=qstart, queryEnd=qend,
                            targetStart=tstart, targetEnd=tend,
                            evalue=evalue, score=bitscore
                        )
                    )
                if indexTarget:
                    blastDict.setdefault(tid, [])
                    blastDict[tid].append(
                        SearchResult(
                            queryID=qid, targetID=tid,
                            identity=identity, alignedLength=length,
                            mismatch=mismatch, gapopen=gapopen,
                            queryStart=qstart, queryEnd=qend,
                            targetStart=tstart, targetEnd=tend,
                            evalue=evalue, score=bitscore
                        )
                    )
        
        # Sort individual entries in blastDict
        for value in blastDict.values():
            value.sort(key = lambda x: (-x.score, x.evalue)) # sort by bitscore (higher) and evalue (lower)
        
        # Enforce numHits threshold
        if self.numHits != None:
            for key in blastDict.keys():
                blastDict[key] = blastDict[key][0:self.numHits]
        
        self.results = blastDict
    
    def keys(self):
        return self.results.keys()
    
    def items(self):
        return self.results.items()
    
    def __iter__(self):
        return self.results.items()
    
    def __len__(self):
        return len(self.results)
    
    def __getitem__(self, key):
        if key in self.results:
            return self.results[key]
    
    def __contains__(self, value):
        return True if value in self.results else False
    
    def __str__(self):
        return (f"Outfmt6Parser parsed '{self.fileLocation}' with evalue={self.evalue} cutoff and " + 
               f"max. {self.numHits} hits returned per query; for {len(self)} keys there are " + 
               "{0} hits".format(
                   sum([ len(v) for v in self.results.values() ])
               ))
    
    def __repr__(self):
        return "<Outfmt6Parser object; fileLocation='{0}', evalue={1}, numHits={2}>".format(
            self.fileLocation, self.evalue, self.numHits
        )

def blast_filter(args):
    results = Outfmt6Parser(args.outfmt6File, evalue=args.evalue, numHits=args.maxhits, identity=args.identity)
    results.parse()
    
    with GzCapableWriter(args.outputFileName) as fileOut:
        for qid, searchResultList in results.items():
            for searchResult in searchResultList:
                fileOut.write(searchResult.to_outfmt6() + "\n")

def blast_to_homologs(args):
    # Parse files 1 and 2
    results1 = Outfmt6Parser(args.inputFile1, evalue=args.evalue, numHits=1)
    results1.parse()
    results2 = Outfmt6Parser(args.inputFile2, evalue=args.evalue, numHits=1)
    results2.parse()
    
    # Check for reciprocity
    reciprocalHits = []
    for qid, searchResultList in results1.items():
        tid = searchResultList[0].targetID
        if tid in results2:
            tidMatch = results2[tid][0].targetID # the tid of the tid should be the qid if files are reciprocal
            if tidMatch == qid:
                reciprocalHits.append((qid, tid))
    
    # Output file
    with GzCapableWriter(args.outputFileName) as fileOut:
        for qid, tid in reciprocalHits:
            fileOut.write(f"{qid}\t{tid}\n")
    
    # Warn user if no hits were found
    if len(reciprocalHits) == 0:
        print("WARNING: No reciprocal best BLAST hits were found. " + 
              "Check your files and arguments to ensure this makes sense.")
