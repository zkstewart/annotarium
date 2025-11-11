#! python3

import os, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, GzCapableWriter

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
    def __init__(self, fileLocation, evalue=None, numHits=None):
        self.fileLocation = fileLocation
        self.evalue = evalue
        self.numHits = numHits
        self.results = {}
        
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
                
                # Skip if evalue isn't significant
                if self.evalue != None and evalue > self.evalue: # self.evalue might differ between BLAST run and parsing
                    continue
                
                # Store result
                if indexQuery:
                    blastDict.setdefault(qid, [])
                    blastDict[qid].append({
                        "qid": qid, "tid": tid, "identity": identity,
                        "length": length, "mismatch": mismatch, "gapopen": gapopen,
                        "qstart": qstart, "qend": qend, "tstart": tstart,
                        "tend": tend, "evalue": evalue, "bitscore": bitscore
                    })
                if indexTarget:
                    blastDict.setdefault(tid, [])
                    blastDict[tid].append({
                        "qid": qid, "tid": tid, "identity": identity,
                        "length": length, "mismatch": mismatch, "gapopen": gapopen,
                        "qstart": qstart, "qend": qend, "tstart": tstart,
                        "tend": tend, "evalue": evalue, "bitscore": bitscore
                    })
        
        # Sort individual entries in blastDict
        for value in blastDict.values():
            value.sort(key = lambda x: (-x["bitscore"], x["evalue"])) # sort by bitscore (higher) and evalue (lower)
        
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

def blast_to_homologs(args):
    # Parse files 1 and 2
    results1 = Outfmt6Parser(args.inputFile1, evalue=args.evalue, numHits=1)
    results1.parse()
    results2 = Outfmt6Parser(args.inputFile2, evalue=args.evalue, numHits=1)
    results2.parse()
    
    # Check for reciprocity
    reciprocalHits = []
    for qid, hits in results1.items():
        tid = hits[0]["tid"]
        if tid in results2:
            tidMatch = results2[tid][0]["tid"] # the tid of the tid should be the qid if files are reciprocal
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
