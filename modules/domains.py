#! python3

import os, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from coordinates import OverlapResolver

class DomainFeature:
    '''
    self.score points to self.evalue, and is given to conform to the interface expected
    by OverlapResolver.
    '''
    def __init__(self, id, start, end, evalue):
        self.id = id
        self.start = start
        self.end = end
        self.evalue = evalue
        
        # Set helper attributes
        self.isDomainFeature = True
    
    @property
    def evalue(self):
        return self._evalue
    
    @evalue.setter
    def evalue(self, value):
        if not isinstance(value, float) and not isinstance(value, int):
            raise TypeError(f"evalue must be a 'float' or 'int', not '{type(value).__name__}'")
        
        if value < 0:
            raise ValueError("evalue must be a value >= 0")
        
        self._evalue = value
    
    @property
    def score(self):
        return self.evalue
    
    @score.setter
    def score(self, value):
        self.evalue = value
    
    @property
    def start(self):
        return self._start
    
    @start.setter
    def start(self, value):
        if not isinstance(value, int):
            raise TypeError(f"start must be an 'int', not '{type(value).__name__}'")
        
        if value < 0:
            raise ValueError("start must be a value >= 0")
        
        self._start = value
    
    @property
    def end(self):
        return self._end
    
    @end.setter
    def end(self, value):
        if not isinstance(value, int):
            raise TypeError(f"end must be an 'int', not '{type(value).__name__}'")
        
        if value < 0:
            raise ValueError("end must be a value >= 0")
        
        self._end = value
    
    def to_list(self):
        return [self.id, self.start, self.end, self.evalue]
    
    def __repr__(self):
        return "<DomainFeature object; id={0}; start={1}, end={2}, evalue={3}>".format(
            self.id, self.start, self.end, self.evalue
        )

class DomainsFormatError(Exception):
    pass

class Domains:
    '''
    Attributes:
        domDict -- the main data structure for recording domain prediction results
                   with structure like:
                   {
                        "proteinID1": [
                            ["domainID1", start, end, evalue],
                            ["domainID2", start, end, evalue],
                            ...
                        ],
                        "proteinID2": [ [...] ]
                    }
        evalue -- a float indicating the evalue to use for filtering any domain predictions
                  to omit during file parsing
    '''
    def __init__(self, domDict=None, evalue=None):
        self.domDict = domDict
        self.evalue = evalue
        
        # Set helper attributes
        self.isDomains = True
    
    @property
    def domDict(self):
        return self._domDict
    
    @domDict.setter
    def domDict(self, value):
        if value is None:
            self._domDict = None
            return
        
        if not isinstance(value, dict):
            raise TypeError(f"domDict must be a 'dict', not '{type(value).__name__}'")
        
        self._domDict = value
    
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
    
    def parse_domtblout(self, fileLocation):
        '''
        Parses the domtblout file from hmmsearch with optional evalue filtering.
        
        Sets:
            self.domDict -- a dict with structure like:
                            {
                                "proteinID1": [
                                    ["domainID1", start, end, evalue],
                                    ["domainID2", start, end, evalue],
                                    ...
                                ],
                                "proteinID2": [ [...] ]
                            }
        '''
        domDict = {}
        with open(fileLocation, "r") as fileIn:
            for line in fileIn:
                l = line.rstrip("\r\n\t'\" ")
                
                # Skip unnecessary lines
                if line.startswith("#") or line == "":
                    continue
                
                # Parse line for relevant details
                sl = l.split()
                try:
                    pid = sl[0]
                    did = sl[3]
                    evalue = float(sl[12])
                    dstart = int(sl[17])
                    dend = int(sl[18])
                except:
                    raise DomainsFormatError(f"'{fileLocation}' is not a valid .domtblout file; line which caused " +
                                             f"this error is '{l}'")
                
                # Skip if evalue is not significant
                if self.evalue and evalue > self.evalue:
                    continue
                
                # Format a DomainFeature object
                feature = DomainFeature(did, dstart, dend, evalue)
                
                # Add into feature dictionary
                if pid not in domDict:
                    domDict[pid] = [feature]
                else:
                    domDict[pid].append(feature)
        self.domDict = domDict
    
    def parse_parsed_domtblout(self, fileLocation):
        '''
        Parses a 'parsed domtblout' file as produced by .write_parsed_domtblout()
        with optional evalue filtering.
        
        Parameters:
            fileLocation -- a string indicating the location of the 'parsed domtblout' file
                            to read in
        Sets:
            self.domDict -- a dict with structure like:
                            {
                                "proteinID1": [
                                    ["domainID1", start, end, evalue],
                                    ["domainID2", start, end, evalue],
                                    ...
                                ],
                                "proteinID2": [ [...] ]
                            }
        '''
        domDict = {}
        with open(fileLocation, "r") as fileIn:
            for line in fileIn:
                l = line.rstrip("\r\n\t'\" ")
                
                # Skip unnecessary lines
                if line.startswith("#") or line == "":
                    continue
                
                # Parse line
                sl = l.split("\t")
                try:
                    pid, domainPredictions = sl[0], sl[1:]
                except:
                    raise DomainsFormatError(f"'{fileLocation}' is not a valid .domains.tsv file; line which caused " +
                                             f"this error is '{l}'")
                
                # Iterate over parsed domain predictions
                for domain in domainPredictions:
                    try:
                        did, dstart, dend, evalue = eval(domain)
                    except:
                        raise DomainsFormatError(f"'{fileLocation}' is not a valid .domains.tsv file; line which caused " +
                                                f"this error is '{l}'")
                    
                    # Skip domain if evalue is not significant
                    if self.evalue and evalue > self.evalue:
                        continue
                    
                    # Format a DomainFeature object
                    feature = DomainFeature(did, dstart, dend, evalue)
                    
                    # Add into feature dictionary
                    if pid not in domDict:
                        domDict[pid] = [feature]
                    else:
                        domDict[pid].append(feature)
        self.domDict = domDict
    
    def resolve_overlaps(self, ovlCutoff=0.25):
        '''
        Calls the Coordinates.OverlapResolver class to resolve overlaps within
        self.domDict. Replaces self.domDict with the resolved version.
        '''
        resolver = OverlapResolver(ovlCutoff=ovlCutoff)
        resolvedDomDict = {}
        for pid, featureList in self.domDict.items():
            resolvedDomDict[pid] = resolver.resolve(featureList)
        self.domDict = resolvedDomDict
    
    def write_parsed_domtblout(self, outputFileLocation):
        '''
        Writes the domain predictions stored in self.domDict to file as a
        'parsed domtblout' format file. This can be read back in by
        parse_parsed_domtblout().
        
        Parameters:
            outputFileLocation -- a string indicating the location to write
                                  the output file; must not already exist
        '''
        if os.path.exists(outputFileLocation):
            raise FileExistsError(f".write_parsed_domtblout() will not overwrite '{outputFileLocation}'")
        
        with open(outputFileLocation, "w") as fileOut:
            for pid, domainFeatures in self.domDict.items():
                formattedPredictions = "\t".join(map(str, [ feature.to_list() for feature in domainFeatures ]))
                fileOut.write(f"{pid}\t{formattedPredictions}\n")
    
    def __getitem__(self, key):
        return self.domDict[key]
    
    def __len__(self):
        return len(self.domDict)
    
    def __iter__(self):
        return self.domDict.items()
    
    def __contains__(self, value):
        return value in self.domDict
    
    def keys(self):
        return self.domDict.keys()
    
    def values(self):
        return self.domDict.values()
    
    def items(self):
        return self.domDict.items()
    
    def __repr__(self):
        if self.domDict:
            return "<Domains object; hasDomDict=True; numProteins={0}, evalue={1}>".format(
                len(self.domDict), self.evalue
            )
        else:
            return "<Domains object; hasDomDict=False, evalue={0}>".format(
                self.evalue
            )

def domains_resolve(args):
    # Parse the domains file
    domains = Domains(evalue=args.evalue)
    if args.isDomainsTsvFormat:
        domains.parse_parsed_domtblout(args.domainsFile)
    else:
        domains.parse_domtblout(args.domainsFile)
    
    # Resolve overlaps
    domains.resolve_overlaps(ovlCutoff=args.overlapCutoff)
    
    # Write output
    domains.write_parsed_domtblout(args.outputFileName)
