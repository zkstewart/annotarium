#! python3

import os, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, GzCapableWriter

class FormatError(Exception):
    pass

class Member:
    '''
    Simple container for a cluster member. 'id' attribute is required and any
    arbitrary values can be specified as a keyword arg.
    '''
    def __init__(self, id, **kwargs):
        self.id = id
        for key, value in kwargs.items():
            self.__dict__[key] = value
        
        # Set helper attributes
        self.isMember = True
    
    def __eq__(self, other):
        if not hasattr(other, "id"):
            return self.id == other
        return self.id == other.id
    
    def __lt__(self, other):
        if not hasattr(other, "id"):
            self.id < other
        return self.id < other.id
    
    def __hash__(self):
        return hash(self.id) # a member is defined ONLY by its id; other kwargs are just toppings
    
    def __repr__(self):
        kwKeys = [ key for key, value in self.__dict__.items() if not key in ["id", "isMember"] ]
        
        return "<Member object;id='{0}'{1}>".format(
            self.id,
            "" if len(kwKeys) == 0 else ";" + ";".join([ f"{key}={self.__dict__[key]}" for key in kwKeys ])
        )

class Cluster:
    def __init__(self, id, members=None):
        self.id = id
        self.members = members
        
        # Set helper attributes
        self.isCluster = True
    
    @property
    def members(self):
        return self._members
    
    @members.setter
    def members(self, value):
        '''
        This setter is intentionally not mandating the members be a set of Member objects, although
        that class does exist for the purpose of parsing a clustering file.
        '''
        if value == None:
            self._members = set()
            return
        
        if not isinstance(value, set):
            try:
                originalLength = len(value)
            except TypeError:
                raise TypeError(f"Cluster.members cannot handle type '{type(value).__name__}' as it has no len operator")
            
            try:
                value = set(value)
            except TypeError:
                raise TypeError("Cluster.members must be a 'set' or a hashable iterable")
            
            if len(value) != originalLength:
                raise ValueError(f"Cluster.members was given an iterable with non-unique values; conversion to 'set' would " + 
                                 "remove members with potentially unintended consequences")
        self._members = value
    
    def add(self, value):
        if value in self.members:
            raise ValueError(f"Cannot add '{value}' to self.members as it already exists")
        self.members.add(value)
    
    def replace(self, value):
        if not value in self.members:
            raise ValueError(f"Cannot replace '{value}' as it does not already exist in self.members")
        self.members.remove(value)
        self.members.add(value)
    
    def format(self, formatType):
        '''
        Waypoint function for choosing a specific formatter for writing output.
        '''
        ACCEPTED_FORMATS = ["binge", "cdhit", "corset", "mmseqs", "pantools", "2columnleft", "2columnright"]
        if not formatType in ACCEPTED_FORMATS:
            raise ValueError(f"formatType '{formatType}' is not among the accepted formats '{ACCEPTED_FORMATS}'")
        
        if formatType in ["mmseqs", "2columnleft"]:
            return self._format_2column("left")
        elif formatType in ["corset", "2columnright"]:
            return self._format_2column("right")
        elif formatType == "binge":
            return self._format_binge()
        elif formatType == "cdhit":
            return self._format_cdhit()
        elif formatType == "pantools":
            return self._format_pantools()
        else:
            raise NotImplementedError(f"Formatter for '{formatType}' is not yet implemented by .format()")
    
    def _format_2column(self, keyColumn):
        ACCEPTED_KEYCOLUMN = ["left", "right"]
        if not keyColumn in ACCEPTED_KEYCOLUMN:
            raise ValueError(f".parse_2column expects keyColumn to be in '{ACCEPTED_KEYCOLUMN}'; got '{keyColumn}' instead")
        
        out = []
        for member in self:
            if keyColumn == "left":
                out.append(f"{self.id}\t{member.id}")
            else:
                out.append(f"{member.id}\t{self.id}")
        return "\n".join(out)
    
    def _format_binge(self):
        out = []
        for member in self:
            if not hasattr(member, "binned"):
                raise ValueError(f"Member '{member}' in '{self.id}' cluster cannot be BINge formatted as it lacks .binned attribute")
            binnedStr = "binned" if member.binned else "unbinned"
            out.append(f"{self.id}\t{member.id}\t{binnedStr}")
        return "\n".join(out)
    
    def _format_cdhit(self):
        originalClustID = self.id.replace("_", " ")
        out = [f">{originalClustID}"]
        
        ongoingCount = -1
        for member in self:
            if not hasattr(member, "seqLength") or not hasattr(member, "seqDetails"):
                raise ValueError(f"Member '{member}' in '{self.id}' cluster cannot be CD-HIT formatted as it lacks " + 
                                 ".seqLength and/or .seqDetails attribute(s)")
            
            ongoingCount += 1
            out.append(f"{ongoingCount}\t{member.seqLength}, >{member.id}... {member.seqDetails}")
        return "\n".join(out)
    
    def _format_pantools(self):
        out = f"{self.id}:"
        for member in self:
            out += f" {member.id}"
        return out
    
    def __contains__(self, value):
        return value in self.members
    
    def __len__(self):
        return len(self.members)
    
    def __iter__(self):
        yield from sorted(self.members)
    
    def __len__(self):
        return len(self.members)
    
    def __repr__(self):
        return "<Cluster object;id='{0}';num_members={1}>".format(
            self.id, len(self.members)
        )

class ClusteringResult:
    def __init__(self, clustersFile, fileFormat, clusters=None):
        self.clustersFile = clustersFile
        self.fileFormat = fileFormat
        self.clusters = clusters
        
        # Set helper attributes
        self.isClusteringResult = True
        
        # Run init command
        self._parse()
    
    @property
    def clusters(self):
        return self._clusters
    
    @clusters.setter
    def clusters(self, value):
        if value == None:
            self._clusters = {}
            return
        if not isinstance(value, dict):
            raise TypeError(f"ClusteringResult.clusters must be a 'dict' not '{type(value).__name__}'")
        
        self._clusters = value
    
    @property
    def clustersFile(self):
        return self._clustersFile
    
    @clustersFile.setter
    def clustersFile(self, value):
        try:
            fullPath = os.path.abspath(value)
            isFile = os.path.isfile(fullPath)
        except TypeError:
            raise TypeError(f"Value of type '{type(value).__name__}' is not os.path compatible; " +
                            "clustersFile should be a 'str' or 'Path'")
        if not isFile:
            raise FileNotFoundError(f"clustersFile value '{value}' does not exist or is not a file")
        self._clustersFile = value
    
    @property
    def fileFormat(self):
        return self._fileFormat
    
    @fileFormat.setter
    def fileFormat(self, value):
        ACCEPTED_FORMATS = ["binge", "cdhit", "corset", "mmseqs", "pantools", "2columnleft", "2columnright"]
        if not value in ACCEPTED_FORMATS:
            raise ValueError(f"fileFormat '{value}' is not among the accepted formats '{ACCEPTED_FORMATS}'")
        self._fileFormat = value
    
    def _parse(self):
        '''
        Waypoint function for choosing a specific parser depending on self.fileFormat. If you are parsing
        a non-specific 2column file, you should call the .parse_2column directly to ensure that hasHeader
        is set appropriately, otherwise this will assume the standard of mmseqs or corset to NOT have a header.
        '''
        if self.fileFormat in ["mmseqs", "2columnleft"]:
            self.parse_2column("left")
        elif self.fileFormat in ["corset", "2columnright"]:
            self.parse_2column("right")
        elif self.fileFormat == "binge":
            self.parse_binge()
        elif self.fileFormat == "cdhit":
            self.parse_cdhit()
        elif self.fileFormat == "pantools":
            self.parse_pantools()
        else:
            raise NotImplementedError(f"Parser for '{self.fileFormat}' is not yet handled by ._parse()")
    
    def parse_2column(self, keyColumn, hasHeader=False):
        '''
        Parses a clustering format where the left or right column are cluster IDs/keys, and the opposite column
        are member IDs. This function is named generically to allow for application to any format which
        adopts this standard.
        
        Parameters:
            keyColumn -- a string indicating whether the key is the "left" or the "right" column
            hasHeader -- (OPTIONAL) a boolean indicating whether the file has a header; in practice this means
                         the first line will be ignored
        Sets:
            self.clusters -- a dictionary where keys are strings identifying the cluster and values are
                             Cluster objects (which themselves index Member objects)
        '''
        # Validate keyColumn
        ACCEPTED_KEYCOLUMN = ["left", "right"]
        if not keyColumn in ACCEPTED_KEYCOLUMN:
            raise ValueError(f".parse_2column expects keyColumn to be in '{ACCEPTED_KEYCOLUMN}'; got '{keyColumn}' instead")
        
        # Parse file
        self.clusters = {}
        skipHeader = hasHeader
        with read_gz_file(self.clustersFile) as fileIn:
            for line in fileIn:
                l = line.rstrip("\r\n\t;, ")
                
                # Handle header
                if skipHeader: # only True on the first iter if hasHeader == True
                    skipHeader = False
                    continue
                # Skip comments
                if l.startswith("#"):
                    continue
                # Handle content
                if l != "": # ignore empty lines
                    delimiter = "\t" if "\t" in l else "," if "," in l else None
                    if delimiter == None:
                        raise FormatError(f"'{self.clustersFile}' is not TSV or CSV formatted; offending line is '{l}'")
                    
                    sl = l.split(delimiter)
                    if len(sl) != 2:
                        raise FormatError(f"'{self.clustersFile}' does not have a 2-column format; offending line is '{l}'")
                    
                    if keyColumn == "left":
                        clustID, seqID = sl
                    else:
                        seqID, clustID = sl
                    
                    self.clusters.setdefault(clustID, Cluster(clustID))
                    self.clusters[clustID].add(Member(seqID))
    
    def parse_corset(self):
        '''
        Corset has a format where cluster IDs/keys are in the right column; this function acts as a
        shortcut to .parse_2column("right")
        '''
        return parse_2column("right", hasHeader=False)
    
    def parse_mmseqs(self):
        '''
        MMseqs2 has a format where cluster IDs/keys are in the left column; this function hence acts as a
        shortcut to .parse_2column("left")
        '''
        return parse_2column("left", hasHeader=False)
    
    def parse_binge(self):
        self.clusters = {}
        
        # Parse file
        self.clusters = {}
        skipHeader = True
        with read_gz_file(self.clustersFile) as fileIn:
            for line in fileIn:
                l = line.rstrip("\r\n\t;, ")
                
                # Skip comments
                if l.startswith("#"): # comment line precedes the header line
                    continue
                # Handle header
                if skipHeader: # only True on the first iter if hasHeader == True
                    skipHeader = False
                    continue
                # Handle content
                if l != "": # ignore empty lines
                    delimiter = "\t"
                    sl = l.split(delimiter)
                    if len(sl) != 3:
                        raise FormatError(f"'{self.clustersFile}' does not have a 3-column BINge format; offending line is '{l}'")
                    
                    clustID, seqID, clustType = sl
                    
                    self.clusters.setdefault(clustID, Cluster(clustID))
                    self.clusters[clustID].add(Member(seqID, binned=True if clustType == "binned" else False))
    
    def parse_cdhit(self):
        self.clusters = {}
        
        # Parse file
        self.clusters = {}
        with read_gz_file(self.clustersFile) as fileIn:
            for line in fileIn:
                l = line.rstrip("\r\n\t;, ")
                
                # Handle content
                if l != "": # ignore empty lines
                    delimiter = "\t"
                    sl = l.split(delimiter)
                    
                    # Handle cluster definition
                    if len(sl) == 1:
                        clustID = l[1:].replace(" ", "_") # drop the '>' prefix and remove whitespace
                        self.clusters.setdefault(clustID, Cluster(clustID))
                    # Handle cluster member
                    elif len(sl) == 2:
                        memberNum, seqDetails = sl
                        seqLength, seqDetails = seqDetails.split(", >") # we don't care about seqLength; variable named only for clarity of the process
                        seqID, seqDetails = seqDetails.split("... ")
                        self.clusters[clustID].add(Member(seqID, seqLength=seqLength, seqDetails=seqDetails, representative=True if seqDetails.endswith("*") else False))
                    else:
                        raise FormatError(f"'{self.clustersFile}' does not have 1 or 2 column CD-HIT TSVformat ; offending line is '{l}'")
    
    def parse_pantools(self):
        self.clusters = {}
        
        # Parse file
        self.clusters = {}
        with read_gz_file(self.clustersFile) as fileIn:
            for line in fileIn:
                l = line.rstrip("\r\n\t;, ")
                
                # Skip comments
                if l.startswith("#"): # first comment line contains parameter details
                    continue
                # Handle content
                if l != "": # ignore empty lines
                    delimiter = "\t"
                    clustID, seqsData = l.split(": ")
                    self.clusters.setdefault(clustID, Cluster(clustID))
                    
                    # Parse out memembers
                    for seqID in seqsData.split(" "):
                        originalID, fileNumber = seqID.split("#")
                        self.clusters[clustID].add(Member(seqID, originalSeqID=originalID, fileNumber=f"#{fileNumber}"))
    
    def keys(self):
        return self.clusters.keys()
    
    def values(self):
        return self.clusters.values()
    
    def items(self):
        return self.clusters.items()
    
    def __contains__(self, value):
        if hasattr(value, "isCluster") and value.isCluster:
            return value.id in self.clusters
        else:
            return value in self.clusters
    
    def __getitem__(self, key):
        return self.clusters[key]
    
    def __len__(self):
        return len(self.clusters)
    
    def __iter__(self):
        yield from self.clusters.values()
    
    def __repr__(self):
        return "<ClusteringResult object;clustersFile='{0}';fileFormat='{1}';num_clusters={2}>".format(
            self.clustersFile, self.fileFormat, len(self.clusters)
        )

def cluster_reformat(args):
    # Parse in filtering files
    shouldFilter = args.listFiles != [] or args.fastaFiles != []
    if shouldFilter:
        retainIfContains = set()
        for listFile in args.listFiles:
            with read_gz_file(listFile) as fileIn:
                for line in fileIn:
                    l = line.strip()
                    retainIfContains.add(l)
        for fastaFile in args.fastaFiles:
            with read_gz_file(fastaFile) as fileIn:
                for line in fileIn:
                    if line.startswith(">"):
                        seqID = line[1:].strip().split(" ")[0]
                        retainIfContains.add(seqID)
    
    # Load and reformat the clustering result
    clusterResult = ClusteringResult(args.clusterFile, fileFormat=args.inputFileFormat)
    with GzCapableWriter(args.outputFileName) as fileOut:
        # Write header if applicable
        if args.outputFileFormat == "binge":
            fileOut.write("#BINge clustering information file\ncluster_num\tsequence_id\tcluster_type\n")
        
        # Write content
        for cluster in clusterResult:
            if shouldFilter:
                if not any([ member.id in retainIfContains for member in cluster ]):
                    continue # if no overlap was found between this cluster and retainIfContains
            
            fileOut.write(cluster.format(args.outputFileFormat) + "\n")
