#! python3

import os, math
import pandas as pd
from ncls import NCLS
from collections import Counter

from .parsing import read_gz_file, write_conditionally
from .fasta import FASTATarium, Sequence

class GFF3Feature:
    IMMUTABLE = ["ID", "ftype"] # these attributes should never change once set
    HEADER_FORMAT = ["contig", "source", "ftype", "start", "end", "score", "strand", "frame", "attributes"]
    CONTROLLED_ATTRIBUTES = ["ID", "Parent"] # these attributes are controlled by the GFF3Feature class
    
    def __init__(self, ID, ftype, start=None, end=None, strand=None, contig=None,
                 source=None, score=None, frame=None, attributes=None, children=None, parents=None):
        self.ID = ID
        self.ftype = ftype
        self.start = start
        self.end = end
        self.strand = strand
        self.contig = contig
        self.source = source
        self.score = score
        self.frame = frame
        self.attributes = attributes
        
        self._children = []
        self.children = children
        self.parents = parents if isinstance(parents, set) \
                       else set(parents) if isinstance(parents, list) \
                       else set([parents]) if isinstance(parents, str) \
                       else set()
        self.isGFF3Feature = True # flag for easier type checking
    
    @staticmethod
    def make_ftype_case_appropriate(ftype):
        if ftype.lower() == "gene":
            return "gene"
        elif ftype.lower() == "mrna":
            return "mRNA"
        elif ftype.lower() == "exon":
            return "exon"
        elif ftype.lower() == "cds":
            return "CDS"
        elif ftype.lower() == "lnc_rna":
            return "lnc_RNA"
        elif ftype.lower() == "product":
            return "Product"
        else:
            return ftype
    
    @property
    def start(self):
        if self._start != None:
            return self._start
        else:
            if len(self.children) == 0:
                return None
            return min([ child.start for child in self.children ])
    
    @start.setter
    def start(self, value):
        if value != None:
            try:
                value = int(value)
            except ValueError:
                raise ValueError(f"Start value of '{value}' is not an integer and is invalid for GFF3 formatting")
            if value < 1:
                raise ValueError(f"Start value '{value}' cannot be zero or negative; GFF3 positions are 1-based")
            self._start = value
        else:
            self._start = None
    
    @property
    def end(self):
        if self._end != None:
            return self._end
        else:
            if len(self.children) == 0:
                return None
            return max([ child.end for child in self.children ])
    
    @end.setter
    def end(self, value):
        if value != None:
            try:
                value = int(value)
            except ValueError:
                raise ValueError(f"End value of '{value}' is not an integer and is invalid for GFF3 formatting")
            if value < 1:
                raise ValueError(f"End value '{value}' cannot be zero or negative; GFF3 positions are 1-based")
            self._end = value
        else:
            self._end = None
    
    @property
    def strand(self):
        if self._strand != None and self._strand != ".":
            return self._strand
        else:
            childStrands = [ child.strand for child in self.children if child.strand != None and child.strand != "." ]
            if len(childStrands) == 0:
                return "+" # default strand if no children have a strand
            else:
                mostCommonStrand = Counter(childStrands).most_common(1)[0][0] # basic majority vote
                return mostCommonStrand
    
    @strand.setter
    def strand(self, value):
        ACCEPTED_STRANDS = ["+", "-", "."] # "." might represent unknown strand
        if value != None:
            if not value in ACCEPTED_STRANDS:
                raise ValueError(f"Strand value '{value}' is not recognised; should be one of {ACCEPTED_STRANDS}")
            self._strand = value
        else:
            self._strand = None
    
    @property
    def source(self):
        if self._source != None:
            return self._strand
        else:
            return "."
    
    @source.setter
    def source(self, value):
        if value == None:
            pass
        elif not isinstance(value, str):
            raise ValueError(f"GFF3 source must be a string, not {type(value)}")
        self._source = value
    
    @property
    def score(self):
        if self._score != None:
            return self._score
        else:
            return "."
    
    @source.setter
    def score(self, value):
        self._score = value
    
    @property
    def frame(self):
        if self._frame != None:
            return self._frame
        else:
            return "."
    
    @frame.setter
    def frame(self, value):
        ACCEPTED_FRAMES = ["0", "1", "2", "."] # "." might represent unknown or irrelevant frame
        if value != None:
            if not value in ACCEPTED_FRAMES:
                raise ValueError(f"Frame value '{value}' is not recognised; should be one of {ACCEPTED_FRAMES}")
            self._frame = value
        else:
            self._frame = None
    
    @property
    def attributes(self):
        # Get the base of this attributes dict
        if self._attributes == None:
            attrDict = {}
        else:
            attrDict = self._attributes
        
        # Set controlled attributes
        attrDict["ID"] = self.ID
        if len(self.parents) != 0:
            attrDict["Parent"] = ",".join(self.parents)
        return attrDict
    
    @attributes.setter
    def attributes(self, value):
        if value == None:
            self._attributes = None
        else:
            if not isinstance(value, dict):
                raise TypeError(f"Attributes value must a dict, not '{type(value)}'")
            cleanedValue = { k:v for k,v in value.items() if not k in GFF3Feature.CONTROLLED_ATTRIBUTES }
            
            self._attributes = cleanedValue
    
    @property
    def children(self):
        return self._children
    
    @children.setter
    def children(self, value):
        if isinstance(value, list):
            for child in value:
                self.add_child(child) # ensure each child is added properly
        else:
            self.add_child(value) # type is not validated, but we assume it's a GFF3Feature object
    
    def add_child(self, childFeature):
        '''
        Adds a child feature to this feature's children list.
        
        Parameters:
            childFeature -- a GFF3Feature object to add as a child of this feature.
        '''
        childFeature.parents.add(self.ID) # ensure the child knows its parent
        self.children.append(childFeature)
        self.__dict__.setdefault(childFeature.ftype, [])
        self.__dict__[childFeature.ftype].append(childFeature)
    
    def find_with_children(self, attribute, foundChildren=None):
        '''
        Finds ("catches") all children contained under this feature which have the specified attribute.
        
        Parameters:
            attribute -- a string indicating an attribute that a child feature should have.
                         For example, to retrieve all "CDS" children at any level under this
                         feature, you would provide "CDS" as the attribute.
            foundChildren -- a list that stores values during recursion, or None for the top-level
                             recursion.
        '''
        if foundChildren is None:
            foundChildren = []
        
        if len(self.children) != 0:
            for child in self.children:
                if hasattr(child, attribute):
                    foundChildren.append(child)
                else:
                    return child.find_with_children(attribute, foundChildren)
        return foundChildren
    
    def length(self, attribute):
        '''
        Sums the length of all .attribute children contained under this feature.
        
        Parameters:
            attribute -- a string indicating the attribute under which children are indexed
        Returns
            length -- the summed length of children
        '''
        if not hasattr(self, attribute):
            return None
        thisLength = 0
        for child in getattr(self, attribute):
            thisLength += (child.end - child.start + 1)
        return thisLength
    
    def format(self, recursion=None):
        '''
        This method will attempt to render a GFF3-correct format of the
        data this GFF3Feature contains. Several assumptions are made which, if you
        haven't done anything truly weird, will hold true.
        '''
        if recursion == None:
            recursion = []
        
        # Depth-first formatting of details
        recursion.append(str(self))
        for child in self.children:
            child.format(recursion)
        return "\n".join(recursion) + "\n"
    
    def as_sequence(self, fastaObj):
        '''
        Making use of the class objects defined in this repository, obtains
        the corresponding sequence portion identified by this GFF3Feature.
        
        Parameters:
            fastaObj -- a FASTATarium or Sequence object
        Returns:
            seqObj -- a Sequence object with start and end defined by
                      self.start and self.end (and optionally self.contig)
        '''
        if hasattr(fastaObj, "isFASTATarium") and fastaObj.isFASTATarium:
            return fastaObj(self.contig, self.start-1, self.end) # counteract 1-based GFF3 numbering
        elif hasattr(fastaObj, "isSequence") and fastaObj.isSequence:
            return fastaObj[self.start-1:self.end] # counteract as well
        else:
            raise TypeError(f"Cannot obtain '{type(fastaObj)}' type as sequence")
    
    def as_gene_model(self, fastaObj, sequenceType="CDS"):
        '''
        Making use of the class objects defined in this repository, this function
        will order children of CDS or exon type and format this as a Sequence
        object corresponding to the annotation this feature represents.
        '''
        ACCEPTED_TYPES = ["CDS", "exon"]
        if not sequenceType in ACCEPTED_TYPES:
            raise ValueError(f"Feature cannot be ordered by '{sequenceType}'; must be in the list '{ACCEPTED_TYPES}'")
        if not hasattr(self, sequenceType):
            raise ValueError(f"'{self.ID}' lacks any children with '{sequenceType}' type")
        
        # Get the sequence as a string
        startingFrame = None
        sequence = ""
        for subFeature in sorted(getattr(self, sequenceType), key = lambda x: (x.start, x.end)):
            sequence += str(subFeature.as_sequence(fastaObj))
            if startingFrame == None: # capture the first frame for +ve strand features
                startingFrame = subFeature.frame
        if self.strand == "-":
            'Sort is from genomic left->right; first exon of a -ve strand feature is the rightmost'
            startingFrame = subFeature.frame
        
        # Convert to Sequence object
        seqObj = Sequence(self.ID, sequence, startingFrame)
        
        # Reverse complement if necessary
        if self.strand == "-":
            seqObj = seqObj.reverse_complement()
        
        return seqObj
    
    def __str__(self):
        return "\t".join([
            str(getattr(self, x)) if x != "attributes" else
            ";".join([ f"{k}={v}" for k,v in getattr(self, x).items() ])
            for x in GFF3Feature.HEADER_FORMAT
        ])
    
    def __repr__(self):
        reprPairs = []
        attrsToShow = ["ID", "ftype", "contig", "coords", "strand", "parents"]
        
        for attr in attrsToShow:
            if attr == "coords":
                reprPairs.append("coords=[{0}, {1}]".format(self.start, self.end))
            elif attr == "strand":
                reprPairs.append("strand={0}".format(self.strand if self.strand != None else "."))
            else:
                reprPairs.append("{0}={1}".format(attr, self.__dict__[attr]))
        
        return "<{0};{1}>".format(";".join(reprPairs),
                                  f"children=[{', '.join([child.ID for child in self.children])}]")

class GFF3Tarium:
    PARENT_INFERENCE = {
        "CDS": "mRNA",
        "exon": "mRNA",
        "mRNA": "gene",
        "lnc_RNA": "gene",
        "Product": "gene" # Product is a special case, but we treat it as a gene parent
        # "gene": None  # Gene is the top-level feature, no parent should be inferred
    }
    
    @staticmethod
    def format_attributes(attributes):
        '''
        Despite how overengineered this may seem, it is necessary to handle cases where multiple redundant
        semi-colons exist, or when e.g., chemical names are embedded which may contain semi colons, not as a 
        delimiter, but as part of the value.
        '''
        splitAttributes = []
        for a in attributes.strip("\r\n\t; ").split("="):
            if ";" in a:
                splitAttributes += a.rsplit(";", maxsplit=1)
            else:
                splitAttributes.append(a)
        
        return { splitAttributes[i]: splitAttributes[i+1] for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2) }
    
    def __init__(self, fileLocation):
        self.fileLocation = fileLocation
        self.ftypes = {}
        self.features = {}
        self.contigs = set()
        
        self.ncls = None
        self._nclsType = None
        self._nclsIndex = None
        
        self.parse_gff3(self.fileLocation)
        self.isGFF3Tarium = True # flag for easier type checking
    
    @property
    def fileLocation(self):
        return self._fileLocation
    
    @fileLocation.setter
    def fileLocation(self, value):
        if not isinstance(value, str):
            raise ValueError("File location must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"GFF3 file '{value}' is not a file")
        
        self._fileLocation = value
    
    def longest_feature(self, feature):
        '''
        Receives a feature (either by its .ID or by its actual self) and finds
        the longest feature indexed under (and including) itself according to summed
        CDS (first preference) or exon length.
        
        Note that by checking its own self, this function can be used on any feature
        to simply find the longest representative. If you input a 'gene' feature, this
        function will likely return the longest mRNA child feature. If you input an 
        mRNA child feature, this function will likely return itself.
        
        Parameters:
            feature -- a GFF3Feature object OR the .ID of an object indexed by
                       this GFF3Tarium object.
        Returns:
            longestFeature -- a GFF3Feature object indexed by this GFF3Tarium object
                              which has the longest summed CDS or exon length.
        '''
        if isinstance(feature, str):
            feature = self[feature]
        
        # Determine the features to use for length calculation
        if hasattr(feature, "CDS"):
            featureLengths = [[feature, feature.length("CDS")]]
        elif hasattr(feature, "exon"):
            featureLengths = [[feature, feature.length("exon")]]
        else:
            childFeatures = feature.find_with_children("CDS")
            if len(childFeatures) != 0:
                featureLengths = [ [f, f.length("CDS")] for f in childFeatures ]
            else:
                childFeatures = feature.find_with_children("exon")
                if len(childFeatures) != 0:
                    featureLengths = [ [f, f.length("exon")] for f in childFeatures ]
                else:
                    featureLengths = []
        
        # Find the longest feature
        longest = [None, 0]
        for feature, length in featureLengths:
            if length > longest[1]:
                longest = [feature, length]
        return longest[0]
    
    def _get_unique_feature_id(self, inputID):
        ongoingCount = 1
        featureID = inputID
        while featureID in self.features:
            featureID = f"{inputID}.{ongoingCount}"
            if not featureID in self.features:
                break
            ongoingCount += 1
        return featureID
    
    def parse_gff3(self, gff3File):
        '''
        Parses a GFF3 file and populates the graph with features.
        
        Parameters:
            gff3 -- a GFF3 object to parse and populate this graph with.
        '''
        # Reset the graph
        self.fileLocation = gff3File
        self.ftypes = {}
        self.features = {}
        self.contigs = set()
        
        # Parse the GFF3 file into a graph structure
        lineCount = 0
        with read_gz_file(self.fileLocation) as fileIn:
            for line in fileIn:
                lineCount += 1
                sl = line.strip("\r\n\t;'\" ").split("\t")
                
                # Skip filler and comment lines
                if line.startswith("#") or len(sl) != 9:
                    continue
                
                # Extract information from this line
                contig, source, ftype, start, end, \
                    score, strand, frame, attributes = sl
                start = int(start)
                end = int(end)
                ftype = GFF3Feature.make_ftype_case_appropriate(ftype)
                attributes = GFF3Tarium.format_attributes(attributes)
                
                # Establish or populate tracking containers
                self.ftypes.setdefault(ftype, [])
                self.contigs.add(contig)
                
                # Get the ID attribute
                if not "ID" in attributes:
                    featureID = f"{ftype}.{len(self.ftypes[ftype]) + 1}"
                else:
                    featureID = attributes["ID"]
                
                # Get the parent ID(s) attribute
                if not "Parent" in attributes:
                    parentIDs = []
                else:
                    parentIDs = [ x.strip() for x in attributes["Parent"].split(",") ]
                
                # Create a feature object
                feature = GFF3Feature(ID=featureID, ftype=ftype,
                                      start=start, end=end, strand=strand,
                                      source=source, score=score, frame=frame, attributes=attributes,
                                      contig=contig, children=[], parents=parentIDs)
                
                # Index the feature if it doesn't already exist
                if not featureID in self.features:
                    self.add(feature)
                # Specifically handle exons or CDS which are allowed to have duplicated IDs
                elif ftype in ["exon", "CDS"]:
                    feature.ID = self._get_unique_feature_id(featureID)
                    self.add(feature)
                # Handle other duplicated feature types
                else:
                    "We assume that the GFF3 is unsorted if we reach this point, so we are detailing an existing feature"
                    feature = self.features[featureID]
                    
                    # Check that inferred details are correct
                    if feature.ftype != ftype:
                        raise ValueError(f"Unsorted GFF3 issue: Feature ID '{featureID}' has a different type '{feature.ftype}' than previously inferred '{ftype}'")
                    if feature.contig != contig:
                        raise ValueError(f"Unsorted GFF3 issue: Feature ID '{featureID}' has a different contig '{feature.contig}' than previously inferred '{contig}'")
                    
                    # Update feature details
                    feature.start = start
                    feature.end = end
                    feature.strand = strand
                    feature.parents.update(parentIDs) # add parents to existing set
                    self.add(feature) # re-add to ensure parents are updated correctly
    
    def add(self, feature):
        # Store feature within the graph
        if not feature.ID in self.features: # only if it doesn't already exist
            self.ftypes.setdefault(feature.ftype, []) # necessary if first occurrence of a subfeature preceeds its parent type
            self.ftypes[feature.ftype].append(feature.ID)
            self.features[feature.ID] = feature
        
        # Update graph features with parent-child relationships
        for parentID in feature.parents:
            # Associate the feature with its existing parents
            if parentID in self.features:
                self.features[parentID].add_child(feature)
            # Create a placeholder for the parent if it doesn't exist
            else:
                if feature.ftype in GFF3Tarium.PARENT_INFERENCE:
                    parentFeature = GFF3Feature(parentID, GFF3Tarium.PARENT_INFERENCE[feature.ftype],
                                                contig=feature.contig,
                                                children=[feature])
                    self.add(parentFeature)
                else:
                    raise ValueError("Your GFF3 is not sorted in top-down hierarchical order which has caused an error; " +
                                     f"I encountered a {feature.ftype} with ID '{feature.ID}' that has a parent '{parentID}' which has " + 
                                     f"not yet appeared in your GFF3 file. I am unsure what parent type to infer for " +
                                     f"'{feature.ftype}' features, so I cannot continue parsing. Sort your GFF3 file in " +
                                     "conventional top-down hierarchical order before trying again.")
    
    def create_ncls_index(self, typeToIndex=["gene"]):
        '''
        Creates an indexed NCLS structure that can be used to find range overlaps
        for the feature types of interest.
        
        Associates the created index to the .ncls field of this object instance.
        A hidden ._nclsIndex dictionary links the ncls indices to feature objects.
        
        Parameters:
            typeToIndex -- a string (case-sensitive) indicating the entry type
                           to index OR an iterable of strings indicating multiple
                           types to index.
        '''
        if isinstance(typeToIndex, str):
            typeToIndex = [typeToIndex]
        
        for indexType in typeToIndex:
            assert indexType in self.ftypes, \
                "'{0}' not found as a feature type within the parsed GFF3 ('{1}')".format(indexType, self.fileLocation)
        
        nclsIndex = {}
        starts, ends, ids = [], [], []
        
        # Add features of the specified type to our NCLS structure
        ongoingCount = 0
        for indexType in typeToIndex:
            for featureID in self.ftypes[indexType]:
                feature = self.features[featureID]
                starts.append(feature.start)
                ends.append(feature.end + 1) # NCLS indexes 0-based like a range so +1 to make this more logically compliant with gff3 1-based system
                ids.append(ongoingCount)
                nclsIndex[ongoingCount] = feature
                ongoingCount += 1
        
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        ncls = NCLS(starts.values, ends.values, ids.values)
        
        # Associate it to this instance
        self.ncls = ncls
        self._nclsType = typeToIndex
        self._nclsIndex = nclsIndex
    
    def ncls_finder(self, start, stop, field, value):
        '''
        Queries the NCLS structure to find Features that exist within the given
        start->stop range. Specifying the field and value will narrow results
        to only those that have a Feature .field with an equal (==) value.
        
        Parameters:
            start -- an integer indicating the start position of the feature to check
                     for overlaps
            end -- an integer indicating the end positon of the feature to check for
                   overlaps; this should be 1-based in GFF3 style e.g., a first
                   position of a feature would be start=1, end=1.
            field -- a string (case-sensitive) indicating the field of the Feature
                     object that we want to check. For example, if you want to find
                     features that overlap the given start->stop range on the contig
                     "X", you'd provide "contig" as the field so this function knows
                     to check the Feature.contig field for the value of "X".
            value -- a string (case-sensitive) indicating the value of the Feature
                     field we want to find. As in the above example, if you want to
                     find the value "X" within the .contig field, you'd provide "X" as
                     the value here.
        Returns:
            features -- a list containing Features that overlap the specified range.
                        These Features are NOT deepcopied, so handle them carefully.
        '''
        assert self.ncls != None and self._nclsIndex != None, \
            "Run create_ncls_index before you call this method!"
        
        overlaps = self.ncls.find_overlap(start, stop+1) # Although our ncls is already 1-based, find_overlap acts as a range. We need to +1 to keep everything logically 1-based.
        
        features = []
        for result in overlaps: # result == [start, end, index]
            feature = self._nclsIndex[result[2]]
            if feature.__dict__[field] == value:
                features.append(feature)
        
        # Return list
        return features
    
    def qc(self, typesToCheck=None):
        '''
        Runs a quality control check on the GFF3Tarium object to ensure that all
        features are properly linked.
        
        Prints a warning if any features are found that have no parents or children.
        
        Parameters:
            typesToCheck -- an iterable of strings indicating the feature types to check
                            for dangling features. If None, all feature types are checked.
        '''
        danglingFeatures = {}
        for feature in self:
            if typesToCheck == None or feature.ftype in typesToCheck:
                if len(feature.parents) == 0 and len(feature.children) == 0:
                    danglingFeatures.setdefault(feature.ftype, 0)
                    danglingFeatures[feature.ftype] = 1
        
        if len(danglingFeatures) != 0:
            print(f"WARNING: Parsing '{self.fileLocation}' resulted in dangling features with no parents or children, " +
                  "which is likely due to an unsorted or incorrectly formatted GFF3 file. This may cause issues with " +
                  "psQTL's functionality.")
            for ftype, count in danglingFeatures.items():
                print(f"# {count} '{ftype}' feature{'s have' if count > 1 else ' has'} no parents or children")
    
    def __getitem__(self, key):
        return self.features[key]
    
    def __len__(self):
        return len(self.features)
    
    def __iter__(self):
        return iter(self.features.values())
    
    def __contains__(self, item):
        return item.ID in self.features
    
    def has_key(self, key):
        return key in self.features
    
    def keys(self):
        return self.features.keys()
    
    def values(self):
        return self.features.values()
    
    def items(self):
        return self.features.items()
    
    def __repr__(self):
        return "<GFF3Tarium object;file='{0}';num_contigs={1};{2}>".format(
            self.fileLocation,
            len(self.contigs),
            ";".join(["num_{0}={1}".format(key, len(self.ftypes[key])) for key in self.ftypes.keys()])
        )

def gff3_stats(args):
    raise NotImplementedError("gff3 stats mode not yet implemented")
    return None

def gff3_to_fasta(args):
    fasta = FASTATarium(args.fastaFile)
    gff3 = GFF3Tarium(args.gff3File)
    
    with write_conditionally(args.outputFileNames["exon"]) as exonOut, write_conditionally(args.outputFileNames["cds"]) as cdsOut, write_conditionally(args.outputFileNames["protein"]) as protOut:
        for featureType in args.features:
            if not featureType in gff3.ftypes:
                raise ValueError(f"'{featureType}' not found within '{args.gff3File}'")
            
            for parentFeatureID in gff3.ftypes[featureType]:
                # Pick out the longest representative for this feature
                feature = gff3.longest_feature(parentFeatureID)
                
                # Format the sequence(s)
                if "exon" in args.types:
                    exonSequence = feature.as_gene_model(fasta, "exon")
                    exonSequence.description = parentFeatureID # propagate the original ID, not the representative
                else:
                    exonSequence = None
                
                if "CDS" in args.types or "protein" in args.types:
                    cdsSequence = feature.as_gene_model(fasta, "CDS")
                    cdsSequence.description = parentFeatureID
                    
                    if "protein" in args.types:
                        proteinSequence = cdsSequence.translate(args.translationTable)
                        proteinSequence.description = parentFeatureID
                    else:
                        proteinSequence = None
                else:
                    cdsSequence = None
                    proteinSequence = None
                
                # Write to file
                if args.outputFileNames["exon"]:
                    exonOut.write(exonSequence.format())
                if args.outputFileNames["cds"]:
                    cdsOut.write(cdsSequence.format())
                if args.outputFileNames["protein"]:
                    protOut.write(proteinSequence.format())

def gff3_to_tsv(args):
    def get_value(feature, key):
        try:
            return getattr(feature, key)
        except:
            try:
                return feature.attributes[key]
            except:
                return None
    
    gff3 = GFF3Tarium(args.gff3File)
    
    # Perform validations that require gff3 to be parsed
    if not args.forEach in gff3.ftypes:
        raise ValueError(f"-forEach value '{value}' is not a feature type in your GFF3")
    
    # Run the query operation
    mapping = {}
    for featureID in gff3.ftypes[args.forEach]:
        feature = gff3[featureID]
        
        # Get value for -map key
        mapValue = get_value(feature, args.map)
        if mapValue is None:
            raise KeyError(f"-map value '{args.map}' not found for a '{args.forEach}' feature with start={feature.start} end={feature.end}")
        mapping.setdefault(mapValue, { k:{} for k in args.to }) # use dict as a sorted set
        
        # Get values for all -to keys
        for toKey in args.to:
            toValue = get_value(feature, toKey)
            mapping[mapValue][toKey].setdefault(None if toValue is None else str(toValue), None)
    
    # Collapse key:value pairings
    for mapKey in mapping.keys():
        collapsedValue = {}
        for toKey, valuesDict in mapping[mapKey].items():
            values = list(valuesDict.keys())
            if values == [None]:
                collapsedValue[toKey] = args.nullChar # replace None with the nullChar
            else:
                collapsedValue[toKey] = args.sepChar.join([ x for x in values if not x is None ]) # purge any None values
        mapping[mapKey] = collapsedValue
    
    # Convert to pandas dataframe to quickly tabulate
    idMapDF = pd.DataFrame.from_dict(mapping, orient="index")
    idMapDF.to_csv(args.outputFileName, sep="\t", index_label=args.map, header=not args.noHeader)
    