#! python3

import os, math, re, sys
import pandas as pd
from ncls import NCLS
from collections import Counter

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import read_gz_file, write_conditionally, parse_annotation_table
from fasta import FASTATarium, Sequence
from coordinates import Coordinates

class GFF3Feature:
    IMMUTABLE = ["ID", "ftype"] # these attributes should never change once set
    HEADER_FORMAT = ["contig", "source", "ftype", "start", "end", "score", "strand", "frame", "attributes"]
    CONTROLLED_ATTRIBUTES = ["ID", "Parent"] # these attributes are controlled by the GFF3Feature class
    
    def __init__(self, ID, ftype, start=None, end=None, strand=None, contig=None,
                 source=None, score=None, frame=None, attributes=None, children=None, parents=None,
                 isInferred=False):
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
        self.parents = parents
        self.isInferred = isInferred # flag for checking if this feature was parsed or inferred
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
            return self._source
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
    
    @score.setter
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
        attrDict = {}
        
        # Set controlled attributes
        attrDict["ID"] = self.ID
        if len(self.parents) != 0:
            attrDict["Parent"] = ",".join(self.parents)
        
        # Add other attributes
        if self._attributes != None:
            for key, value in self._attributes.items():
                if not key in GFF3Feature.CONTROLLED_ATTRIBUTES:
                    attrDict[key] = value
        
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
    
    @property
    def parents(self):
        return self._parents
    
    @parents.setter
    def parents(self, value):
        self._parents = value if isinstance(value, set) \
                        else set(value) if isinstance(value, list) \
                        else set([value]) if isinstance(value, str) \
                        else set()
    
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
    
    def find_all_children(self, foundChildren=None):
        '''
        Returns _all_ children contained under this feature.
        '''
        if foundChildren is None:
            foundChildren = []
        
        if len(self.children) != 0:
            for child in self.children:
                foundChildren.append(child)
                child.find_all_children(foundChildren)
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
    
    def format(self, alreadyFound, recursion=None):
        '''
        This method will attempt to render a GFF3-correct format of the
        data this GFF3Feature contains. Several assumptions are made which, if you
        haven't done anything truly weird, will hold true.
        
        alreadyFound should always be empty for the parent-level feature
        
        Parameters:
            alreadyFound -- a set maintained by the calling function which tracks
                            feature .ID values that have already been found/formatted.
                            Necessary to accommodate multiparent features and prevent
                            duplicate outputs.
        '''
        # Recursion management
        if recursion == None:
            recursion = []
        if self.ID in alreadyFound:
            return
        alreadyFound.add(self.ID)
        
        # Depth-first formatting of details
        recursion.append(str(self))
        for child in self.children:
            child.format(alreadyFound, recursion)
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
        attrsToShow = ["ID", "ftype", "contig", "start", "end", "strand", "parents"]
        
        for attr in attrsToShow:
            # if attr == "coords":
            #     reprPairs.append("coords=[{0}, {1}]".format(self.start, self.end))
            # elif attr == "strand":
            #     reprPairs.append("strand={0}".format(self.strand if self.strand != None else "."))
            # elif attr == "parents":
            #     reprPairs.append("parents={0}".format(self.parents))
            # else:
            #     reprPairs.append("{0}={1}".format(attr, self.__dict__[attr]))
            reprPairs.append("{0}={1}".format(attr, getattr(self, attr)))
        
        return "<{0};{1}>".format(";".join(reprPairs),
                                  f"children=[{', '.join([child.ID for child in self.children])}]")

class UnsortedGff3Error(Exception):
    pass
class DuplicateFeatureError(Exception):
    pass

class GFF3Tarium:
    PARENT_INFERENCE = {
        "CDS": "mRNA",
        "exon": "mRNA",
        "mRNA": "gene",
        "lnc_RNA": "gene",
        "Product": "gene" # Product is a special case, but we treat it as a gene parent
        # "gene": None  # Gene is the top-level feature, no parent should be inferred
    }
    SEMICOLON_REGEX = re.compile(r";{2,}")
    
    @staticmethod
    def format_attributes(attributes):
        '''
        Despite how overengineered this may seem, it is necessary to handle cases where multiple redundant
        semi-colons exist, or when e.g., chemical names are embedded which may contain semi colons, not as a 
        delimiter, but as part of the value.
        '''
        # Run an initial cleaning of the string
        attributes = attributes.strip("\r\n\t; ")
        multiSemiColons = sorted(set(GFF3Tarium.SEMICOLON_REGEX.findall(attributes)), key=len, reverse=True)
        for multipleSemiColonString in multiSemiColons:
            attributes = attributes.replace(multipleSemiColonString, ";")
        
        # Parse string into pairs of key:value attributes
        splitAttributes = []
        for a in attributes.split("="):
            if ";" in a:
                splitAttributes += a.rsplit(";", maxsplit=1)
            else:
                splitAttributes.append(a)
        
        return { splitAttributes[i]: splitAttributes[i+1] for i in range(0, len(splitAttributes)-(len(splitAttributes)%2), 2) }
    
    def __init__(self, fileLocation, deduplicate=False):
        self.fileLocation = fileLocation
        self.ftypes = {} # stores ftype:[featureIDs...]
        self.parentFtypes = set() # stores keys in self.ftypes that point to parent-level features
        self.features = {} # stores featureID:GFF3Feature
        self.contigs = set() # stores unique feature .contig values
        
        self.ncls = None
        self._nclsType = None
        self._nclsIndex = None
        self._nclsMax = None
        
        self.parse_gff3(self.fileLocation, deduplicate)
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
    
    def _format_new_feature_id(self, ftype):
        '''
        Parameters:
            ftype -- a string of a feature type that exists in self.ftypes
        Returns:
            newFeatureID -- a (most likely to be) unique feature ID, although
                            this is not guaranteed for weird GFF3s that are
                            trying to break things as edge cases
        '''
        return f"{ftype}.{len(self.ftypes[ftype]) + 1}"
    
    def _get_unique_feature_id(self, inputID, separator="."):
        ongoingCount = 1
        featureID = inputID
        while featureID in self.features:
            featureID = f"{inputID}{separator}{ongoingCount}"
            if not featureID in self.features:
                break
            ongoingCount += 1
        return featureID
    
    def parse_gff3(self, gff3File, deduplicate=False):
        '''
        Parses a GFF3 file and populates the graph with features.
        
        Parameters:
            gff3File -- a string pointing to a GFF3 file location
                        to parse and populate this object with.
            deduplicate -- (OPTIONAL) a boolean indicating whether
                           duplicate feature IDs should raise an error
                           (False) or if we should automatically
                           deduplicate feature IDs (True) to avoid errors
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
                    featureID = self._format_new_feature_id(ftype)
                else:
                    featureID = attributes["ID"]
                
                # Get the parent ID(s) attribute
                if not "Parent" in attributes:
                    parentIDs = []
                    self.parentFtypes.add(ftype)
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
                    existingFeature = self.features[featureID]
                    
                    # Handle inferred features
                    if existingFeature.isInferred:
                        # Check that inferred details are correct
                        if existingFeature.ftype != ftype:
                            raise UnsortedGff3Error(f"Feature ID '{featureID}' has a different type '{ftype}' than previously inferred '{existingFeature.ftype}'")
                        if existingFeature.contig != contig:
                            raise UnsortedGff3Error(f"Feature ID '{featureID}' has a different contig '{contig}' than previously inferred '{existingFeature.contig}'")
                        
                        # Update feature details
                        existingFeature.start = start
                        existingFeature.end = end
                        existingFeature.strand = strand
                        existingFeature.source = source
                        existingFeature.score = score
                        existingFeature.frame = frame
                        existingFeature.attributes = attributes
                        existingFeature.parents.update(parentIDs) # add parents to existing set
                        existingFeature.isInferred = False # turn off flag since we found the feature
                        self.add(existingFeature) # re-add to ensure parents are updated correctly
                    
                    # Handle truly duplicated features
                    else:
                        if deduplicate:
                            newFeatureID = self._format_new_feature_id(ftype)
                            feature.ID = newFeatureID # although ID should be "immutable" it is not added into this graph yet so this is acceptable
                            self.add(feature)
                        else:
                            raise DuplicateFeatureError(f"Feature ID '{featureID}' occurs more than once in file '{self.fileLocation}'")
    
    def add(self, feature):
        # Store a new feature within the graph
        if not feature.ID in self.features:
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
                                                contig=feature.contig, children=[feature],
                                                isInferred=True)
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
        nclsMax = 0
        ongoingCount = 0
        for indexType in typeToIndex:
            for featureID in self.ftypes[indexType]:
                feature = self.features[featureID]
                starts.append(feature.start)
                ends.append(feature.end + 1) # NCLS indexes 0-based like a range so +1 to make this more logically compliant with gff3 1-based system
                ids.append(ongoingCount)
                nclsIndex[ongoingCount] = feature
                ongoingCount += 1
                nclsMax = max(nclsMax, feature.end + 1)
        
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        ncls = NCLS(starts.values, ends.values, ids.values)
        
        # Associate it to this instance
        self.ncls = ncls
        self._nclsType = typeToIndex
        self._nclsIndex = nclsIndex
        self._nclsMax = nclsMax # used to know the upper bound for querying all ranges indexed
    
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
                  "downstream annotarium functionality.")
            for ftype, count in danglingFeatures.items():
                print(f"# {count} '{ftype}' feature{'s have' if count > 1 else ' has'} no parents or children")
    
    @property
    def parents(self):
        '''
        Iterates through this GFF3Tarium object to yield all features that are
        at the parent level. Beginning with self.parentFtypes, a set which
        only contains ftypes with at least one parent-level feature, we loop
        through all features under those ftype(s) and yield features which
        lack any parents.
        '''
        for ftype in sorted(self.parentFtypes):
            for featureID in self.ftypes[ftype]:
                feature = self[featureID]
                if len(feature.parents) == 0: # we don't trust the GFF3 file to have an ftype ALWAYS be a parent or a child
                    yield feature
    
    def get_feature_parents(self, feature, parents=None):
        '''
        Climbs up any parent(s) to get to the top-level of a given feature. Will return
        the input feature if it is already at the top-level.
        '''
        if parents is None:
            parents = []
        
        if len(feature.parents) == 0:
            parents.append(feature)
        else:
            for parentID in feature.parents:
                return self.get_feature_parents(self[parentID])
        
        return parents
    
    def write(self, outputFileName, typesToWrite=None, idsToWrite=None):
        '''
        Writes this GFF3 object to file. For typesToWrite, you should specify the
        highest level feature type for each feature you are interested in, if you want
        to avoid output duplication. For example, to output all mRNA features, you probably
        want to specify the 'gene' parent-level feature.
        
        Parameters:
            outputFileName -- a string indicating the location to write to; overwriting is
                              not permitted
            typesToWrite -- None to indicate that all feature types should be output, OR
                            a list of strings indicating self.type values to output.
            idsToWrite -- None to indicate that all features should be output, OR
                          a list of strings indicating feature IDs to output.
        '''
        # Validate arguments
        if os.path.exists(outputFileName):
            raise FileExistsError(f"'{outputFileName}' already exists; GFF3Tarium will not .write() to here")
        if typesToWrite != None and idsToWrite != None:
            raise ValueError("GFF3Tarium.write() can't handle typesToWrite and idsToWrite both being set")
        
        # Get the contig order to write
        contigOrder = sorted(self.contigs)
        
        # Get feature IDs in orderable fashion
        found = set() # prevents us from writing the same parent twice with multi-parent features
        sortOrder = []
        if typesToWrite is None and idsToWrite is None: # get all parent-level features
            for feature in self.parents:
                sortOrder.append((contigOrder.index(feature.contig), feature.start, feature.ID))
        elif idsToWrite != None: # just use the given IDs
            for featureID in idsToWrite:
                parentFeatures = self.get_feature_parents(self[featureID]) # climb up to the parent-level if provided a child ID
                for parentFeature in parentFeatures:
                    if not parentFeature.ID in found:
                        sortOrder.append((contigOrder.index(parentFeature.contig), parentFeature.start, parentFeature.ID))
                        found.add(parentFeature.ID)
        else: # limit ourselves to the feature types provided
            for ftype in typesToWrite:
                if not ftype in self.ftypes:
                    raise ValueError(f"typesToWrite received '{ftype}' which is not part of this GFF3Tarium object")
                
                for feature in [ self[fid] for fid in self.ftypes[ftype] ]:
                    for parentFeature in self.get_feature_parents(feature): # climb up to the parent-level if provided child ftypes
                        if not parentFeature.ID in found:
                            sortOrder.append((contigOrder.index(parentFeature.contig), parentFeature.start, parentFeature.ID))
                            found.add(parentFeature.ID)
        sortOrder.sort(key = lambda x: (x[0], x[1]))
        
        # Write ordered features to file
        with open(outputFileName, "w") as fileOut:
            haveWritten = set()
            for _, _, featureID in sortOrder:
                output = self[featureID].format(haveWritten)
                if output: # there may theoretically be an edge case where .format() returns None; this handles it
                    fileOut.write(output)
    
    def reset_id(self, feature, newIDPrefix, merging=False):
        '''
        Takes a feature, which may be part of this object or not, and sets its .ID attributes to avoid
        conflict with any features in this object.
        
        Parameters:
            feature -- a GFF3Feature object which may or may not be part of this GFF3Tarium object.
            newIDPrefix -- a string to set as the prefix for setting new .ID values for the input features
                           with cascading down to children.
            merging -- (OPTIONAL) a boolean indicating whether this feature is being merged from another
                       GFF3Tarium object into this one. This flag will change how we set parent values.
        '''
        newID = self._get_unique_feature_id(newIDPrefix, separator="_")
        feature.ID = newID
        for i, child in enumerate(feature.children):
            if merging:
                child.parents = set([newID]) # set new parents value
            else:
                child.parents = set([newID] + [ p for p in child.parents if p != feature.ID ]) # wipe previous parent ID, add new one
            self.reset_id(child, f"{newID}.{child.ftype}{i+1}")
    
    def merge_feature(self, feature):
        '''
        Merges a feature from another GFF3Tarium object into this one. Functionally, this lets us apply
        slightly different logic to what it seen in .add() to accommodate this scenario better. Specifically,
        we expect the merged feature to be parent-level, so we want to add its children into the graph
        rather than inferring previously-unseen parents as might occur in an unsorted GFF3.
        '''
        if feature.ID in self.features:
            raise KeyError(f"'{feature.ID}' cannot merge into this GFF3Tarium object as its ID is a duplicate")
        
        # Store the new feature within the graph
        self.parentFtypes.add(feature.ftype)
        self.contigs.add(feature.contig)
        
        self.ftypes.setdefault(feature.ftype, [])
        self.ftypes[feature.ftype].append(feature.ID)
        self.features[feature.ID] = feature
        
        # Update graph features with parent-child relationships
        for childFeature in feature.children:
            childFeature.parents = feature.ID # this is necessary despite reset_id setting the parents value; unsure why
            self.merge_feature(childFeature)
    
    def merge_gff3(self, otherGff3, isoPct=0.3, dupePct=0.6):
        '''
        Merges another GFF3Tarium object into this one. Both GFF3s should have NCLS indexing applied
        in a manner which makes sense for finding overlaps.
        
        Percentage cutoffs dictate how merging occurs.
        - Features which overlap < the isoform percentage will be treated as new features. 
        - Features which overlap >= the isoform percentage cutoff but < the duplicate percentage will
        be handled as alternative isoforms.
        - Features which overlap >= the duplicate percentage will be excluded.
        
        Parameters:
            otherGff3 -- another GFF3Tarium instance
            isoPct -- (OPTIONAL) a float from 0->1 inclusive indicating the lower cutoff for a feature
                      to be considered a potential isoform variant
            dupePct -- (OPTIONAL) a float from 0->1 inclusive indicating the lower cutoff for a feature
                       to be deemed duplicated and omitted from merging
        Returns:
            isoforms -- a dictionary where keys are parent IDs from this GFF3 object, and values
                        are lists of feature IDs that were merged as isoforms into this GFF3
            additions -- a list of sequence IDs that were merged as new features into this GFF3
        '''
        if isoPct < 0 or isoPct > 1:
            raise ValueError(f"isoPct={isoPct} is not a valid value for merge_gff3()")
        if dupePct < 0 or dupePct > 1:
            raise ValueError(f"dupePct={dupePct} is not a valid value for merge_gff3()")
        
        # Establish data storage for knowing what changes have been made
        isoforms = {}
        additions = []
        
        for parent2Feature in otherGff3.parents:
            p2Coordinates = Coordinates([ (e.start, e.end) for f in parent2Feature.find_with_children("exon") for e in f.exon ])
            
            # Skip empty features
            if len(p2Coordinates) == 0:
                continue
            
            # Store a new feature if there are no overlaps
            parent1Features = self.ncls_finder(parent2Feature.start, parent2Feature.end, "contig", parent2Feature.contig)
            if len(parent1Features) == 0:
                additions.append(parent2Feature.ID)
                self.reset_id(parent2Feature, parent2Feature.ID, merging=True)
                self.merge_feature(parent2Feature)
                continue
            
            # Handle overlapping features
            overlapResults = {}
            for p1Feature in parent1Features:
                p1Coordinates = Coordinates([ (e.start, e.end) for f in p1Feature.find_with_children("exon") for e in f.exon ])
                
                # Calculate overlap percentage for these child features
                firstPct, secondPct = p1Coordinates.overlap_percent(p2Coordinates)                
                overlapPct = max([firstPct, secondPct])
                
                # Categorise what we might do according to duplication cutoffs
                if overlapPct > dupePct:
                    overlapResults[p1Feature.ID] = "duplicate"
                    continue
                elif overlapPct >= isoPct: # isoform sweetspot
                    overlapResults[p1Feature.ID] = "isoform"
                    pass
                else:
                    overlapResults[p1Feature.ID] = "new"
            
            # Skip undesired features
            if any([ v == "duplicate" for v in overlapResults.values() ]): # do not merge a feature with any overlaps
                continue
            if sum([ v == "isoform" for v in overlapResults.values() ]) > 1: # do not merge competing isoforms
                continue
            
            # Merge isoforms
            if any([ v == "isoform" for v in overlapResults.values() ]):
                p1FeatureID = [ k for k,v in overlapResults.items() if v == "isoform" ][0] # there will only be 1 match
                p2Features = parent2Feature.find_with_children("exon")
                for p2Feature in p2Features:
                    isoforms.setdefault(p1FeatureID, [])
                    isoforms[p1FeatureID].append(p2Feature.ID) # logging
                    
                    # Integrate this new feature into this GFF3Tarium object
                    p2Feature.parents = p1FeatureID
                    self.reset_id(p2Feature, p2Feature.ID, merging=True)
                    self.merge_feature(p2Feature)
                    
                    # Add this new feature as a child of the parent
                    p1Feature = self[p1FeatureID]
                    p1Feature.add_child(p2Feature)
            
            # Add new features
            else:
                additions.append(parent2Feature.ID) # logging
                self.reset_id(parent2Feature, parent2Feature.ID, merging=True)
                self.merge_feature(parent2Feature)
        return isoforms, additions
    
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

def gff3_merge(args):
    # Parse GFF3 with NCLS indexing
    gff3_1 = GFF3Tarium(args.gff3File)
    gff3_1.create_ncls_index(typeToIndex=list(gff3_1.parentFtypes))
    
    gff3_2 = GFF3Tarium(args.gff3File2)
    gff3_2.create_ncls_index(typeToIndex=list(gff3_2.parentFtypes))
    
    # Merge and write output GFF3
    isoforms, additions = gff3_1.merge_gff3(gff3_2, isoPct=args.isoformPercent, dupePct=args.duplicatePercent)
    gff3_1.write(args.outputFileName)
    
    # Print out basic statistics
    print((f"# '{args.gff3File2}' with {sum( 1 for p in gff3_2.parents )} parent-level " + 
           f"features was merged into '{args.gff3File}' which has {sum( 1 for p in gff3_1.parents )} " +
           "parent-level features"))
    print(f"# {len(additions)} new parent features were added, and {len(isoforms)} isoforms were added")
    
    # Optionally emit merge details if requested
    if args.outputDetailsName != None:
        with open(args.outputDetailsName, "w") as fileOut:
            fileOut.write(f"# {len(additions)} new features added")
            if len(additions) > 0:
                fileOut.write(", including:\n")
                for featureID in additions:
                    fileOut.write(f"{featureID}\n")
            else:
                fileOut.write("\n")
            
            fileOut.write(f"# {len(isoforms)} new isoforms added as children")
            if len(isoforms) > 0:
                fileOut.write(", including:\n")
                for geneID, newFeatureIDs in isoforms.items():
                    for newFeatureID in newFeatureIDs:
                        fileOut.write(f"{geneID} <- {newFeatureID}\n")
            else:
                fileOut.write("\n")

def gff3_filter(args):
    # Parse list file (if applicable)
    selectionValues = []
    if args.listFile != None:
        with open(args.listFile, "r") as fileIn:
            for line in fileIn:
                selectionValues.append(line.strip())
    
    # Mix in any values specified on command-line
    selectionValues.extend(args.values)
    
    # Remove duplicates and set None if no selection is to occur
    selectionValues = set(selectionValues)
    if len(selectionValues) == 0:
        selectionValues = None
    
    # Parse GFF3 with NCLS indexing
    gff3 = GFF3Tarium(args.gff3File)
    gff3.create_ncls_index(typeToIndex=list(gff3.parentFtypes))
    
    # Set upper bounds for any regions with end==None
    if args.regions != None:
        for region in args.regions:
            if region["end"] == None:
                region["end"] = gff3._nclsMax
    
    # Filter based on selection criteria
    passedIDs = []
    for parentFeature in gff3.parents:
        # Handle region selection
        if args.regions != None: # ignore region selection if == None
            isSelected = any([
                Coordinates.isOverlapping(parentFeature.start, parentFeature.end,
                                          region["start"], region["end"])
                for region in args.regions
                if parentFeature.contig == region["contig"]
            ])
            
            if args.retrieveOrRemove == "retrieve" and isSelected:
                passedIDs.append(parentFeature.ID)
                continue # passes selection criteria
            elif isSelected:
                continue # filter and remove
        
        # Handle value selection
        if selectionValues != None:
            isSelected = False
            
            for feature in [parentFeature] + parentFeature.find_all_children():
                for attribute in GFF3Feature.HEADER_FORMAT:
                    if attribute == "attributes":
                        for value in feature.attributes.values():
                            if value in selectionValues:
                                isSelected = True
                    else:
                        if getattr(feature, attribute) in selectionValues:
                            isSelected = True
                            break
            
            if args.retrieveOrRemove == "retrieve" and isSelected:
                passedIDs.append(parentFeature.ID)
                pass # passes selection criteria
            elif isSelected:
                continue # filter and remove
    
    # Exit if no features are selected
    if len(passedIDs) == 0:
        raise ValueError("No features would be written to file with your filtering criteria")
    
    # Write to output if we pass the previous selection checks
    gff3.write(args.outputFileName, idsToWrite=passedIDs)

def gff3_annotate(args):
    gff3 = GFF3Tarium(args.gff3File)
    
    # Parse and annotate all relevant attributes
    for dataDict in parse_annotation_table(args.tableFile):
        # Format the attributes to annotate within our GFF3
        attributeAnnotations = {}
        for column, attribute, delimiter in args.columnAttributeDelimiter:
            try:
                value = dataDict[column]
            except KeyError:
                raise KeyError(f"'{column}' as part of {column, attribute, delimiter} trio is not a table column")
            
            if delimiter != None:
                value = value.split(delimiter)[0].strip()
            attributeAnnotations[attribute] = value
        
        # Annotate these attributes into the GFF3
        try:
            feature = gff3[dataDict["leftcolumn"]]
        except KeyError:
            tableID = dataDict['leftcolumn']
            raise KeyError(f"'{tableID}' from annotation table was not found in your GFF3")
        for key, value in attributeAnnotations.items():
            feature._attributes[key] = value
    
    # Write to file
    gff3.write(args.outputFileName)

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

def gff3_to_gff3(args):
    gff3 = GFF3Tarium(args.gff3File, deduplicate=True)
    gff3.write(args.outputFileName)
