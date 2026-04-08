#! python3

import os, sys
import pandas as pd
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from gff3tarium import GFF3Feature, GFF3Tarium
from parsing import read_gz_file, write_conditionally, parse_annotation_table, GzCapableWriter
from fasta import FASTATarium
from coordinates import Coordinates

def gff3_stats(args):
    def print_and_write(text, fileHandle):
        print(text)
        fileHandle.write(text + "\n")
    
    gff3 = GFF3Tarium(args.gff3File)
    
    # Emit details about this GFF3 while writing to file
    with GzCapableWriter(args.outputFileName) as fileOut:
        print_and_write(f"# Num contigs with annotations = {len(gff3.contigs)}", fileOut)
        
        print_and_write("# Parent features", fileOut)
        for ftype in sorted(gff3.parentFtypes):
            print_and_write("## Num '{0}' = {1}".format(ftype, len(gff3.ftypes[ftype])), fileOut)
        
        otherFtypes = set(gff3.ftypes.keys()).difference(gff3.parentFtypes)
        print_and_write("# Child features", fileOut)
        if len(otherFtypes) == 0:
            print_and_write("## No child features found", fileOut)
        else:
            for ftype in sorted(otherFtypes):
                print_and_write("## Num '{0}' = {1}".format(ftype, len(gff3.ftypes[ftype])), fileOut)

def gff3_merge(args):
    # Parse GFF3 with NCLS indexing
    gff3_1 = GFF3Tarium(args.gff3File)
    gff3_1.create_ncls_index(typeToIndex=list(gff3_1.parentFtypes))
    
    gff3_2 = GFF3Tarium(args.gff3File2)
    gff3_2.create_ncls_index(typeToIndex=list(gff3_2.parentFtypes))
    
    # Store details on each file before modification
    file1ParentNum = sum( 1 for p in gff3_1.parents )
    file2ParentNum = sum( 1 for p in gff3_2.parents )
    
    # Merge and write output GFF3
    isoforms, additions, rejections = gff3_1.merge_gff3(gff3_2, isoPct=args.isoformPercent, dupePct=args.duplicatePercent)
    gff3_1.write(args.outputFileName)
    
    # Print out basic statistics
    print(f"# File 1: '{args.gff3File}' has {file1ParentNum} parent-level features")
    print(f"# File 2: '{args.gff3File2}' has {file2ParentNum} parent-level features")
    
    print(f"# {len(additions)} new parent features were added into file 1")
    print(f"# {len(isoforms)} isoforms were added into file 1 parent-level features")
    print(f"# {len(rejections)} features from file 2 were rejected due to overlap with existing file 1 features")
    
    print(f"# The output file (file 2 merged into file 1) has {sum( 1 for p in gff3_1.parents )} parent-level features")
    
    # Optionally emit merge details if requested
    if args.outputDetailsName != None:
        with GzCapableWriter(args.outputDetailsName) as fileOut:
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
            
            fileOut.write(f"# {len(rejections)} features were rejected")
            if len(rejections) > 0:
                fileOut.write(", including:\n")
                for geneID, overlappingIDs in rejections.items():
                    formattedIDs = ", ".join(overlappingIDs)
                    fileOut.write(f"{geneID} \\ {formattedIDs}\n")
            else:
                fileOut.write("\n")

def gff3_pcr(args):
    # Load in requisite data
    fasta = FASTATarium(args.fastaFile)
    gff3 = GFF3Tarium(args.gff3File)
    variants = gff3.get_variants()
    
    # Obtain the feature
    try:
        originalFeature = gff3[args.modelIdentifier]
    except KeyError:
        raise KeyError(f"-m value '{args.modelIdentifier}' was not found as a feature ID in '{args.gff3File}'")
    
    # Obtain representative isoform if modelIdentifier is parent-level
    feature = gff3.longest_feature(originalFeature)
    if feature.ID != args.modelIdentifier:
        if len(originalFeature.children) > 1:
            print(f"# Note: '{args.modelIdentifier}' had its longest isoform '{feature.ID}' chosen implicitly; " + 
                  "if you want to choose a different isoform, specify that with -m instead")
    
    # Get the PCR model sequence objects
    sequenceObjects = feature.as_pcr_model(fasta, buffer=args.buffer, variantsObj=variants)
    
    # Write output
    with GzCapableWriter(args.outputFileName) as fileOut:
        for sequenceObject in sequenceObjects:
            fileOut.write(sequenceObject.format())

def gff3_rc(args):
    gff3 = GFF3Tarium(args.gff3File)
    
    # Parse genome for contig lengths
    with read_gz_file(args.fastaFile) as fileIn:
        genomeRecords = SeqIO.parse(fileIn, "fasta")
        lengthsDict = { record.id:len(record) for record in genomeRecords }
    
    # Reverse complement all relevant features
    for key, feature in gff3.features.items():
        if args.toRC == True or feature.contig in args.toRC:
            feature.reverse_this(lengthsDict[feature.contig])
    
    # Write modified output
    gff3.write(args.outputFileName)

def gff3_relabel(args):
    # Parse list file (if applicable)
    renameList = []
    if args.listFile != None:
        with read_gz_file(args.listFile) as fileIn:
            for line in fileIn:
                sl = line.strip("\r\n\t").split("\t")
                if not len(sl) == 2:
                    raise ValueError(f"List file is expected to have 2 columns; malformed line is '{line.rstrip()}'")
                renameList.append(sl)    
    
    # Parse GFF3 with NCLS indexing
    gff3 = GFF3Tarium(args.gff3File)
    
    # Iterate over features, make edits where necessary, and write modified output
    with GzCapableWriter(args.outputFileName) as fileOut:
        haveWritten = set()
        for featureID in gff3._sorted_features(): # give None defaults to output all feature IDs
            feature = gff3[featureID]
            
            # Apply contig suffix
            if args.contigSuffix != None:
                feature.contig += args.contigSuffix
                for childFeature in feature.find_all_children():
                    childFeature.contig += args.contigSuffix
            
            # Apply gene suffix
            if args.geneSuffix != None:
                newGeneID = feature.ID + args.geneSuffix
                gff3.update_feature(feature, newGeneID)
            
            # Produce output line
            line = feature.format(haveWritten)
            if line: # .format() may return None with multiparent features; this prevents duplicate child feature writing
                # Apply naive line substitutions
                for original, new in renameList:
                    line = line.replace(original, new)
                
                fileOut.write(line)

def gff3_filter(args):
    # Parse list file (if applicable)
    selectionValues = []
    if args.listFile != None:
        with read_gz_file(args.listFile) as fileIn:
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
                continue # if selected, retrieve this feature
            elif args.retrieveOrRemove == "remove" and not isSelected:
                passedIDs.append(parentFeature.ID)
                continue # if unselected, do not remove this feature
        
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
            elif args.retrieveOrRemove == "remove" and not isSelected:
                passedIDs.append(parentFeature.ID)
                continue # if unselected, do not remove this feature
    
    # Exit if no features are selected
    if len(passedIDs) == 0:
        raise ValueError("No features would be written to file with your filtering criteria")
    
    # Write to output if we pass the previous selection checks
    gff3.write(args.outputFileName, idsToWrite=passedIDs)

def gff3_annotate(args):
    gff3 = GFF3Tarium(args.gff3File)
    
    # Parse and annotate all relevant attributes
    warnedAlready = False
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
            #raise KeyError(f"'{tableID}' from annotation table was not found in your GFF3")
            if not warnedAlready:
                print(f"WARNING: '{tableID}' from annotation table was not found in your GFF3; further warnings will be suppressed")
                warnedAlready = True
            continue
        for key, value in attributeAnnotations.items():
            feature._attributes[key] = value
    
    # Write to file
    gff3.write(args.outputFileName)

def _reformat_miniprot(gff3, identifiers, parentFeature):
    '''
    Helper function for updating a miniprot GFF3 annotation to have full
    GFF3 fields (i.e., a 'gene' parent and 'exon' children). Functionalised
    for reuse in different gff3_mp_* functions.
    
    Parameters:
        gff3 -- a GFF3Tarium object, parsed from the miniprot GFF3 file
        identifiers -- a dictionary maintained by the caller which tracks the
                       number of times we've encountered a particular 'Target'
                       feature (which was then used as the base for our newID)
        parentFeature -- the miniprot mRNA feature which we want to tack 'gene'
                         and 'exon' features onto.
    Returns:
        newGeneParent -- a GFF3Feature of the new 'gene' parent under which the
                         (now renamed) parentFeature resides (and its associated
                         exon children)
    '''
    # Get a new gene ID
    originalID, _start, _end = parentFeature.attributes["Target"].split(" ")
    identifiers.setdefault(originalID, 0)
    identifiers[originalID] += 1
    newID = f"{originalID}.mp.gene{identifiers[originalID]}"
    
    # Create new gene parent feature
    newGeneParent = GFF3Feature("tempminiprot", "gene", start=parentFeature.start, end=parentFeature.end,
                                strand=parentFeature.strand, contig=parentFeature.contig,
                                source=parentFeature.source, score=parentFeature.score,
                                frame=parentFeature.frame, attributes=parentFeature.attributes,
                                children=[parentFeature])
    parentFeature.parents = set([newID])
    gff3.update_feature(newGeneParent, newID, merging=False)
    
    # Add exon features
    for cdsFeature in parentFeature.CDS:
        idPrefix, idSuffix = cdsFeature.ID.rsplit(".", maxsplit=1)
        newExonID = f"{idPrefix}.exon{idSuffix[3:]}"
        newExonFeature = GFF3Feature(newExonID, "exon", start=cdsFeature.start, end=cdsFeature.end,
                                     strand=cdsFeature.strand, contig=cdsFeature.contig,
                                     source=cdsFeature.source, score=cdsFeature.score,
                                     frame=cdsFeature.frame, attributes=cdsFeature.attributes,
                                     parents=cdsFeature.parents)
        parentFeature.add_child(newExonFeature)
    
    # Merge gene parent feature into GFF3
    gff3.update_feature(newGeneParent, newID, merging=True)
    
    return newGeneParent

def gff3_mp_reformat(args):
    gff3 = GFF3Tarium(args.gff3File)
    
    identifiers = {}
    for parentFeature in gff3.parents:
        _reformat_miniprot(gff3, identifiers, parentFeature) # the returned geneFeature is unneeeded
    
    gff3.write(args.outputFileName)

def gff3_mp_resolve(args):
    # Parse GFF3 with NCLS indexing
    gff3 = GFF3Tarium(args.gff3File)
    gff3.create_ncls_index(typeToIndex=list(gff3.parentFtypes))
    
    # Get list of parent features ordered from best (identity > length) to worst
    orderedParents = sorted(
        [
            (parentFeature.ID, float(parentFeature.attributes["Identity"]), parentFeature.length("CDS"))
            for parentFeature in gff3.parents
        ], key = lambda x: (-x[1], -x[2])
    )
    
    # Resolve overlaps
    toWrite = []
    handled = set()
    identifiers = {}
    for parentID, identity, length in orderedParents:
        parentFeature = gff3[parentID]
        
        # Skip if this has been handled
        if parentID in handled:
            continue
        
        # Reformat this miniprot mRNA
        newGeneParent = _reformat_miniprot(gff3, identifiers, parentFeature)
        
        # Note this feature for writing to file
        toWrite.append(newGeneParent.ID) # we don't write the features that end up in the handled set
        
        # Note any overlaps as being handled
        overlaps = gff3.ncls_finder(parentFeature.start, parentFeature.end, "contig", parentFeature.contig)
        handled.update([
            o.ID
            for o in overlaps
            if o.strand == newGeneParent.strand
        ])
    
    # Write ordered output to file
    gff3.write(args.outputFileName, idsToWrite=toWrite)

def gff3_to_fasta(args):
    fasta = FASTATarium(args.fastaFile)
    gff3 = GFF3Tarium(args.gff3File)
    
    warnedOnce = False
    with write_conditionally(args.outputFileNames["exon"]) as exonOut, write_conditionally(args.outputFileNames["cds"]) as cdsOut, write_conditionally(args.outputFileNames["protein"]) as protOut:
        for featureType in args.features:
            if not featureType in gff3.ftypes:
                raise ValueError(f"'{featureType}' not found within '{args.gff3File}'")
            
            for parentFeatureID in gff3.ftypes[featureType]:
                # Pick out the longest representative for this feature
                feature = gff3.longest_feature(parentFeatureID)
                
                # Make sure we can produce each requested sequence type for this feature
                if "exon" in args.types:
                    # Prevent sequence output if we cannot create matching exon/CDS/protein files
                    if "cds" in args.types or "protein" in args.types:
                        if not hasattr(feature, "exon"):
                            raise ValueError(f"'{feature.ID}' cannot produce both 'exon' and 'CDS' sequences as it lacks 'exon' children")
                        if not hasattr(feature, "CDS"):
                            raise ValueError(f"'{feature.ID}' cannot produce both 'exon' and 'CDS' sequences as it lacks 'CDS' children")
                        if not hasattr(feature, "exon"):
                            raise ValueError(f"'{feature.ID}' cannot produce both 'exon' and 'CDS' sequences as it lacks 'exon' children")
                    # Allow skipping of exon output if it would not affect matching of files
                    elif not hasattr(feature, "exon"): # i.e., we only want 'exon' output but we can't produce it
                        if not warnedOnce:
                            print(f"WARNING: '{feature.ID}' is a '{featureType}' feature but it lacks 'exon' children; " +
                                  "this and any similar sequences will be omitted from output")
                            warnedOnce = True
                        continue
                
                elif "cds" in args.types or "protein" in args.types:
                    if not hasattr(feature, "CDS"): # i.e., we only want 'CDS' or 'protein' output but we can't produce it
                        if not warnedOnce:
                            print(f"WARNING: '{feature.ID}' is a '{featureType}' feature but it lacks 'CDS' children; " +
                                  "this and any similar sequences will be omitted from output")
                            warnedOnce = True
                        continue
                
                # Format the sequence(s)
                if "exon" in args.types:
                    exonSequence = feature.as_gene_model(fasta, "exon", variantsObj=gff3.get_variants())
                    exonSequence.description = parentFeatureID # propagate the original ID, not the representative
                else:
                    exonSequence = None
                
                if "cds" in args.types or "protein" in args.types:
                    cdsSequence = feature.as_gene_model(fasta, "CDS", variantsObj=gff3.get_variants())
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
        raise ValueError(f"-forEach value '{args.forEach}' is not a feature type in your GFF3")
    
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
