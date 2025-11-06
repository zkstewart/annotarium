#! python3

import intervaltree
from copy import deepcopy

class Coordinates:
    @staticmethod
    def isOverlapping(start1, end1, start2, end2):
        """Does the range (start1, end1) overlap with (start2, end2)?"""
        return end1 >= start2 and end2 >= start1
    
    @staticmethod
    def cleanCoords(coordsTuples):
        """intervaltree won't handle ranges that rub shoulders"""
        return [ (start, end) for start, end in coordsTuples if start != end ]
    
    def __init__(self, coordsTuples):
        self._tree = intervaltree.IntervalTree.from_tuples(Coordinates.cleanCoords(coordsTuples))
        self._tree.merge_overlaps()
        self.coords = sorted( (interval.begin, interval.end) for interval in self._tree )
    
    def intersection(self, otherCoordinates):
        overlaps = [
            (max(start1, start2), min(end1, end2))
            for start1, end1 in self.coords
            for start2, end2 in otherCoordinates.coords
            if Coordinates.isOverlapping(start1, end1, start2, end2)
        ]
        overlaps = [ (start, end) for start, end in overlaps if start != end ] 
        
        return Coordinates(overlaps)
    
    def overlap_length(self, otherCoordinates):
        return len(self.intersection(otherCoordinates))
    
    def overlap_percent(self, otherCoordinates):
        overlapLen = self.overlap_length(otherCoordinates)
        if overlapLen == 0:
            return 0.0, 0.0
        
        thisLen = len(self)
        otherLen = len(otherCoordinates)
        thisOverlapPct = (thisLen - (thisLen - overlapLen)) / thisLen
        otherOverlapPct = (otherLen - (otherLen - overlapLen)) / otherLen
        return thisOverlapPct, otherOverlapPct
    
    def __iter__(self):
        yield from self.coords
    
    def __len__(self):
        return sum( end - start + 1 for start, end in self.coords )
    
    def __repr__(self):
        return f"<Coordinates object; coords={list(self.coords)}"

class OverlapResolver:
    '''
    Notes for future development:
    - this is legacy code with some minor refactoring to accept objects with .id/.start/.end/.score rather
      than a list of those values; it is ugly at times but it has seen a lot of use (and hence testing) and
      scrapping the whole codebase seems wasteful.
    - sharedPos and the whole system of doing set(range()) is incredibly inefficient and should be
      refactored if this code is used for any heavy lifting.
    
    Attributes:
        ovlCutoff -- a float or int as a ratio from 0 to 1 (inclusive) representing a percentage value.
                     Controls how domain overlapping is handling; values up to ovlCutoff %
                     are handled by trimming, and values exceeding ovlCutoff % have the lower score
                     feature culled.
        minScore -- a float or int representing the minimum score a feature must have; values
                    lower than this are culled. Or give None to prevent any filtering during overlap
                    resolution.
    '''
    def __init__(self, ovlCutoff=0.25, minScore=None):
        self.ovlCutoff = ovlCutoff
        self.minScore = minScore
        
        # Set helper attribute
        self.isOverlapResolver = True
    
    @staticmethod
    def validate_feature_list(featureList):
        if not isinstance(featureList, list):
            raise TypeError(f"featureList must be a 'list', not '{type(featureList).__name__}'")
        for feature in featureList:
            try:
                feature.start
                feature.end
                feature.score
            except:
                raise TypeError(f"featureList must contain values with .start, .end. and .score attributes; " +
                                f"{feature} lacks these")
    
    @staticmethod
    def find_middle(inputList):
        '''
        See https://stackoverflow.com/questions/38130895/find-middle-of-a-list
        '''
        middle = float(len(inputList))/2
        if middle % 2 != 0:
            return [inputList[int(middle - .5)]]
        else:
            return [inputList[int(middle)], inputList[int(middle-1)]]
    
    @staticmethod
    def split_middle(sharedPos, featureGroup, y):
        '''
        Updates overlapping features by setting their end (leftmost feature) and
        start (rightmost feature) values to be adjacent down their centre.
        Essentially, it splits the difference of an overlapping region as fairly
        as possible.
        '''
        splitPos = list(sharedPos)
        splitPos.sort()
        middle = OverlapResolver.find_middle(splitPos)
        if len(middle) == 1:
            featureGroup[y].end = middle[0]
            featureGroup[y+1].start = middle[0]+1
        else:
            featureGroup[y].end = middle[1]
            featureGroup[y+1].start = middle[0]
        return featureGroup
    
    @staticmethod
    def join_features(featureGroup, y):
        '''
        "Merges" two features together (occurring at index y and y+1) by joining their
        coordinate ranges and maintaining the worse scoring value. Retaining the worse
        score is important to ensure it does not "domino" its score up by consuming
        other features.
        '''
        featureGroup[y].end = max(featureGroup[y].end, featureGroup[y+1].end)
        featureGroup[y].score = max(featureGroup[y].score, featureGroup[y+1].score)
        del featureGroup[y+1]
        return featureGroup
    
    @property
    def ovlCutoff(self):
        return self._ovlCutoff
    
    @ovlCutoff.setter
    def ovlCutoff(self, value):
        if not isinstance(value, float) and not isinstance(value, int):
            raise TypeError(f"ovlCutoff must be a 'float' or 'int', not '{type(value).__name__}'")
        
        if value < 0:
            raise ValueError("ovlCutoff must be a value >= 0")
        if value > 1:
            raise ValueError("ovlCutoff must be a value <= 1")
        
        self._ovlCutoff = value
    
    @property
    def minScore(self):
        return self._minScore
    
    @minScore.setter
    def minScore(self, value):
        if value is None:
            self._minScore = None
            return
        
        if not isinstance(value, float) or isinstance(value, int):
            raise TypeError(f"minScore must be a 'float' or 'int', not '{type(value).__name__}'")
        
        self._minScore = value
    
    def _ovl_resolver(self, featureList):
        '''
        Parameters:
            featureList -- a list of any Feature-type objects, which minimally have the following attributes:
                id -- a string labelling this feature e.g., a PFAM domain identifier
                start -- an integer indicating the 1-based start position
                end -- an integer indicating the 1-based (inclusive) end position
                score -- an integer or float where LOWER values are better e.g., E-values
        '''
        OverlapResolver.validate_feature_list(featureList)
        featureList = deepcopy(featureList)
        
        featureList.sort(key = lambda x: (x.score, x.start, x.end))
        while True:
            # Flow control: if no possibility for overlap, break
            if len(featureList) == 1:
                break
            
            # Flow control: if no overlaps remain, break
            overlapping = False
            for y in range(len(featureList)-1):
                feature1 = featureList[y]
                for z in range(y+1, len(featureList)):
                    feature2 = featureList[z]
                    if Coordinates.isOverlapping(feature1.start, feature1.end, feature2.start, feature2.end):
                        overlapping = True
                        break
            if not overlapping:
                break
            
            # Resolve overlaps through looping structure
            featureList = self._seed_looping_structure(featureList)
        
        featureList.sort(key = lambda x: (x.start, x.end))
        return featureList
    
    def _seed_looping_structure(self, featureList):
        # Set up
        origFeatureList = deepcopy(featureList) # This lets us compare our new domain against the original to make sure we haven't excessively cut and trimmed it
        ORIG_CUTOFF = 0.60  # Arbitrary; this means that, if the trimmed feature is less than "a bit above" the length of the original, we drop it entirely
        # Main loop
        for y in range(len(featureList)-1):
            z = y + 1
            while True:
                # Exit condition
                if z >= len(featureList): # len(featureList)-1 would correspond to the final entry, this means we've gone at least one step further beyond
                    break
                
                # Get the features being compared
                feature1 = featureList[y]
                origFeature1 = origFeatureList[y]
                
                feature2 = featureList[z]
                origFeature2 = origFeatureList[z]
                
                # Trim 1-bp overlaps [note the consistent design choice: the most N-proximal domain gets the extra AA position for no particular reason, we just need to handle this]
                if feature1.end == feature2.start:
                    feature2.start =  feature2.start + 1
                if feature2.end == feature1.start:
                    feature1.start =  feature1.start + 1
                
                # If there is overlap, resolve this
                if Coordinates.isOverlapping(feature1.start, feature1.end, feature2.start, feature2.end):
                    # Get details of sequence overlap
                    sharedPos = set(range(max(feature1.start, feature2.start), min(feature1.end, feature2.end) + 1))
                    ovlLen = len(sharedPos)
                    feature1Perc = ovlLen / (feature1.end - feature1.start + 1)
                    feature2Perc = ovlLen / (feature2.end - feature2.start + 1)
                    bestScore = min(feature1.score, feature2.score)
                    
                    # Handle slight mutual overlaps by trimming based on lowest/best score
                    if feature1Perc < self.ovlCutoff and feature2Perc < self.ovlCutoff:
                        ## Identical scores [ mutual trimming ]
                        if feature1.score == feature2.score:
                            posList = list(sharedPos)
                            posList.sort()
                            midPoint = OverlapResolver.find_middle(posList)
                            if feature1.start < feature2.start:
                                feature1.end = midPoint[0]
                                feature2.start = midPoint[0] + 1
                            else:
                                feature2.end = midPoint[0]
                                feature1.start = midPoint[0] + 1
                        
                        ## Different scores [ trim higher/worse scoring feature ]
                        elif bestScore == feature1.score:
                            if feature1.start < feature2.start:
                                feature2.start = feature1.end + 1
                            else:
                                feature2.end = feature1.start - 1
                        else:
                            if feature1.start < feature2.start:
                                feature1.end = feature2.start - 1
                            else:
                                feature1.start = feature2.end + 1
                        
                        # Validate program behaviour; assure that domains being compared are still equivalent
                        assert origFeature1.id == feature1.id
                        assert origFeature2.id == feature2.id
                        
                        # If we've trimmed one of these domains too much, drop it
                        changed = False
                        if (feature2.end - feature2.start + 1) / (origFeature2.end - origFeature2.start + 1) < ORIG_CUTOFF: # Need to handle z first lest we upset the ordering
                            del featureList[z] # feature2
                            del origFeatureList[z]
                            changed = True
                        if (feature1.end - feature1.start + 1) / (origFeature1.end - origFeature1.start + 1) < ORIG_CUTOFF:
                            del featureList[y] # feature1
                            del origFeatureList[y]
                            changed = True
                        if not changed:
                            z += 1  # We've made the current pair compatible, now we can just move onto the next pairing
                    
                    # Handle larger overlaps by deleting based on score
                    else:
                        ## Identical scores [ delete the most C-proximal ]
                        if feature1.score == feature2.score:
                            if feature1.start < feature2.start:
                                del featureList[z] # feature2
                                del origFeatureList[z] # keep the original and modified lists equivalent
                            else:
                                del featureList[y] # feature1
                                del origFeatureList[y]
                        
                        ## Different scores [ delete the higher/worse score ]
                        elif bestScore == feature1.score:
                            del featureList[z] # feature2
                            del origFeatureList[z]
                        else:
                            del featureList[y] # feature1
                            del origFeatureList[y]
                        
                        # We make no changes to our z value since we deleted a sequence
                
                # If there is no overlap, continue the loop by changing our z value
                else:
                    z += 1
        return featureList
    
    def resolve(self, featureList):
        '''
        Represents the main entry point to resolving the overlap of features by coordinates
        (.start, .end) and score (minimal == best). The logic was initially developed for handling
        HMMER domain predictions, but the idea can apply to any non-unique region with a score
        attached to it.
        '''
        EXTENSION_CUTOFF = 20 # This is arbitrary; seems to work well, don't see any reason why this should be variable by the user
        OverlapResolver.validate_feature_list(featureList)
        
        # Collapse overlaps of identical domains
        resolvedFeatureList = []
        uniqueFeatureIDs = list(set([ f.id for f in featureList ]))
        for groupID in uniqueFeatureIDs:
            # Obtain all features that have the same identifier
            featureGroup = []
            for feature in featureList:
                if feature.id == groupID:
                    featureGroup.append(feature)
            featureGroup.sort(key = lambda x: (x.start, x.end))  # Technically this should not be needed - the HMMER domtblout file is pre-sorted - but it's useful to put here _just in case_, and to make it clear that this script operates on the basis of this sorting
            
            # Begin collapsing process
            overlapping = True
            while True:
                # Exit condition
                if len(featureGroup) == 1 or overlapping == False:
                    break
                
                # Iterate over features
                for y in range(len(featureGroup)-1):
                    feature1 = featureGroup[y]
                    feature2 = featureGroup[y+1]
                    
                    # Skip non-overlaps and the last pair of comparisons
                    ## Note: don't understand the logic here as of 06/11/2025
                    if feature2.start > feature1.end and y != len(featureGroup)-2: # we want to skip this if it's the last pair since that will allow us to reach the final "else" condition and exit out of the loop
                        continue
                    # Quickly resolve 1bp overlaps
                    elif feature2.start == feature1.end:
                        feature2.start =  feature2.start + 1 # Consistent design choice: the most N-proximal domain gets the extra AA position for no particular reason, we just need to do _something_
                        continue
                    # Resolve larger overlaps
                    elif feature2.start < feature1.end:
                        # Calculate overlap proportion
                        feature1Len = feature1.end - feature1.start + 1
                        feature2Len = feature2.end - feature2.start + 1
                        sharedPos = set(range(max(feature1.start, feature2.start), min(feature1.end, feature2.end) + 1)) # +1 to offset Python counting up-to but not including the last value in a range
                        ovlLen = len(sharedPos)
                        r1Perc = ovlLen / (feature1Len + 1)
                        r2Perc = ovlLen / (feature2Len + 1)
                        highest = max(r1Perc, r2Perc)
                        lowest = min(r1Perc, r2Perc)
                        
                        # Determine the length of the sequence extension of the most-overlapped sequence
                        if highest == 0.50: # magic number
                            longest = max(feature1Len, feature2Len)
                            if longest == feature1Len:
                                extension = feature2Len - ovlLen
                            else:
                                extension = feature1Len - ovlLen
                        elif highest == r1Perc:
                            extension = feature1Len - ovlLen
                        else:
                            extension = feature2Len - ovlLen
                        
                        # Handle scenario 1: small overlap of both sequences (TRIM TO MINIMISE SCORE) 
                        if highest <= 0.20: # magic number
                            if feature1.score < feature2.score:
                                # Trim y+1
                                feature2.start = feature1.end+1
                            elif feature2.score < feature1.score:
                                # Trim y
                                feature1.end = feature2.start-1
                            else:
                                # If the two E-value are identical, we just split down the middle!
                                featureGroup = OverlapResolver.split_middle(sharedPos, featureGroup, y)
                            continue
                        
                        # Handle scenario 2: intermediate overlap with significant sequence extension beyond the overlap region (SPLIT MIDDLE)
                        elif extension > EXTENSION_CUTOFF and lowest <= 0.80: # magic number
                            featureGroup = OverlapResolver.split_middle(sharedPos, featureGroup, y)
                            continue
                        
                        # Handle scenario 3: >= intermediate overlap with a short sequence extension beyond the overlap region (JOIN)
                        else:
                            featureGroup = OverlapResolver.join_features(featureGroup, y)
                            break
                    
                    # Trigger exit condition if no overlaps remain
                    else: # We need the y != check above since we need to set an exit condition when no more overlaps are present. The if/elif will always trigger depending on whether there is/is not an overlap UNLESS it's the second last entry and there is no overlap. In this case we finally reach this else clause, and we trigger an exit.
                        overlapping = False
                        break
            
            # Add corrected individual features to resolvedFeatureList list
            resolvedFeatureList += featureGroup
        
        # Resolve non-identical domain predictions
        if len(resolvedFeatureList) != 1:
            resolvedFeatureList = self._ovl_resolver(resolvedFeatureList) # we've merged, joined, and trimmed identical features above. Now, we're looking at different domains.
        
        return resolvedFeatureList
