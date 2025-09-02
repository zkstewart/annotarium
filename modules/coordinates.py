#! python3

import intervaltree

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
