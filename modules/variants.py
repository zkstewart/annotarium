#! python3

class Variants:
    def __init__(self):
        self.variants = {}
        self.isVariants = True
    
    def add(self, variant):
        if not hasattr(variant, "isVariant") and variant.isVariant:
            return TypeError(f"Variant class cannot add '{type(variant).__name__}' type object")
        self.variants.setdefault(variant.contig, set())
        self.variants[variant.contig].add(variant)
    
    def __setitem__(self, idx, value):
        self.variants[idx] = value
    
    def __getitem__(self, value):
        return sorted(self.variants[value], key = lambda x: (x.start, x.end))
    
    def __contains__(self, value):
        if isinstance(value, str):
            return value in self.variants
        
        if not hasattr(value, "isVariant") and not value.isVariant:
            return False
        elif value.contig in self.variants:
            return value in self.variants[value.contig]
        else:
            return False
    
    def __len__(self):
        return len(self.variants)
    
    def __str__(self):
        return str(self.variants)
    
    def __repr__(self):
        return "<Variants object;num_contigs={0};num_variants={1}".format(
            len(self.variants),
            sum([ len(v) for v in self.variants.values() ])
        )

class Variant:
    '''
    Class is expected to follow VCF conventions such that positions are 1-based and inclusive
    e.g., start==5 and end==5 points to the fifth nucleotide along a sequence.
    '''
    EDIT_TYPES = ["substitution", "insertion", "deletion"]
    
    def __init__(self, contig, start, end, editType, residues):
        self.contig = contig
        self.start = start
        self.end = end
        self.editType = editType
        self.residues = residues
        self.isVariant = True
    
    @property
    def start(self):
        return self._start
    
    @start.setter
    def start(self, value):
        if not isinstance(value, int):
            raise TypeError("Variant .start must be an int")
        if value == 0:
            raise ValueError("Variant .start must be 1-based; 0 is not a valid position")
        if value < 0:
            raise ValueError("Variant .start must be a positive integer")
        self._start = value
     
    @property
    def end(self):
        return self._end
    
    @end.setter
    def end(self, value):
        if not isinstance(value, int):
            raise TypeError("Variant .end must be an int")
        if value == 0:
            raise ValueError("Variant .end must be 1-based; 0 is not a valid position")
        if value < 0:
            raise ValueError("Variant .end must be a positive integer")
        self._end = value
    
    @property
    def editType(self):
        return self._editType
    
    @editType.setter
    def editType(self, value):
        if not isinstance(value, str):
            raise TypeError("Variant .editType must be a string")
        if not value in Variant.EDIT_TYPES:
            raise ValueError(f"Variant .editType must be one of {Variant.EDIT_TYPES}")
        self._editType = value
    
    @property
    def residues(self):
        return self._residues
    
    @residues.setter
    def residues(self, value):
        if value == None:
            if not self.editType == "deletion":
                raise TypeError("Variant .residues cannot be None if .editType is not 'deletion'")
        elif not isinstance(value, str):
            raise TypeError("Variant .residues must be a string")
        self._residues = value
    
    def __hash__(self):
        return hash((self.contig, self.start, self.end, self.editType, self.residues))
    
    def __eq__(self, otherValue):
        if hasattr(otherValue, "isVariant"):
            return hash(self) == hash(otherValue)
        else:
            return False
    
    def __repr__(self):
        return f"<Variant object;start={self.start};end={self.end};editType={self.editType};residues={self.residues}>"
