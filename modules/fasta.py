#! python3

import os
from Bio.Data import CodonTable

from .parsing import read_gz_file

class Sequence:
    COMPLEMENT = {"A": "T", "a": "t",
                  "T": "A", "t": "a",
                  "C": "G", "c": "g",
                  "G": "C", "g": "c"}
    '''
    Parameters:
        description -- a string with or without whitespace
        sequence -- a string of nucleotides or amino acids
        frame -- (OPTIONAL) an integer used for translation if relevant
    '''
    def __init__(self, description, sequence, frame=None):
        self.description = description
        self.sequence = sequence
        self.frame = frame
        self.isSequence = True
    
    @property
    def description(self):
        return self._description
    
    @description.setter
    def description(self, value):
        if not isinstance(value, str):
            raise TypeError("description must be a string")
        value = value.strip(">\r\n\t ")
        self._description = value
        
        # Format and set self.id
        self._id = value.split(" ")[0]
    
    @property
    def id(self):
        return self._id
    
    @property
    def sequence(self):
        return self._sequence
    
    @sequence.setter
    def sequence(self, value):
        if not isinstance(value, str):
            raise TypeError("sequence must be a string")
        self._sequence = value
    
    @property
    def frame(self):
        if self._frame == None:
            return 0
        return self._frame
    
    @frame.setter
    def frame(self, value):
        if value == None or value == ".":
            self._frame = None
        else:
            try:
                value = int(value)
            except:
                raise TypeError(f"frame value '{value}' is neither an integer nor '.'")
            
            if value < 0 or value > 2:
                raise ValueError("frame must be in the range of 0->2 (inclusive)") 
            self._frame = value
    
    def reverse_complement(self):
        '''
        Returns the reverse complement of this sequence. The Sequence class
        is naive to what type of molecule is contains, so you should ensure you
        only run this on nucleotides.
        '''
        reverseSequence = "".join(Sequence.COMPLEMENT.get(nuc, nuc) for nuc in reversed(self.sequence))
        return Sequence(self.description, reverseSequence, self.frame)
    
    def translate(self, tableNum=1):
        '''
        Translates this sequence into a protein. It is up to you to make this sequence is
        a nucleotide, as this function is naive to molecule type.
        
        Parameters:
            tableNum -- (OPTIONAL) an integer for the NCBI codon table to use during translation
        Returns:
            seqObj -- a new Sequence object containing this object's sequence, but translated
                      into a protein according to the translation table indicated
        '''
        # Get our codons table for translation
        try:
            table = CodonTable.ambiguous_dna_by_id[tableNum]
            codons = table.forward_table.forward_table
            codons.update({ codon:"*" for codon in table.stop_codons})
        except:
            raise ValueError(f"Translation table '{tableNum}' is not recognised by Biopython as a valid table number")
        
        # Translate the sequence
        protein = []
        for i in range(self.frame, len(self.sequence), 3):
            codon = self.sequence[i:i+3]
            protein.append(codons.get(codon.upper(), "X"))
        
        # Return the translation as a new object
        return Sequence(self.description, "".join(protein), self.frame)
    
    def format(self, lineLength=None):
        '''
        Returns this Sequence object in FASTA format for e.g., writing to file.
        
        Parameters:
            lineLength -- an integer indicating how many characters per line (for
                          multiline FASTA format) or None to have one line per sequence
        Returns:
            fastaString -- a FASTA-formatted string ending with newline
        '''
        sequence = self.sequence
        if lineLength != None:
            sequence = "\n".join([sequence[i:i+lineLength] for i in range(0, len(sequence), lineLength)])
        return f">{self.description}\n{sequence}\n"
    
    def __getitem__(self, value):
        if isinstance(value, int):
            return self.sequence[value] # just return the individual residue as a string
        elif isinstance(value, slice):
            return Sequence(self.description, self.sequence[value]) # slice returns a new object
    
    def __str__(self):
        return self.sequence
    
    def __repr__(self):
        return "<Sequence object;description='{0}';id='{1}';sequence='{2}'>".format(
            self.description,
            self.id,
            self.sequence[0:50] + "..." if len(self.sequence) > 51 else self.sequence
        )

class Records:
    def __init__(self):
        self.dict = {}
        self.isRecords = True
    
    def add(self, sequenceObj):
        if not hasattr(sequenceObj, "id"):
            raise TypeError(f"Cannot index '{type(sequenceObj)}' as a FASTA record")
        if sequenceObj.id in self.dict:
            raise ValueError(f"Disallowed attempt to index a duplicate '{sequenceObj.id}' FASTA record")
        self.dict[sequenceObj.id] = sequenceObj
    
    def __getitem__(self, key):
        return self.dict[key]
    
    def __len__(self):
        return len(self.dict)
    
    def __iter__(self):
        return iter(self.dict.values())
    
    def __contains__(self, value):
        return value in self.dict
    
    def keys(self):
        return self.dict.keys()
    
    def values(self):
        return self.dict.values()
    
    def items(self):
        return self.dict.items()
    
    def __repr__(self):
        return "<Records object;num_records={0}>".format(
            len(self.dict)
        )

class FastaFormatError(Exception):
    pass
class IncompleteRecordError(Exception):
    pass

class FASTATarium:
    @staticmethod
    def fasta_parser(handle):
        '''
        Parameters:
            handle -- an active file handle to read from
        Yields:
            seqObj -- a Sequence object as defined in this file
        '''
        description = None
        sequence = []
        for line in handle:
            # Handle first iteration
            if description is None:
                if not line.startswith(">"):
                    raise FastaFormatError(f"'{handle.name}' does not begin with the expected '>' character")
                else:
                    description = line[1:].strip()
                    continue
            
            # Handle description line
            if line.startswith(">"):
                if len(sequence) == 0:
                    raise IncompleteRecordError(f"'{description}' has no associated sequence")
                yield Sequence(description, "".join(sequence))
                
                description = line[1:].strip()
                sequence = []
            # Handle sequence line
            else:
                sequence.append(line.strip())
        
        # Handle end of file
        if len(sequence) == 0:
            raise IncompleteRecordError(f"'{description}' has no associated sequence")
        yield Sequence(description, "".join(sequence))
    
    def __init__(self, fastaFile):
        self.records = Records()
        self.fastaFile = fastaFile
        self.isFASTATarium = True
    
    @property
    def fastaFile(self):
        return self._fastaFile
    
    @fastaFile.setter
    def fastaFile(self, value):
        if not isinstance(value, str):
            raise TypeError("fastaFile must be a string")
        if not os.path.isfile(value):
            raise FileNotFoundError(f"FASTA '{value}' is not a file.")
        
        self._fastaFile = value
        self.parse_file()
    
    def parse_file(self):
        self.records = Records()
        with read_gz_file(self.fastaFile) as fileIn:
            for seqObj in FASTATarium.fasta_parser(fileIn):
                self.add(seqObj)
    
    def add(self, seqObj):
        self.records.add(seqObj)
    
    def __getitem__(self, value):
        if hasattr(value, "isGFF3Feature") and value.isGFF3Feature:
            return self(value.contig, value.start-1, value.end) # counteract 1-based GFF3 numbering
        else:
            return self.records[value]
    
    def __call__(self, key, start, end):
        return self[key][start:end]
    
    def __len__(self):
        return len(self.records)
    
    def __iter__(self):
        return iter(self.records.values())
    
    def __contains__(self, value):
        return value in self.records
    
    def ids(self):
        return self.records.keys()
    
    def sequences(self):
        return self.records.values()
    
    def __repr__(self):
        return "<FASTATarium object;file='{0}';num_records={1}>".format(
            self.fastaFile,
            len(self.records)
        )

def fasta_stats(args):
    return None
