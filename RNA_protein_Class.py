


class RNASequence:
    CODON_TABLE = {
    "AUG": "Met",

    "UUU": "Phe", "UUC": "Phe",
    "UUA": "Leu", "UUG": "Leu",

    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",

    "UAU": "Tyr", "UAC": "Tyr",
    "UAA": "STOP", "UAG": "STOP",

    "UGU": "Cys", "UGC": "Cys",
    "UGA": "STOP",
    "UGG": "Trp",

    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",

    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",

    "CAU": "His", "CAC": "His",
    "CAA": "Gln", "CAG": "Gln",

    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",

    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile",

    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",

    "AAU": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys",

    "AGU": "Ser", "AGC": "Ser",
    "AGA": "Arg", "AGG": "Arg",

    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",

    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",

    "GAU": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu",

    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"}
    def __init__(self,sequence: str):
                
        self.__sequence = sequence
        self.validate_seq() 
        
    def validate_seq(self) -> None: 
        seq = self.__sequence.upper().replace(" ", "").replace('\n', "")
        valid_bases = {'A', 'U', 'G', 'C'}

        if len(seq) == 0:
            raise ValueError('Sequence can\'t be empty') 
            
        for base in seq:
            if base not in valid_bases:
                raise ValueError(f'{base} not in {valid_bases}')
            
    
        self.__sequence = seq 
     
    def get_codons(self) -> list[str]:
        seq = self.__sequence
        
        return [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
    
    def translate(self) -> list[str]:
        codons = self.get_codons()
        protein = [ ]
        
        for codon in codons:
            amino = RNASequence.CODON_TABLE.get(codon)
            
            if amino == 'STOP':
                break
            
            protein.append(amino)
            
        return protein 
        
    
    def get_start_codon(self):
        return self.__sequence.find('AUG')
    
    def has_start_codon(self):
        if self.get_start_codon() != -1:
        
            return True
        
        return False 
    def has_stop_codon(self):
        if self.get_stop_codon() != -1:
            
            return True
        
        return False
    
    def get_stop_codon(self):
        seq = self.__sequence
        stop_codons = {'UAA', 'UAG', 'UGA'}
        
        if self.has_start_codon():
            for i in range(self.get_start_codon(), len(seq), 3):
                codon = seq[i:i+3]
                
                if len(codon) < 3 : 
                    break
                if codon in stop_codons:
                    raise Exception('Missense mutation found')  
                    return i
        return i     


        
    def get_coding_region(self):
        
        start = self.get_start_codon()
        
        stop = self.get_stop_codon()
        
        if not start or not stop:
            Exception('No coding region') 
            return ""
        return self.__sequence[start:stop]
    
    def amino_acids_count(self):
        protein = self.translate()
        amino_dict = {} 
        for amino in protein:
            if amino not in amino_dict:
                amino_dict[amino] = 1
            else:
                amino_dict[amino] += 1
        return amino_dict
    
    def predict_protein_functionality(self):
        
        protein = self.translate()
        
        if self.get_start_codon() == False and self.get_stop_codon() == False and len(portein) < 30:
            return 'This sequence is not producing a functional protein'
        
        elif self.get_start_codon() == False :
            return 'Start codon not found, RNA sequence is possibly frameshited'
        
        elif self.get_stop_codon() == False:
            return 'Can not find stop codon, sequence could be functional'

        elif len(protein) > 30:
            return 'Protein is less than 30 bases likely to be small peptide chain or nonfunctional protein.'
        
        
        return 'Sequence contains start, stop, and meets basic length requirements\nwill likely be a functional protein.'
    
    
        
            
    
rna = RNASequence('AUGGCCUACGUGGAGGCUAA')

print(rna.translate())
print(rna.amino_acids_count())
print(rna.has_start_codon())

print(rna.get_stop_codon()) 
    
    
        
        
        
        