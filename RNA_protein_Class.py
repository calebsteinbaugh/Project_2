class RNASequence:
    """
    Represents an RNA sequence and provides methods for validation,

    codon parsing, translation, and basic protein analysis.
    """
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

        "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
    }

    def __init__(self, sequence: str) -> None:
        """
        Represents an RNA sequence and provides methods for validation,

        codon parsing, translation, and basic protein analysis.
        """
        self.__sequence = sequence
        self.validate_seq()

    def validate_seq(self) -> None:
        """
        Validate and clean the RNA sequence.

        Ensures:

        - Sequence is not empty

        - Only contains A, U, G, C

        Raises:

            ValueError: If sequence is invalid.
        """
        seq = self.__sequence.upper().replace(" ", "").replace("\n", "")
        valid_bases = {"A", "U", "G", "C"}

        if len(seq) == 0:
            raise ValueError("Sequence can't be empty.")

        for base in seq:
            if base not in valid_bases:
                raise ValueError(f"{base} not in {valid_bases}.")

        self.__sequence = seq

    def get_codons(self) -> list[str]:
        """
        Split RNA sequence into codons (triplets).

        Returns:

            list[str]: List of codons.
        """
        seq = self.__sequence
        return [
            seq[i:i + 3]
            for i in range(0, len(seq), 3)
            if len(seq[i:i + 3]) == 3
        ]

    def translate(self) -> list[str]:
         """
        Translate RNA sequence into amino acids until a stop codon.

        Returns:

            list[str]: Amino acid sequence.
        """
        protein = []

        for codon in self.get_codons():
            amino = RNASequence.CODON_TABLE.get(codon)

            if amino == "STOP":
                break

            if amino is not None:
                protein.append(amino)

        return protein

    def get_start_codon(self) -> int:
        """
        Find index of first start codon (AUG).

        Returns:

            int: Index of AUG, or -1 if not found.
        """
        return self.__sequence.find("AUG")

    def has_start_codon(self) -> bool:
        """
        Check if sequence contains a start codon.

        Returns:

            bool: True if AUG exists, False otherwise.
        """
        return self.get_start_codon() != -1

    def get_stop_codon(self) -> int:
         """
        Find index of first stop codon after the start codon.

        Returns:

            int: Index of stop codon, or -1 if not found.
        """
        seq = self.__sequence
        stop_codons = {"UAA", "UAG", "UGA"}

        start_codon = self.get_start_codon()

        if start_codon == -1:
            return -1

        for i in range(start, len(seq), 3):
            codon = seq[i:i + 3]

            if len(codon) < 3:
                break

            if codon in stop_codons:
                return i

        return -1

    def has_stop_codon(self) -> bool:
          """
        Check if sequence contains a valid stop codon.

        Returns:

            bool: True if stop codon exists, False otherwise.
        """
        return self.get_stop_codon() != -1

    def get_coding_region(self) -> str:
         """
        Extract coding region from start codon to stop codon.

        Returns:

            str: Coding region, or empty string if invalid.
        """
        start_codon = self.get_start_codon()
        stop_codon = self.get_stop_codon()

        if start_codon == -1 or stop_codon == -1:
            return ""

        return self.__sequence[start_codon:stop_codon]

    def amino_acids_count(self) -> dict[str, int]:
        """
        Count frequency of amino acids in translated protein.

        Returns:

            dict[str, int]: Amino acid counts.
        """
        protein = self.translate()
        amino_dict = {}

        for amino in protein:
            if amino not in amino_dict:
                amino_dict[amino] = 1
            else:
                amino_dict[amino] += 1

        return amino_dict

    def predict_protein_functionality(self) -> str:
        """
        Provide a basic prediction of protein functionality.

        Returns:

            str: Interpretation of the RNA sequence.
        """
        protein: list[str] = self.translate()

        if not self.has_start_codon():
            return "Start codon not found. RNA sequence is unlikely to translate properly."

        if not self.has_stop_codon():
            return "Stop codon not found. Protein may be incomplete or extended."

        if len(protein) == 0:
            return "No amino acid product detected."

        if len(protein) < 30:
            return "Protein is less than 30 amino acids. It may be a small peptide or nonfunctional fragment."

        return (
            "Sequence contains a start codon, stop codon, and meets basic length requirements.\n"
            "It will likely produce a functional protein."
        )
