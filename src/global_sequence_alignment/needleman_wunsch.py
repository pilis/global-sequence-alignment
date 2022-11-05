from typing import Dict, List, Set


class ScoringFunction:
    """Abstract class for scoring functions"""

    def __init__(self, gap_penalty: int):
        self.gap_penalty = gap_penalty

    def score(self, gap_lenght: int) -> int:
        raise NotImplementedError


class LinearGapPenalty(ScoringFunction):
    """Linear gap penalty scoring function"""

    def __init__(self, gap_penalty: int):
        super().__init__(gap_penalty)
        pass

    def score(self, gap_length: int) -> int:
        return gap_length * self.gap_penalty


class AffineGapPenalty(ScoringFunction):
    """Affine gap penalty scoring function"""

    def __init__(self, gap_penalty: int, gap_extension_penalty: int):
        super().__init__(gap_penalty)
        self.gap_extension_penalty = gap_extension_penalty

    def score(self, gap_length: int) -> int:
        return self.gap_penalty + (gap_length - 1) * self.gap_extension_penalty


class InvalidSymbolError(Exception):
    """Exception raised when an invalid symbol is found"""

    pass


class SubstitutionMatrix:
    def __init__(self, transitions: Dict[str, Set[str]]):
        self.transitions = transitions

    def is_equal(self, symbol_1: str, symbol_2: str) -> bool:
        if symbol_1 not in self.transitions:
            raise InvalidSymbolError(
                f"Symbol {symbol_1} is not in the substitution matrix"
            )
        if symbol_2 not in self.transitions:
            raise InvalidSymbolError(
                f"Symbol {symbol_2} is not in the substitution matrix"
            )
        substitutes = self.transitions[symbol_1]
        return symbol_2 in substitutes


NUCLEOTIDE_SUBSTITUTIONS = {"A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"}}


class NucleotideSubstitutionMatrix(SubstitutionMatrix):
    """Nucleotide substitution matrix"""

    def __init__(self):
        super().__init__(transitions=NUCLEOTIDE_SUBSTITUTIONS)
        pass


# TODO: Implement substitution matrix for proteins


class Alignment:
    def __init__(self):
        pass

    # TODO: Think how to model alignment between two sequences


class ScoringMatrix:
    def __init__(
        self,
        sequence_1: str,
        sequence_2: str,
        scoring_function: ScoringFunction,
        substitution_matrix: SubstitutionMatrix,
    ):
        self.sequence_1 = sequence_1
        self.sequence_2 = sequence_2
        self.scoring_function = scoring_function
        self.substitution_matrix = substitution_matrix

        self.scoring_matrix = self._init_scoring_matrix()

    def _init_scoring_matrix(self) -> List[List[None]]:  # type: ignore
        """Initialize 2D matrix for holding scores"""
        horizontal_length = len(self.sequence_1) + 1
        vertical_length = len(self.sequence_2) + 1

        # Initialize matrix with None values
        scoring_matrix = [[None] * horizontal_length for _ in range(vertical_length)]

        gap_penalty = self.scoring_function.gap_penalty
        # Initialize first row with gap penalties times index
        for i in range(horizontal_length):
            scoring_matrix[0][i] = i * gap_penalty  # type: ignore
        # Initialize first column with gap penalties times index
        for j in range(vertical_length):
            scoring_matrix[j][0] = j * gap_penalty  # type: ignore

        return scoring_matrix

    def fill(self):
        """Fill 2D matrix with scores"""
        pass

    def traceback(self) -> List[Alignment]:
        """Traceback 2D matrix to find optimal alignments"""
        pass


class NeedlemanWunsch:
    def __init__(
        self, scoring_function: ScoringFunction, substitution_matrix: SubstitutionMatrix
    ) -> None:
        self.scoring_function = scoring_function
        self.substitution_matrix = substitution_matrix

    def align(self, sequence_1, sequence_2) -> List[Alignment]:
        """Align two sequences using the Needleman-Wunsch algorithm"""
        # scoring_matrix = ScoringMatrix(sequence_1, sequence_2)
        # return scoring_matrix.traceback()
        pass
