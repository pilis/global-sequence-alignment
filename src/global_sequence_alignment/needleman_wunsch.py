import enum
import logging
from typing import List, Tuple, Union

GAP_PENALTY = -1
GAP_EXTENSION_PENALTY = -1


class ScoringFunction:
    """Abstract class for scoring functions"""

    def __init__(self, gap_penalty: int = GAP_PENALTY):
        self.gap_penalty = gap_penalty

    def score(self, gap_lenght: int) -> int:
        raise NotImplementedError


class ConstantGapPenalty(ScoringFunction):
    """Constant gap penalty scoring function"""

    def __init__(self, gap_penalty: int = GAP_PENALTY):
        super().__init__(gap_penalty)

    def score(self, gap_lenght: int) -> int:
        return self.gap_penalty


class LinearGapPenalty(ScoringFunction):
    """Linear gap penalty scoring function"""

    def __init__(self, gap_penalty: int = GAP_PENALTY):
        super().__init__(gap_penalty)
        pass

    def score(self, gap_length: int) -> int:
        return gap_length * self.gap_penalty


class AffineGapPenalty(ScoringFunction):
    """Affine gap penalty scoring function"""

    def __init__(
        self,
        gap_penalty: int = GAP_PENALTY,
        gap_extension_penalty: int = GAP_EXTENSION_PENALTY,
    ):
        super().__init__(gap_penalty)
        self.gap_extension_penalty = gap_extension_penalty

    def score(self, gap_length: int) -> int:
        return self.gap_penalty + (gap_length - 1) * self.gap_extension_penalty


class InvalidSymbolError(Exception):
    """Exception raised when an invalid symbol is found"""

    pass


class SubstitutionMatrix:
    def __init__(self, scores, symbol_to_index: List[str]):
        self.scores = scores
        self.symbol_to_index = symbol_to_index

    def get_score(self, symbol_1: str, symbol_2: str) -> bool:
        if symbol_1 not in self.symbol_to_index:
            raise InvalidSymbolError(
                f"Symbol {symbol_1} is not in the substitution matrix"
            )
        if symbol_2 not in self.symbol_to_index:
            raise InvalidSymbolError(
                f"Symbol {symbol_2} is not in the substitution matrix"
            )
        symbol_1_index = self.symbol_to_index.index(symbol_1)
        symbol_2_index = self.symbol_to_index.index(symbol_2)
        score = self.scores[symbol_1_index][symbol_2_index]
        return score


NUCLEOTIDE_SCORES = [[1, -1, -1, -1], [-1, 1, -1, -1], [-1, -1, 1, -1], [-1, -1, -1, 1]]


class NucleotideSubstitutionMatrix(SubstitutionMatrix):
    """Nucleotide substitution matrix"""

    def __init__(self):
        super().__init__(scores=NUCLEOTIDE_SCORES, symbol_to_index=["A", "C", "G", "T"])
        pass


# BLOSUM62 scores from https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
PROTEIN_SCORES = [
    [4, 1, 2, 2, 0, 1, 1, 0, 2, 1, 1, 1, 1, 2, 1, 1, 0, 3, 2, 0, 2, 1, 0, 4],
    [-1, 5, 0, 2, 3, 1, 0, 2, 0, 3, 2, 2, 1, 3, 2, 1, 1, 3, 2, 3, 1, 0, 1, 4],
    [-2, 0, 6, 1, 3, 0, 0, 0, 1, 3, 3, 0, 2, 3, 2, 1, 0, 4, 2, 3, 3, 0, 1, 4],
    [-2, 2, 1, 6, 3, 0, 2, 1, 1, 3, 4, 1, 3, 3, 1, 0, 1, 4, 3, 3, 4, 1, 1, 4],
    [0, 3, 3, 3, 9, 3, 4, 3, 3, 1, 1, 3, 1, 2, 3, 1, 1, 2, 2, 1, 3, 3, 2, 4],
    [-1, 1, 0, 0, 3, 5, 2, 2, 0, 3, 2, 1, 0, 3, 1, 0, 1, 2, 1, 2, 0, 3, 1, 4],
    [-1, 0, 0, 2, 4, 2, 5, 2, 0, 3, 3, 1, 2, 3, 1, 0, 1, 3, 2, 2, 1, 4, 1, 4],
    [0, 2, 0, 1, 3, 2, 2, 6, 2, 4, 4, 2, 3, 3, 2, 0, 2, 2, 3, 3, 1, 2, 1, 4],
    [-2, 0, 1, 1, 3, 0, 0, 2, 8, 3, 3, 1, 2, 1, 2, 1, 2, 2, 2, 3, 0, 0, 1, 4],
    [-1, 3, 3, 3, 1, 3, 3, 4, 3, 4, 2, 3, 1, 0, 3, 2, 1, 3, 1, 3, 3, 3, 1, 4],
    [-1, 2, 3, 4, 1, 2, 3, 4, 3, 2, 4, 2, 2, 0, 3, 2, 1, 2, 1, 1, 4, 3, 1, 4],
    [-1, 2, 0, 1, 3, 1, 1, 2, 1, 3, 2, 5, 1, 3, 1, 0, 1, 3, 2, 2, 0, 1, 1, 4],
    [-1, 1, 2, 3, 1, 0, 2, 3, 2, 1, 2, 1, 5, 0, 2, 1, 1, 1, 1, 1, 3, 1, 1, 4],
    [-2, 3, 3, 3, 2, 3, 3, 3, 1, 0, 0, 3, 0, 6, 4, 2, 2, 1, 3, 1, 3, 3, 1, 4],
    [-1, 2, 2, 1, 3, 1, 1, 2, 2, 3, 3, 1, 2, 4, 7, 1, 1, 4, 3, 2, 2, 1, 2, 4],
    [1, 1, 1, 0, 1, 0, 0, 0, 1, 2, 2, 0, 1, 2, 1, 4, 1, 3, 2, 2, 0, 0, 0, 4],
    [0, 1, 0, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 1, 1, 5, 2, 2, 0, 1, 1, 0, 4],
    [-3, 3, 4, 4, 2, 2, 3, 2, 2, 3, 2, 3, 1, 1, 4, 3, 2, 11, 2, 3, 4, 3, 2, 4],
    [-2, 2, 2, 3, 2, 1, 2, 3, 2, 1, 1, 2, 1, 3, 3, 2, 2, 2, 7, 1, 3, 2, 1, 4],
    [0, 3, 3, 3, 1, 2, 2, 3, 3, 3, 1, 2, 1, 1, 2, 2, 0, 3, 1, 4, 3, 2, 1, 4],
    [-2, 1, 3, 4, 3, 0, 1, 1, 0, 3, 4, 0, 3, 3, 2, 0, 1, 4, 3, 3, 4, 1, 1, 4],
    [-1, 0, 0, 1, 3, 3, 4, 2, 0, 3, 3, 1, 1, 3, 1, 0, 1, 3, 2, 2, 1, 4, 1, 4],
    [0, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 0, 0, 2, 1, 1, 1, 1, 1, 4],
    [-4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 1],
]

PROTEIN_SYMBOL_TO_INDEX = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
    "B",
    "Z",
    "X",
    "*",
]


class ProteinSubstitutionMatrix(SubstitutionMatrix):
    """Protein substitution matrix"""

    def __init__(self):
        super().__init__(scores=PROTEIN_SCORES, symbol_to_index=PROTEIN_SYMBOL_TO_INDEX)
        pass


class Alignment:
    def __init__(self, sequence_1: str, sequence_2: str):
        self.sequence_1 = sequence_1
        self.sequence_2 = sequence_2

    # TODO: Implement __str__ method
    def __eq__(self, other: object) -> bool:
        return (
            self.sequence_1 == other.sequence_1 and self.sequence_2 == other.sequence_2  # type: ignore
        )

    def __repr__(self) -> str:
        return f"{self.sequence_1}\n{self.sequence_2}"


class TracebackDirection(enum.Enum):
    DIAGONAL = 0
    UPPER = 1
    SIDE = 2


class ScoringMatrix:
    def __init__(
        self,
        sequence_1: str,
        sequence_2: str,
        scoring_function: ScoringFunction,
        substitution_matrix: SubstitutionMatrix,
        match_score: int = 1,
        mismatch_score: int = -1,
    ):
        self.sequence_1 = sequence_1
        self.sequence_2 = sequence_2
        self.scoring_function = scoring_function
        self.substitution_matrix = substitution_matrix
        self.match_score = match_score
        self.mismatch_score = mismatch_score

        self._init_matrices()

    def _init_matrices(self) -> Tuple[List[List[None]], List[List[TracebackDirection]]]:  # type: ignore
        """Initialize 2D matrix for holding scores and traceback directions"""
        logging.info("Prepaing to initialize matrices")
        horizontal_length = len(self.sequence_1) + 1
        vertical_length = len(self.sequence_2) + 1

        # Initialize matrices with None values
        scoring_matrix = [[None] * horizontal_length for _ in range(vertical_length)]
        traceback_matrix = [[None] * horizontal_length for _ in range(vertical_length)]

        gap_penalty = self.scoring_function.gap_penalty
        # Initialize first row with gap penalties times index
        for i in range(horizontal_length):
            scoring_matrix[0][i] = i * gap_penalty  # type: ignore
        # Initialize first column with gap penalties times index
        for j in range(vertical_length):
            scoring_matrix[j][0] = j * gap_penalty  # type: ignore

        self.scoring_matrix = scoring_matrix
        self.traceback_matrix = traceback_matrix
        logging.info(
            "Initialized scoring and traceback matrices both of size %dx%d totalling to %d cells",
            vertical_length,
            horizontal_length,
            vertical_length * horizontal_length,
        )

    def fill(self):
        """Fill 2D matrix with scores"""
        horizontal_length = len(self.sequence_1) + 1
        vertical_length = len(self.sequence_2) + 1

        cells_computed = 0
        total_cells = vertical_length * horizontal_length
        last_percentage = 0
        for j in range(1, vertical_length):
            for i in range(1, horizontal_length):
                # Check symbol equality
                symbol_1 = self.sequence_1[i - 1]
                symbol_2 = self.sequence_2[j - 1]
                symbol_score = self.substitution_matrix.get_score(symbol_1, symbol_2)

                diagonal_value = self.scoring_matrix[j - 1][i - 1]
                upper_value = self.scoring_matrix[j - 1][i]
                side_value = self.scoring_matrix[j][i - 1]

                # Compute score
                diagonal_score = diagonal_value + symbol_score
                upper_score = upper_value + self.scoring_function.score(upper_value)
                side_score = side_value + self.scoring_function.score(side_value)

                # Set score
                scores = [diagonal_score, upper_score, side_score]
                max_score = max(scores)
                self.scoring_matrix[j][i] = max_score

                # Set traceback directions
                traceback_directions = []
                if scores[0] == max_score:
                    traceback_directions.append(TracebackDirection.DIAGONAL)
                if scores[1] == max_score:
                    traceback_directions.append(TracebackDirection.UPPER)
                if scores[2] == max_score:
                    traceback_directions.append(TracebackDirection.SIDE)
                self.traceback_matrix[j][i] = traceback_directions
            cells_computed += horizontal_length

            new_percentage = int(cells_computed / total_cells * 100)
            if new_percentage > last_percentage:
                logging.info("Computed %d%% of cells", new_percentage)
                last_percentage = new_percentage

    def get_optimal_score(self) -> int:
        """Get optimal score from the bottom right corner of the matrix"""
        optimal_score = self.scoring_matrix[-1][-1]
        if optimal_score is None:
            raise ValueError("Matrix is not filled")
        return optimal_score

    def get_alignments(self) -> List[Alignment]:
        """Traceback 2D matrix to find optimal alignments"""
        logging.info("Starting extracting alignments")
        sequence_1_alignment = ""
        sequence_2_alignment = ""
        horizontal_length = len(self.sequence_1)
        vertical_length = len(self.sequence_2)
        j = vertical_length
        i = horizontal_length
        while i > 0 and j > 0:
            traceback_directions = self.traceback_matrix[j][i]
            traceback_direction = traceback_directions[0]  # type: ignore
            # TODO: Handle multiple traceback directions
            if traceback_direction == TracebackDirection.DIAGONAL:
                sequence_1_alignment += self.sequence_1[i - 1]
                sequence_2_alignment += self.sequence_2[j - 1]
                i -= 1
                j -= 1
            elif traceback_direction == TracebackDirection.UPPER:
                sequence_1_alignment += "-"
                sequence_2_alignment += self.sequence_2[j - 1]
                j -= 1
            elif traceback_direction == TracebackDirection.SIDE:
                sequence_1_alignment += self.sequence_1[i - 1]
                sequence_2_alignment += "-"
                i -= 1
            else:
                raise ValueError("Invalid traceback direction")

        # TODO: Reverse alignments
        sequence_1_alignment = sequence_1_alignment[::-1]
        sequence_2_alignment = sequence_2_alignment[::-1]

        return [Alignment(sequence_1_alignment, sequence_2_alignment)]

    def __str__(self) -> str:
        """String representation of the matrix"""
        horizontal_length = (len(self.sequence_1) + 1) * 2
        vertical_length = (len(self.sequence_2) + 1) * 2
        output = [
            ["  " for _ in range(horizontal_length)] for __ in range(vertical_length)
        ]
        for row_idx, row in enumerate(self.scoring_matrix):
            for column_idx, cell in enumerate(row):
                value = str(cell).rjust(2)
                output[row_idx * 2 + 1][column_idx * 2 + 1] = value
        for row_idx, row in enumerate(self.traceback_matrix):
            for column_idx, traceback_directions in enumerate(row):
                if traceback_directions is None:
                    continue
                for traceback_direction in traceback_directions:
                    if traceback_direction == TracebackDirection.DIAGONAL:
                        output[row_idx * 2][column_idx * 2] = " ↖"
                    elif traceback_direction == TracebackDirection.UPPER:
                        output[row_idx * 2][column_idx * 2 + 1] = " ↑"
                    elif traceback_direction == TracebackDirection.SIDE:
                        output[row_idx * 2 + 1][column_idx * 2] = " ←"
        # Stringify output
        output_stringified = "\n".join(["".join(row) for row in output])
        return output_stringified


SCORING_FUNCTIONS = {
    "constant": ConstantGapPenalty,
    "linear": LinearGapPenalty,
    "affine": AffineGapPenalty,
}

SUBSTITUTION_MATRICES = {
    "nucleotide": NucleotideSubstitutionMatrix,
    "protein": ProteinSubstitutionMatrix,
}


class NeedlemanWunsch:
    def __init__(
        self,
        scoring_function: Union[str, ScoringFunction] = "constant",
        substitution_matrix: Union[str, SubstitutionMatrix] = "nucleotide",
    ) -> None:
        # Setup
        if isinstance(scoring_function, str):
            if scoring_function not in SCORING_FUNCTIONS:
                raise ValueError("Invalid scoring function")
            self.scoring_function = SCORING_FUNCTIONS[scoring_function]()
        else:
            self.scoring_function = scoring_function

        if isinstance(substitution_matrix, str):
            if substitution_matrix not in SUBSTITUTION_MATRICES:
                raise ValueError("Invalid substitution matrix")
            self.substitution_matrix = SUBSTITUTION_MATRICES[substitution_matrix]()
        else:
            self.substitution_matrix = substitution_matrix  # type: ignore

    def align(
        self, sequence_1, sequence_2
    ) -> Tuple[List[Alignment], int, ScoringMatrix]:
        """Align two sequences using the Needleman-Wunsch algorithm"""
        scoring_matrix = ScoringMatrix(
            sequence_1, sequence_2, self.scoring_function, self.substitution_matrix
        )
        scoring_matrix.fill()
        alignments = scoring_matrix.get_alignments()
        optimal_score = scoring_matrix.get_optimal_score()
        return alignments, optimal_score, scoring_matrix
