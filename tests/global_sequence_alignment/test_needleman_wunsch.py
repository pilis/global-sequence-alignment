from unittest import TestCase
from unittest.mock import MagicMock

from global_sequence_alignment.needleman_wunsch import (
    Alignment,
    ConstantGapPenalty,
    InvalidSymbolError,
    NucleotideSubstitutionMatrix,
    ScoringMatrix,
    TracebackDirection,
)


class TestNucleotideSubstitutionMatrix(TestCase):
    def test_same_valid_symbols(self):
        symbol_1 = "A"
        symbol_2 = "A"

        substitution_matrix = NucleotideSubstitutionMatrix()
        is_equal = substitution_matrix.is_equal(symbol_1, symbol_2)

        self.assertTrue(is_equal)

    def test_different_valid_symbols(self):
        symbol_1 = "A"
        symbol_2 = "C"

        substitution_matrix = NucleotideSubstitutionMatrix()
        is_equal = substitution_matrix.is_equal(symbol_1, symbol_2)

        self.assertFalse(is_equal)

    def test_invalid_symbol(self):
        symbol_1 = "X"
        symbol_2 = "Y"

        substitution_matrix = NucleotideSubstitutionMatrix()
        with self.assertRaises(InvalidSymbolError):
            substitution_matrix.is_equal(symbol_1, symbol_2)


class TestScoringMatrix(TestCase):
    def test_building_scoring_matrix_for_empty_sequences(self):
        sequence_1 = ""
        sequence_2 = ""
        scoring_function = MagicMock(gap_penalty=-1)
        substitution_matrix = MagicMock()

        scoring_matrix = ScoringMatrix(
            sequence_1, sequence_2, scoring_function, substitution_matrix
        )

        self.assertEqual(scoring_matrix.scoring_matrix, [[0]])
        self.assertEqual(scoring_matrix.traceback_matrix, [[None]])

    def test_building_scoring_matrix_for_nonempty_sequences(self):
        sequence_1 = "ABC"
        sequence_2 = "ABC"
        scoring_function = MagicMock(gap_penalty=-1)
        substitution_matrix = MagicMock()

        scoring_matrix = ScoringMatrix(
            sequence_1, sequence_2, scoring_function, substitution_matrix
        )

        expected_scoring_matrix = [
            [0, -1, -2, -3],
            [-1, None, None, None],
            [-2, None, None, None],
            [-3, None, None, None],
        ]
        self.assertEqual(scoring_matrix.scoring_matrix, expected_scoring_matrix)

    def test_filling_for_single_symbol_sequences(self):
        sequence_1 = "A"
        sequence_2 = "A"
        scoring_function = MagicMock(gap_penalty=-1)
        scoring_function.score.return_value = -1
        substitution_matrix = MagicMock()
        substitution_matrix.is_equal.return_value = True

        scoring_matrix = ScoringMatrix(
            sequence_1, sequence_2, scoring_function, substitution_matrix
        )

        scoring_matrix.fill()

        expected_scoring_matrix = [[0, -1], [-1, 1]]
        self.assertEqual(scoring_matrix.scoring_matrix, expected_scoring_matrix)

        expected_traceback_matrix = [
            [None, None],
            [None, [TracebackDirection.DIAGONAL]],
        ]
        self.assertEqual(scoring_matrix.traceback_matrix, expected_traceback_matrix)

    def test_filling_for_multiple_symbol_sequences(self):
        sequence_1 = "ATC"
        sequence_2 = "ATC"
        scoring_function = ConstantGapPenalty(gap_penalty=-1)
        substitution_matrix = NucleotideSubstitutionMatrix()

        scoring_matrix = ScoringMatrix(
            sequence_1, sequence_2, scoring_function, substitution_matrix
        )

        scoring_matrix.fill()

        expected_scoring_matrix = [
            [0, -1, -2, -3],
            [-1, 1, 0, -1],
            [-2, 0, 2, 1],
            [-3, -1, 1, 3],
        ]
        self.assertEqual(scoring_matrix.scoring_matrix, expected_scoring_matrix)

        expected_traceback_matrix = [
            [None, None, None, None],
            [
                None,
                [TracebackDirection.DIAGONAL],
                [TracebackDirection.SIDE],
                [TracebackDirection.SIDE],
            ],
            [
                None,
                [TracebackDirection.UPPER],
                [TracebackDirection.DIAGONAL],
                [TracebackDirection.SIDE],
            ],
            [
                None,
                [TracebackDirection.UPPER],
                [TracebackDirection.UPPER],
                [TracebackDirection.DIAGONAL],
            ],
        ]
        self.assertEqual(scoring_matrix.traceback_matrix, expected_traceback_matrix)

    def test_getting_alignment_for_empty_sequences(self):
        sequence_1 = ""
        sequence_2 = ""
        scoring_function = MagicMock(gap_penalty=-1)
        substitution_matrix = MagicMock()

        scoring_matrix = ScoringMatrix(
            sequence_1, sequence_2, scoring_function, substitution_matrix
        )

        alignments = scoring_matrix.get_alignments()

        expected_alignment = Alignment("", "")
        self.assertEqual(alignments, [expected_alignment])

    def test_getting_alignment_for_multiple_symbol_sequences(self):
        sequence_1 = "ATC"
        sequence_2 = "ATC"
        scoring_function = ConstantGapPenalty(gap_penalty=-1)
        substitution_matrix = NucleotideSubstitutionMatrix()

        scoring_matrix = ScoringMatrix(
            sequence_1, sequence_2, scoring_function, substitution_matrix
        )

        scoring_matrix.fill()

        alignments = scoring_matrix.get_alignments()

        expected_alignment = Alignment("ATC", "ATC")
        self.assertEqual(alignments, [expected_alignment])

    def test_getting_alignment_for_multiple_symbol_sequences_with_differences(self):
        sequence_1 = "GA"
        sequence_2 = "G"
        scoring_function = ConstantGapPenalty(gap_penalty=-1)
        substitution_matrix = NucleotideSubstitutionMatrix()

        scoring_matrix = ScoringMatrix(
            sequence_1, sequence_2, scoring_function, substitution_matrix
        )

        scoring_matrix.fill()

        alignments = scoring_matrix.get_alignments()

        expected_alignment = Alignment("GA", "G-")
        self.assertEqual(alignments, [expected_alignment])


# class TestNeedlemanWunsch(TestCase):
#     def test_equal_sequences(self):
#         sequence_1 = "ACGT"
#         sequence_2 = "ACGT"

#         needleman_wunsch = NeedlemanWunsch()
#         alignments = needleman_wunsch.compare(sequence_1, sequence_2)

#         self.assertEqual(len(alignments), 1)
