from unittest import TestCase
from unittest.mock import MagicMock

from global_sequence_alignment.needleman_wunsch import (
    InvalidSymbolError,
    NucleotideSubstitutionMatrix,
    ScoringMatrix,
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


# class TestNeedlemanWunsch(TestCase):
#     def test_equal_sequences(self):
#         sequence_1 = "ACGT"
#         sequence_2 = "ACGT"

#         needleman_wunsch = NeedlemanWunsch()
#         alignments = needleman_wunsch.compare(sequence_1, sequence_2)

#         self.assertEqual(len(alignments), 1)
