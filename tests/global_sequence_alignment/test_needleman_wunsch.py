from unittest import TestCase

from global_sequence_alignment.needleman_wunsch import (
    InvalidSymbolError,
    NucleotideSubstitutionMatrix,
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


# class TestNeedlemanWunsch(TestCase):
#     def test_equal_sequences(self):
#         sequence_1 = "ACGT"
#         sequence_2 = "ACGT"

#         needleman_wunsch = NeedlemanWunsch()
#         alignments = needleman_wunsch.compare(sequence_1, sequence_2)

#         self.assertEqual(len(alignments), 1)
