import click

from global_sequence_alignment.needleman_wunsch import NeedlemanWunsch


@click.command()
@click.argument("sequence_1")
@click.argument("sequence_2")
@click.option("--scoring_function", default="constant")
@click.option("--substitution_matrix", default="nucleotide")
def main(
    sequence_1: str, sequence_2: str, scoring_function: str, substitution_matrix: str
) -> None:
    """Run Needleman-Wunsch algorithm"""
    needleman_wunsch = NeedlemanWunsch(scoring_function, substitution_matrix)
    alignments, optimal_score = needleman_wunsch.align(sequence_1, sequence_2)
    print(f"Optimal score: {optimal_score}")
    for alignment in alignments:
        print(alignment)


if __name__ == "__main__":
    main()
