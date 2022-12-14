import logging
import sys

import click

from global_sequence_alignment.needleman_wunsch import NeedlemanWunsch

root = logging.getLogger()
root.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
root.addHandler(handler)


def read_fasta_file(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
        return "".join([line.strip() for line in lines[1:]])


def write_optimal_alignments_to_file(file_path, alignments):
    with open(file_path, "w") as f:
        for alignment in alignments:
            alignment_serialized = str(alignment)
            f.write(alignment_serialized)


@click.command()
@click.argument("sequence_1")
@click.argument("sequence_2")
@click.option("--scoring_function", default="constant")
@click.option("--substitution_matrix", default="nucleotide")
@click.option(
    "--direct",
    is_flag=True,
    help="If set, instead of reading sequences from files, the sequences are read from the command line",
)
@click.option(
    "--print_scoring_matrix",
    is_flag=True,
    help="If set, the scoring matrix is printed to the console",
)
@click.option("--output-path")
def main(
    sequence_1: str,
    sequence_2: str,
    scoring_function: str,
    substitution_matrix: str,
    direct: bool,
    output_path: str,
    print_scoring_matrix: bool = False,
) -> None:
    """Run Needleman-Wunsch algorithm"""
    if not direct:
        logging.info(f"Reading {sequence_1}")
        sequence_1 = read_fasta_file(sequence_1)
        logging.info(f"Length: {len(sequence_1)}")
        logging.info(f"Reading {sequence_2}")
        sequence_2 = read_fasta_file(sequence_2)
        logging.info(f"Length: {len(sequence_2)}")
    logging.info("Sequences loaded")

    # Execute Needleman-Wunsch algorithm
    needleman_wunsch = NeedlemanWunsch(scoring_function, substitution_matrix)
    alignments, optimal_score, scoring_matrix = needleman_wunsch.align(
        sequence_1, sequence_2
    )

    # Print optimal score
    print("\n")
    print(f"Optimal score: {optimal_score}")

    # Print optimal alignments
    if output_path:
        write_optimal_alignments_to_file(output_path, alignments)
    else:
        print("\n")
        print("Optimal alignments:")
        for alignment_idx, alignment in enumerate(alignments):
            print("Found alignment #{}:".format(alignment_idx + 1))
            print(alignment)
            print("\n")

    # Print scoring matrix
    print("\n")
    print("Scoring and traceback matrix:")
    if print_scoring_matrix:
        print(scoring_matrix)


if __name__ == "__main__":
    main()
