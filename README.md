# Global sequence alignment

This is a simple implementation of the Needleman-Wunsch algorithm for global sequence alignment

## Getting started

Prerequisities:
- Python 3.9.7+
- pre-commit

To create virtual environment:

    python -m venv venv
    source venv/bin/activate

To install Python dependencies:

    pip install -r requirements-dev.txt -r src/requirements.txt

To install pre-commit hooks:

    pre-commit install

To run tests:

    pytest

To run tests with coverage:

    coverage run --branch -m pytest && coverage report --omit="*/tests*"

## Usage

To run the program for nucleotide sequences:

    python src/main.py ./data/homologous_genes/pax6/mouse.fna ./data/homologous_genes/pax6/chicken.fna --output-path=output.txt

To run the program for protein sequences:

    python src/main.py ./data/proteins/insulin/hamster.faa ./data/proteins/insulin/human.faa --substitution_matrix=protein

To run the program with directly provided sequences:

    python src/main.py --direct GATTACA GTCGACGCA

To run the program with directly provided sequences and save the output to a file:

    python src/main.py --direct GATTACA GTCGACGCA --output-path=output.txt
