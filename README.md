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

To run the program:

    python src/main.py ./data/homologous/pax6/mouse.fna ./data/homologous/pax6/chicken.fna

> This is comparison of mouse and chicken pax6 genes. The output is a table with the alignment of the sequences. DNA similarity can be found here: https://www.ncbi.nlm.nih.gov/homologene/?term=Pax6%5Bgene+name%5D+AND+mouse%5Borgn%5D&report=alignmentscores

To run the program for proteins:

    python src/main.py ./data/proteins/insulin/hamster.faa ./data/proteins/insulin/human.faa --substitution_matrix=protein

To run the program with directly provided sequences:

    python src/main.py --direct GATTACA GTCGACGCA
