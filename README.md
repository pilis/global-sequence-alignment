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

    python src/main.py GATTACA GTCGACGCA
