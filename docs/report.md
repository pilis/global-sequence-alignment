# Report

Author: Piotr Pilis <pilispiotr@gmail.com>

This is report for the project 1 of the course Introduction to Bioinformatics at the MiNI faculty of Warsaw University of Technology in the winter semester of 2021/2022.

For more details on how to run the project refer the [README.md](README.md) file.

## Implementation

> Note: The implementation is not optimized for speed. It's optimized for readability and simplicity.

Implementation of the algorithm is done in Python module `global_sequence_alignment`.

The module contain top-level class `NeedlemanWunsch` which wraps the algorithm and provides methods for running it on provided sequences. It's possible to customize the algorithm by providing custom:
- scoring_function: constant (default), linear, affine. In all 3 scoring_functions the gap penalty is set to -1
- substitution_matrix: nucleotide (default) or protein. Substitution Matrix define how well two symbols match. For nucleotide it's a constant value of 1 for match and -1 for mismatch. For protein it's a matrix of values from BLOSUM62 (https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt)

When `align` method is called for 2 sequences it creates `ScoringMatrix` object that first got initialized and then it's filled according to the recursive relation:

- $F_{i-1,j-1} + S_{x[i],y[j]}$ is a score for top-left diagonal cell. It's calculated by taking the value of the top-left cell and adding a score for a match or mismatch of symbols based on a decision from the Substitution matrix (typically "1" or "-1")
- $F_{i-1,j} + G$ is a score for upper cell. It's calculated by taking value from upper cell and adding score from  selected Scoring  function
- $F_{i,j-1} + G$ is a score for side cell calculated same as for upper cell and then we take the maximum value out of those 3

After the matrix is filled we can start backtracking from the bottom-right cell. We take the maximum value from the 3 cells and then we go to the cell that had that value. We repeat this process until we reach the top-left cell. The path we took is the alignment of the sequences.

> Note: The implementation returns one optimal alignment

Example run for made-up sequences:

```bash
python src/main.py --direct GATTACA GTCGACGCA

Optimal score: 0
GAT-TA--CA
G-TCGACGCA
```

Or to save optimal alignment to a file:

```bash
python src/main.py --direct GATTACA GTCGACGCA --output-path=output.txt
```

Next to the implementation I have authored suite of unit tests which can be run with `pytest` command.

On top of Python module I have create simple CLI interface that it's described in the [README.md](README.md) file.

## Homologous genes alignment

To test the implementation I used the sequences of homologous genes from the NCBI database. I used the sequences of mouse and chicken pax6 genes. The output is a table with the alignment of the sequences. DNA similarity can be found here: https://www.ncbi.nlm.nih.gov/homologene/?term=Pax6%5Bgene+name%5D+AND+mouse%5Borgn%5D&report=alignmentscores

> Note: Since mouse.fna and chicken.fna files are fairly large - about 30k symbols each it's not possible to display them in the report. The output is saved in the [output.txt](output.txt) file.

> Time that it took to run the program was too long: I had to wait about 5 minutes to get 1% progress. I think that the reason for that is the size of the sequences. The program is not optimized for speed.

> Scoring and traceback matrices were fairly large: 17577x29516 totalling to 518802732 cells

To run the program:

    python src/main.py ./data/homologous_genes/pax6/mouse.fna ./data/homologous_genes/pax6/chicken.fna --output-path=output.txt

Output:

    2022-11-05 19:31:36,834 - root - INFO - Reading ./data/homologous_genes/pax6/mouse.fna
    2022-11-05 19:31:36,834 - root - INFO - Length: 29515
    2022-11-05 19:31:36,834 - root - INFO - Reading ./data/homologous_genes/pax6/chicken.fna
    2022-11-05 19:31:36,834 - root - INFO - Length: 17576
    2022-11-05 19:31:36,834 - root - INFO - Sequences loaded
    2022-11-05 19:31:36,834 - root - INFO - Prepaing to initialize matrices
    2022-11-05 19:31:42,585 - root - INFO - Initialized scoring and traceback matrices both of size 17577x29516 totalling to 518802732 cells
    2022-11-05 19:33:03,575 - root - INFO - Computed 1% of cells
    ...
    2022-11-05 21:12:48,067 - root - INFO - Computed 99% of cells
    2022-11-05 21:12:54,705 - root - INFO - Starting extracting alignments
    Optimal score: -1413

## Insulin sequences alignment

First I have searched for insuline proteins on NCBI.

Human insuline protein: https://www.ncbi.nlm.nih.gov/gene/3630
Hamster insuline protein: https://www.ncbi.nlm.nih.gov/gene/101823595

Since the sequences are **proteins** I have used `--substitution_matrix=protein` flag to use BLOSUM62 matrix.

To run the program:

    python src/main.py ./data/proteins/insulin/hamster.faa ./data/proteins/insulin/human.faa --substitution_matrix=protein

Result for **constant** scoring function:
```
Optimal score: 520
MTLWMRLLPLLALLVLWEPNPAQAFVNQHLCGSHLVEALYLVCGERGFFYTPKSRRGVEDPQVAQLELGGGPGADDLQTLALEVAQQKRGIVDQCCTSICSLYQLENYCN
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
```

Result for **linear** scoring function:
```
Optimal score: 520
MTLWMRLLPLLALLVLWEPNPAQAFVNQHLCGSHLVEALYLVCGERGFFYTPKSRRGVEDPQVAQLELGGGPGADDLQTLALEVAQQKRGIVDQCCTSICSLYQLENYCN
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
```

Result for **affine** scoring function:
```
Optimal score: 520
MTLWMRLLPLLALLVLWEPNPAQAFVNQHLCGSHLVEALYLVCGERGFFYTPKSRRGVEDPQVAQLELGGGPGADDLQTLALEVAQQKRGIVDQCCTSICSLYQLENYCN
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
```

Since the sequences are identical the optimal score is the same for all 3 scoring functions
