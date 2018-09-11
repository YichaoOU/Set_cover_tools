# The set cover algorithms for motif selection.

### What is motif selection?

The motif selection problem seeks to identify a minimal set of regulatory motifs that characterize sequences of interest (e.g. ChIP-seq binding regions).

### Dependencies

- Pandas

- Emotif Alpha: https://github.com/YichaoOU/Emotif_Alpha

### The set cover algorithms

1. Greedy approach

This algorithm is developed by Rami Al-Ouran. I created a wrapper for his code inside the Greedy folder.

Al-Ouran, Rami, et al. "Discovering gene regulatory elements using coverage-based heuristics." IEEE/ACM transactions on computational biology and bioinformatics 15.4 (2018): 1290-1300.

2. Tabu search

This algorithm is developed by Yating Liu. 

To compile the code, please first download the METSlib framework. Then:

`g++ -I/PATH/metslib-0.5 main.cc`

To run the code: `tabu_MS -in input_matrix -out output_file_name -t 0.2`

3. Relaxed integer linear programming

This algorithm is developed by David Juedes. 

To compile the code, please first download the GLPK framework. Then:

`g++ -std=c++11 RIPL_simplex2_07_31.cc -o msdc -L/PATH/glpk/lib -I/PATH/glpk/include -lglpk`

I created a wrapper for his code: `python MSDC_wrapper.py input_matrix output_file_name pos_k neg_j`

### The nested cross-validation framework

The code is in the Evams (short for Evaluation of Motif Selection) folder
