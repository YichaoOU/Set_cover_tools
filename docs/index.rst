The set cover algorithms for motif selection
============================================

.. toctree::
   :maxdepth: 3
   
   RILP
   Tabu

Motivation
----------

De novo motif discovery algorithms find statistically over-represented sequence motifs that may function as transcription factor binding sites. Current methods often report large numbers of motifs, making it difficult to perform further analyses and validation. The motif selection problem seeks to identify a minimal set of regulatory motifs that characterize sequences of interest (e.g. ChIP-Seq binding regions).

What is motif selection?
----------

The motif selection problem seeks to identify a minimal set of regulatory motifs that characterize sequences of interest (e.g. ChIP-seq binding regions). The output motifs represent putative binding sites for primary transcription factors (ChIP-ed factors) and co-factors.


Our solution to motif selection
----------

The greedy set cover solution
""""""""""

This is our first generation solution. For more details, see: https://github.com/RamiOran/SeqCov

Contrast to the first generation solution, our current methods introduce a background set, which ideally, we wish to (1) cover as much of the foreground as possible, (2) cover as little of the background as possible, and (3) select the smallest set of motifs. See below:

The tabu search method
"""""""""

Tabu search is a metaheuristic local search method. It starts with a randomly generated initial solution then searches its neighborhood. Traditional local search methods, such as hill climbing, update current solution if they find a better solution in the neighborhood and thus result in local optima. In contrast, tabu search alleviates this issue by employing two strategies: (1) tabu search accepts non-improving moves when better moves are unavailable in the neighborhood of current solution and (2) tabu search uses a short-term memory structure, called tabu list, to store recently visited solutions and prevent selecting solutions that are visited previously. 

The relaxed integer linear programming method
"""""""""""""""""""""""

This algorithm contains two steps. The first step is to model the motif selection problem as a linear programming problem and obtain an optimal solution via GLPK. The second step is to solve the integer-version problem through a randomized algorithm. 



