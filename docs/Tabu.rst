Tabu search for motif selection
===============================


Usage
-----


**Step 0. Download the GLPK library**

https://projects.coin-or.org/metslib

**Step 1. Compile the C code**

.. code:: bash

	g++ -I/PATH/metslib-0.5 main.cc

**Step 2. Run the RILP method**

.. code:: bash

	tabu_MS -in input_matrix -out output_file_name -t 0.2

Parameters
----------


`-in` inputFile
`-out` outputFile
`-a` alpha (between [0, 1])
`-b` beta (default is the number of motifs)
`-t` tenure (between (0, 1))
`-max` maxIteration (default 1000)
`-d` delta (between [0,1], default 0.05)
`-p` penalty (default 10000)
`-iter` iterations (default 1, the times of running the whole program. )


















