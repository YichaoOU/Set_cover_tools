Relaxed Integer Linear Programming
==================================


Usage
-----


**Step 0. Download the GLPK library**

https://www.gnu.org/software/glpk/

**Step 1. Compile the C code**

.. code:: bash

	g++ -std=c++11 RIPL_simplex2_07_31.cc -o msdc -L/PATH/glpk/lib -I/PATH/glpk/include -lglpk

**Step 2. Run the RILP method**

.. code:: bash

	./msdc foreground.list background.list motif_mapping.fimo

You can also use the python wrapper:

.. code:: bash

	python MSDC_wrapper.py [input_matrix] [output_file_name] [pos_k] [neg_j]

Below is an example:

.. code:: bash

	python MSDC_wrapper.py AP1.Motif_Selection_Paper.csv AP1_2_2 2 2


Parameters
^^^^^^^^^^

``pos_k``: the maximal uncovered percent of foreground set. Possible values are integers from 1 to 10, where 1 means 10%.

``neg_j``: the maximal covered percent of background set. Possible values are integers from 1 to 10, where 1 means 10%.



















