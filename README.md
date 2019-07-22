# The set cover algorithms for motif selection.

## What is motif discovery?

A common biological question when dealing with ChIP-seq data is to identify transcription factor binding sites. This problem is called motif discovery. Motif discovery algorithms produce a large number of motifs (e.g., MEME) and they don’t evaluate motifs systematically. Therefore, several methods (e.g., PMID: 24335146) have been developed to select high-confidence motifs. 

## What is motif selection?

The motif selection problem seeks to identify a minimal set of regulatory motifs that characterize sequences of interest (e.g. ChIP-seq binding regions). The output motifs represent putative binding sites for primary transcription factors (ChIP-ed factors) and co-factors.

## Getting started

Please see the [documentation](https://set-cover-tools.readthedocs.io/en/latest/index.html) for more details.

The following example has been tested in Ubuntu x64.

**Step 0. Download the GLPK library**

https://www.gnu.org/software/glpk/

**Step 1. Compile the C code**

.. code:: bash

	cd RILP

	g++ -std=c++11 RIPL_simplex2_07_31.cc -o msdc -L/PATH/glpk/lib -I/PATH/glpk/include -lglpk

**Step 2. Run the RILP method**

.. code:: bash

	./msdc foreground.list background.list motif_mapping.fimo

**Output**

The output is a list of selected motifs (i.e., subset of the input motifs) defined in `motif_mapping.fimo`.

## Data folder

This folder contains all raw inputs (i.e., binary matrices) and the Fisher exact test p-values for every motif.

## To replicate our results

Please use the detailed installation steps below

### **Install Anaconda Python 2.7 version**

https://www.anaconda.com/distribution/

### **Create Conda environment**

`conda create -n set_cover python=2.7`

`conda activate set_cover` or `source activate set_cover`

### **Install dependencies**

* **Pandas:** `conda install -c anaconda pandas`
* **Joblib:** `conda install -c anaconda joblib`
* **scikit-learn:** `conda install -c anaconda scikit-learn`

### **Run the nested cross-validation framework**

`cd Evams`

Note that Evams imports `conf.py` from the same directory. You can add `Evams` to your env path by doing `export PATH=$PATH:/path_to_Evams_folder/`

`Evams2 -h`

`Evams2 -jid test -confFile example.conf`

In `example.conf`, please change input to individual files in the Data folder.

Configuration files used in this paper are provided in the **conf_paper** folder.

`-jid` is job ID. An output folder will be created using job ID, in which all results will be put.












