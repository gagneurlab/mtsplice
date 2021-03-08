Multi-tissue Splicing Model (MTSplice)
==========================================

This repository contains code to reproduce the analysis in the MTSplice paper. 
To run the code, you need to first install the dependencies, provided in `environment.yml`. We recommend using conda environment.

```
conda env create -f environment.yml
conda install cyvcf2 cython -y
```

You would also need to install this repository as a python package. At the root directory of this repo:

```
pip install mmsplice==1.0.1
pip install -e .
```

In addition, we use R for plotting. The libraries required are:

```
cowplot
data.table
ggplot2
magrittr
rhdf5
ggrepel
ggpval
```

You also need to download gencode 27 reference sequence being placed in `data/` named as `hg19.fa` and the index `hg19.fa.fai`.

Structure of the repo
---------------------

- `mtsplice/*` python class and functions used in the analysis  
- `notbooks/*:` Jupyter notebooks to run predictions with the MTSplice model
- `config.R` Configuration setting for the plots
- `Figure1-7_*` R script to reproduce plots in the manuscript
- `data/*` processed data