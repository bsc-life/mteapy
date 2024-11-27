# MTEApy

     __  __   ___                      __.....__                 
    |  |/  `.'   `.               .-''         '.               
    |   .-.  .-.   '      .|     /     .-''"'-.  `.             
    |  |  |  |  |  |    .' |_   /     /________\   \     __     
    |  |  |  |  |  |  .'     | |                  |  .:--.'.   
    |  |  |  |  |  | '--.  .-'  \    .-------------' / |   \ |  
    |  |  |  |  |  |    |  |     \    '-.____...---. `" __ | |  
    |__|  |__|  |__|    |  |      `.             .'   .'.''| |  
                        |  '.'      `''-...... -'    / /   | |_ 
                        |   /                        \ \._,\ '/ 
                        `'-'                          `--'  `"  

MTEApy is a Python library for **Metabolic Task Enrichment Analysis** (MTEA) that leverages the use of powerful contraint-based metabolic model frameworks. It uses metabolic tasks to inferr the metabolic states or changes using transcriptomic data.

## Installation

To install MTEApy, you can install it using `pip`:

```sh
pip install mteapy
```

Or you can download this repository and install it locally, again using `pip`:

```sh
git clone https://github.com/bsc-life/mteapy/
pip install -e mteapy/
```

## Overview

MTEApy is comprised of two main contraint-based metabolic modeling frameworks, TIDE and CellFie, implemented in Python (the original source codes are published in Matlab at their respective repositories). Each framework runs using different types of input files.

| Framework | Original Code | Usage |
| --------- | ------------- | ----- |
| **CellFie** [[1](#references)] | [LewisLabUCSD/CellFie](https://github.com/LewisLabUCSD/CellFie) |  To run the CellFie framework you need a gene expression matrix (where rows are genes, columns samples and each cell contains an expression value), a metabolic model, and a metabolic task structure. |
| **TIDE** [[2](#references)] | [csbl/iCardio](https://github.com/csbl/iCardio) | To run the TIDE framework you need a Differential Expression Analysis result file (which should contain a gene name/ID column, a log-FC values column, and a significance column), a metabolic model, and a metabolic task structure. |
| **TIDE-essential** | - | To run the TIDE-essential framework you need the same requirements as TIDE's. | 

MTEApy is designed to be used both as a command-line tool and as a Python module in a Jupyter Notebook or Python script.

### Command-line

If used as a command-line tool, run the command `run-mtea` and specify the desired framework. By default, the metabolic model used by the command is the Human-GEM [[3](#references)] and, therefore, the metabolic tasks are also compatible with Human-GEM.

```sh
run-mtea [-h] [-v] [-c] [-t] [-s] {TIDE-essential,TIDE,CellFie}
```
For more details on the input parameters, run the `-h` or `--help` after any of the commands.

### Python module

If used as a Python module, import the `mteapy` module or directly import the desired functions to compute a framework.

```python
from mteapy.tide import compute_TIDE, compute_TIDEe
from mteapy.cellfie import compute_CellFie
```

## References

1. Richelle, A.; Kellman, B.P.; Wenzel, A.T.; Chiang, A.W.; Reagan, T.; Gutierrez, J.M.; Joshi, C.; Li, S.; Liu, J.K.; Masson, H.; _et al._ Model-based assessment of mammalian cell metabolic functionalities using omics data. _Cell Reports Methods_ **2021**, 1, 100040. https://doi.org/10.1016/j.crmeth.2021.100040.
2. Dougherty, B.V.; Rawls, K.D.; Kolling, G.L.; Vinnakota, K.C.; Wallqvist, A.; Papin, J.A. Identifying functional metabolic shifts in heart failure with the integration of omics data and a heart-specific, genome-scale model. _Cell Reports_ **2021**, 34, 108836. https://doi.org/10.1016/j.celrep.2021.108836.
3. Robinson, J.L.; Kocabaş, P.; Wang, H.; Cholley, P.E.; Cook, D.; Nilsson, A.; Anton, M.; Ferreira, R.; Domenzain, I.; Billa, V.; _et al_. An atlas of human metabolism. _Science Signaling_ **2020**, 13, eaaz1482. https://doi.org/10.1126/scisignal.aaz1482.

***
# Support
Xavier Benedicto:        [xavier.benedicto@bsc.es](mailto:xavier.benedicto@bsc.es)
Miguel Ponce-de-León:    [miguel.ponce@bsc.es](mailto:miguel.ponce@bsc.es)