---
title: Home
layout: home
nav_order: 1
---

# **MTEApy**

***
**MTEApy** is a Python library for **Metabolic Task Enrichment Analysis** (MTEA) that leverages the use of powerful contraint-based metabolic model frameworks. It uses metabolic tasks to inferr the metabolic states or changes using transcriptomic data. This bundle of frameworks has been created to facilitate the access to contraint-based metabolic modeling approaches to researchers without the need of advanced bioinformatic skills in a simple Python library.


## Installation

To install MTEApy, you can do so using `pip`:

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

| Framework | Original Code | Description |
|:--------- |:------------- |:----------- |
| **CellFie** [[1](#references)] | [LewisLabUCSD/CellFie](https://github.com/LewisLabUCSD/CellFie) | Utilises a normalized expression matrix (e.g., TPMs) to compute a gene activity score using user-defined thresholds, and then projects it into metabolic reactions. Using the participating reactions for each metabolic task, a metabolic score is computed which indicates the metabolic activity of the metabolic tasks across samples. |
| **TIDE** [[2](#references)] | [csbl/iCardio](https://github.com/csbl/iCardio) | Utilises a differential expression result and its log-FC values to project them into metabolic reactions. Using the participating reactions for each metabolic task, a metabolic score is computed which indicates the change in metabolic activity for one control-sample. A p-value is assigned to each score after performing a permutation test. |
| **TIDE-essential** | [bsc-life/mteapy](https://github.com/bsc-life/mteapy) | Utilises a differential expression result, its log-FC and essential genes to metabolic tasks to compute a metabolic score which indicates the change in metabolic activity for one control-sample. A p-value is assigned to each score after performing a permutation test. | 

MTEApy is designed to be used both as a command-line tool and as a Python module in a Jupyter Notebook or Python script. By default, the metabolic model used by the command is the Human-GEM [[3](#references)].

{: .note}
Only the **Human-GEM** and its metabolic tasks are implemented, so the framework will only take in **EnsemblIDs** as valid genic nomenclature. We are working to allow for any metabolic model and metabolic tasks to be used for more customisable analyses!

### Command-line

If used as a command-line tool, run the command `run-mtea` and specify the desired framework. 

```sh
run-mtea [-h] [-v] [-c] [-t] [-s] {TIDE-essential,TIDE,CellFie}
```
For more details on the input parameters, run the `-h` or `--help` after any of the commands or see the dedicated pages.

### Python module

If used as a Python module, import the `mteapy` module or directly import the desired wrapper functions to compute a framework.

```python
from mteapy.tide import compute_TIDE, compute_TIDEe
from mteapy.cellfie import compute_CellFie
```

## Citation

> Comming soon!


## References

1. Richelle, A.; Kellman, B.P.; Wenzel, A.T.; Chiang, A.W.; Reagan, T.; Gutierrez, J.M.; Joshi, C.; Li, S.; Liu, J.K.; Masson, H.; _et al._ Model-based assessment of mammalian cell metabolic functionalities using omics data. _Cell Reports Methods_ **2021**, 1, 100040. [https://doi.org/10.1016/j.crmeth.2021.100040](https://doi.org/10.1016/j.crmeth.2021.100040).
2. Dougherty, B.V.; Rawls, K.D.; Kolling, G.L.; Vinnakota, K.C.; Wallqvist, A.; Papin, J.A. Identifying functional metabolic shifts in heart failure with the integration of omics data and a heart-specific, genome-scale model. _Cell Reports_ **2021**, 34, 108836. [https://doi.org/10.1016/j.celrep.2021.108836](https://doi.org/10.1016/j.celrep.2021.108836).
3. Robinson, J.L.; Kocabaş, P.; Wang, H.; Cholley, P.E.; Cook, D.; Nilsson, A.; Anton, M.; Ferreira, R.; Domenzain, I.; Billa, V.; _et al_. An atlas of human metabolism. _Science Signaling_ **2020**, 13, eaaz1482. [https://doi.org/10.1126/scisignal.aaz1482](https://doi.org/10.1126/scisignal.aaz1482).


***
## Contact

- Xavier Benedicto Molina ([xavier.benedicto@bsc.es](mailto:xavier.benedicto@bsc.es))
- Miguel Ponce-de-León ([miguel.ponce@bsc.es](mailto:miguel.ponce@bsc.es))