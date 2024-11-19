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

MTEApy is a Python library for Metabolic Task Enrichment Analysis (MTEA) that leverages the use of powerful contraint-based metabolic model frameworks. It uses metabolic tasks to inferr the metabolic states or changes using transcriptomic data.

## Installation

> [!WARNING] 
> Not implemented yet!

To install MTEApy, you can install it using `pip`:

```sh
pip install mteapy
```

Or you can download this repository and install it locally, again using `pip`:

```sh
git clone https://github.com/bsc-life/mteapy/
pip install -e mteapy/
```

## MTEA Usage

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

### Secretory tasks

MTEApy allows the use of secretory tasks under the secretory module accessible using the option `-s` or `--secretory`, which will run the command using both the metabolic and secretory modules [[4](#references)]. This option is only available for the CellFie and TIDE frameworks, and leverages the use of an updated version of the metabolic model HumanGEM with secretory reactions and an updated task structure.

```sh
run-mtea {TIDE|Cellfie} --secretory [...] expr_file
```

More info on the curation of secretory tasks and the updated metabolic model can be found at the repository [xavibemo/secretory-tasks](https://github.com/xavibemo/secretory-tasks).

## References

1. Richelle, A.; Kellman, B.P.; Wenzel, A.T.; Chiang, A.W.; Reagan, T.; Gutierrez, J.M.; Joshi, C.; Li, S.; Liu, J.K.; Masson, H.; _et al._ Model-based assessment of mammalian cell metabolic functionalities using omics data. _Cell Reports Methods_ **2021**, 1, 100040. https://doi.org/10.1016/j.crmeth.2021.100040.
2. Dougherty, B.V.; Rawls, K.D.; Kolling, G.L.; Vinnakota, K.C.; Wallqvist, A.; Papin, J.A. Identifying functional metabolic shifts in heart failure with the integration of omics data and a heart-specific, genome-scale model. _Cell Reports_ **2021**, 34, 108836. https://doi.org/10.1016/j.celrep.2021.108836.
3. Robinson, J.L.; Kocabaş, P.; Wang, H.; Cholley, P.E.; Cook, D.; Nilsson, A.; Anton, M.; Ferreira, R.; Domenzain, I.; Billa, V.; _et al_. An atlas of human metabolism. _Science Signaling_ **2020**, 13, eaaz1482. https://doi.org/10.1126/scisignal.aaz1482.
4. Masson, H. O., Samoudi, M., Robinson, C. M., Kuo, C.-C., Weiss, L., Shams Ud Doha, K., Campos, A., Tejwani, V., Dahodwala, H., Menard, P., Voldborg, B. G., Robasky, B., Sharfstein, S. T., & Lewis, N. E. Inferring secretory and metabolic pathway activity from omic data with secCellFie. _Metabolic Engineering_, **2024**, 81, 273–285. https://doi.org/10.1016/j.ymben.2023.12.006