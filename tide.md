---
title: TIDE
layout: default
nav_order: 2
---

# **The TIDE framework**
{: .no_toc }
***

## Table of contents
{: .no_toc .text-delta }
1. TOC
{:toc}
***

## Description

The **T**ask **I**nferred from **D**ifferential **E**xpression (TIDE) framework is a contraint-based metabolic modeling framework that was originally published by [Dougherty _et al._, 2021](https://doi.org/10.1016/j.celrep.2021.108836). It leverages the use of mathematical descriptions of metabolic functions (metabolic tasks) and the results of a Differential Expression Analysis to study metabolic perturbations in a case control assay.

## Command options

| Argument | Shortcut | Description | Default |
|:-------- |:-------- |:----------- |:------- |
| `dea_file` | | Filename for a differential expression analysis results file. It should contain at least three columns: genic (string), log-FC (numeric) and significance (numeric, e.g.: p-value, adjusted p-value, FDR). | |
| `--delim` | `-d` | Field delimiter for inputed file. | `\t` |
| `-out` | `-o` | Name (and location) to store the analysis' results. They will be stored in a tab-sepparated file, so filenames should contain the  `.tsv` or `.txt` extensions. | `tide_results.tsv` |
| `--gene_col` | | Name of the column in the inputed file containing gene names/symbols. | `geneID` |
| `--lfc_col` | | Name of the column in the inputed file containing log-FC values. | `log2FoldChange` |
| `--pvalue_col` | | Name of the column in the inputed file containing significance values. Only required if the flag `--mask_lfc_values` is `True`. | `padj` |
| `--alpha` | `-a` | Significance threshold to mask log-FC. Only required if the flag `--mask_lfc_values` is `True`. | `0.05` |
| `--n_permutations` | `-n` | Number of permutations to infer p-values for the metabolic scores. The resolution of the computed p-values will depend on this number. | `1000` |
| `--n_cpus` | | Number of CPUs for parallel execution. | `1` |
| `--or_func` | | Name of the function that will be used to resolve OR relationships in gene-protein-reaction (GPR) rules. Possible values are `absmax`, which will return the absolute maximum value, and `max`, which will return the maximum value. | `absmax` |
| `--mask_lfc_values` | | Flag to indicate whether to mask log-FC values to 0 according to their significance. That is, if a log-FC value is non-significant (determined by the user), they will be masked to 0. | `False` |
| `--random_scores` | | Flag to indicate whether to return the null distribution of random scores used to inferr significance with the results file. | `False` |

## Usage example

### Differential expression analysis

The first thing that the TIDE framework requires is a Differential Expression Analysis (DEA) result. Usually, this kind of data is stored in a tabular format and contains at least three columns: gene names/symbols, expression change values (log-FC) and significancy (p-value). 

A typical DEA result will look like the following:

```
                geneID geneSymbol  log2FoldChange      padj
0      ENSG00000000003     TSPAN6        3.710229  0.259406
1      ENSG00000000005       TNMD       -2.437056  0.485180
2      ENSG00000000419       DPM1        8.749658  0.802934
3      ENSG00000000457      SCYL3      -10.409959  0.051220
4      ENSG00000000460      FIRRM       -0.977916  0.926198
...                ...        ...             ...       ...
```

### Running TIDE

To run the TIDE framework using the command-line, the command `run-mtea TIDE` should be used with the desired arguments. A typical TIDE analysis is run using a range of `1,000` to `10,000` permutations, the `absmax` function to evaluate OR GPR rules, and selecting the `--mask_lfc_values` flag, which will mask non-significant log-FC values to 0.

{: .note}
Only the **Human-GEM** and its metabolic tasks are implemented, so the framework will only take in **EnsemblIDs** as valid genic nomenclature. We are working to allow for any metabolic model and metabolic tasks to be used for more customisable analyses!

```sh
run-mtea TIDE dea_file.tsv \
    -o results/tide_results.tsv \
    -n 1000 \
    --n_cpus 4 \
    --or_func absmax \
    --gene_col geneID \
    --lfc_col log2FoldChange \
    --pvalue_col padj \
    -a 0.05 \
    --mask_lfc_values
```

### Understanding the TIDE results

Once the analysis is run, a tabular file containing the analysis results will be saved into the inputed location. The results file will contain 7 columns: a task ID, the metabolic score, the mean random score obtained during the permutation test, its associated p-value, and three more columns detailing the metabolic task description, metabolic system and subsystem.

```
    task_id     score  random_score  pvalue                                 task_description               metabolic_system            metabolic_subsystem
0      X159  1.150921     -0.201177   0.000                           Linolenate degradation              Lipids Metabolism          Fatty Acid Metabolism
1      X164  1.216932     -0.185784   0.001                         Arachidonate degradation              Lipids Metabolism          Fatty Acid Metabolism
2      X160  1.026041     -0.184708   0.001                            Linoleate degradation              Lipids Metabolism          Fatty Acid Metabolism
3      X107  1.228656     -0.178454   0.001         Conversion of lysine to L-2-Aminoadipate         Amino Acids Metabolism              Lysine Metabolism
4      X162  0.857837     -0.182584   0.001                     gamma-Linolenate degradation              Lipids Metabolism          Fatty Acid Metabolism
..      ...       ...           ...     ...                                              ...                            ...                            ...
```

The results can then be used to explore the metabolic changes of a case-control sample.

## The TIDE-essential framework

{: .highlight}
Under construction! Please, come back soon.