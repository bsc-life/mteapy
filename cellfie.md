---
title: CellFie
layout: default
nav_order: 3
---

# **The CellFie framework**
{: .no_toc }
***

## Table of contents
{: .no_toc .text-delta }
1. TOC
{:toc}

## Description

The **CellFie** framework is a contraint-based metabolic modeling framework that was originally published by [Richelle _et al._, 2021](https://doi.org/10.1016/j.crmeth.2021.100040). It leverages the use of mathematical descriptions of metabolic functions (metabolic tasks) and transcriptomics data to quantify metabolic functions. As opposed to TIDE, the CellFie framework allows for the processing of multiple samples at a time, making it suitable for large datasets and single-cell RNA sequencing.

## Command options

| Argument | Shortcut | Description | Default |
|:-------- |:-------- |:----------- |:------- |
| `expr_file` | | Filename for a normalized gene expression file (e.g., TPM). It should contain at least one column with gene names/symbols. | |
| `--delim` | `-d` | Field delimiter for inputed file. | `\t` |
| `-out` | `-o` | Directory to store the analysis' results. The result file(s) will be stored in the specified directory in a tab-sepparated format (`.tsv`). | `CellFie_results/` |
| `--gene_col` | | Name of the column in the inputed file containing gene names/symbols. | `geneID` |
| `--threshold_type` | | Determines the threshold approach to be used. A `global` approach used the same threshold for all genes whereas a `local` approach uses a different threshold for each gene when computing the gene activity levels. | `local` |
| `--global_threshold_type` | | Whether to use a `value` or a `percentile` of the distribution of all genes as global treshold for all genes. | `percentile` |
| `--global_value` | | Value to use as global threshold according to the `global_threshold_type` option selected. Note that percentile values must be between 0 and 1. | `0.75` |
| `--local_threshold_type` | | Determines the threshold type to be used in a local approach. `minmaxmean`: the threshold for each gene is determined by the mean of expression values across all conditions/samples but must be higher or equal than a lower bound and lower or equal to an upper bound. `mean`: the threshold of a gene is determined as its mean expression across all conditions/samples. | `minmaxmean` |
| `--minmaxmean_threshold_type` | | Whether to use `value` or `percentile` of the distribution of all genes as upper and lower bounds. | `percentile` |
| `--upper_bound` | | Upper bound value to be used according to the `minmaxmean_threshold_type`. Note that percentile values must be between 0 and 1. | `0.75` |
| `--lower_bound` | | Lower bound value to be used according to the `minmaxmean_threshold_type`. Note that percentile values must be between 0 and 1. | `0.25` |
| `--binary_scores` | | Flag to indicate whether to also return the binary metabolic score matrix as a second result file. See the original publication for more details | `False` |

## Usage Example

### Transcriptomics Data

One of the first things that the CellFie framework requires is a normalized gene expression matrix (usually stored as TPMs). Normaly, this type of data contains gene names/symbols as rows, and samples as columns. For the command to run, one of the columns of the matrix must store the information regarding gene names/symbols.

A typical normalized gene expression matrix will look like the following:

```
               geneID         S1         S2         S3         S4
0     ENSG00000000419   6.721972   7.768211   0.111999   0.561086
1     ENSG00000001036   5.880123  10.804611   4.273897   3.703098
2     ENSG00000001084  13.568022  11.912389  21.792070   4.126645
3     ENSG00000001630   9.830659  10.973878  16.052115   3.264040
4     ENSG00000002549  10.312642  10.373970   6.246490   0.597024
...               ...        ...        ...        ...        ...
```

### Running CellFie

To run the CellFie framework using the command-line, the command `run-mtea CellFie` should be used with the desired arguments. A typical CellFie analysis is run using the `minmaxmean` local thresholding strategy, which will be used by default by the command, with a percentile upper and lower bounds of `0.75` and `0.25`.

{: .note}
Only the **Human-GEM** and its metabolic tasks are implemented, so the framework will only take in **EnsemblIDs** as valid genic nomenclature. We are working to allow for any metabolic model and metabolic tasks to be used for more customisable analyses!

```sh
run-mtea CellFie expression_file.tsv \
    -o results/ \
    --gene_col geneID \
    --threshold_type local \
    --local_threshold_type minmaxmean \
    --minmaxmean_threshold_type percentile \
    --upper_bound 0.75 \
    --lower_bound 0.25 \
    --binary_scores
```

### Understanding the CellFie results

Once the analysis is run, one or two results files will be stored in the specified directory.

| Result File | Description |
|:----------- |:----------- |
| `cellfie_scores.tsv` | Main result file containing the metabolic activity score values. Columns represent the samples in the original gene expression file, and rows represent all the different metabolic tasks (stored in the `task_id` column). |
| `cellfie_binary_scores.tsv` | Secondary result file that will only be generated if the flag `--binary_scores` is specified. It has the same structure as the main result file, but contains the binary interpretation of the activity of a metabolic task (`0` if the task is considered inactive, `1` if the task is considered active). | 

A standard run of the CellFie framework should produce a `cellfie_scores.tsv` file similar to the following:

```
    task_id        S1        S2        S3        S4
0        X1  0.076515  0.671443  0.050100  0.733470
1        X2  0.863416  0.561653  1.204112  1.253820
2        X3  1.354543  0.889970  1.586738  2.489626
3        X4  1.195976  1.961420  1.423547  1.644596
4        X5  1.554831  1.785477  1.541452  1.704194
..      ...       ...       ...       ...       ...
```

Metabolic tasks are stored using their internal IDs, and their metadata can easily retrieved at the [task_info/](https://github.com/bsc-life/mteapy/tree/main/task_info) folder at the MTEApy repository.