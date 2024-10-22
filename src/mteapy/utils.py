import re
from multiprocessing import Pool

import pandas as pd
import numpy as np

import argparse

from ast import Name, And, Or, BoolOp, Expression
from cobra.core.gene import GPR
from cobra.core import model


###########################################
# Generic functions
###########################################

def absolute_minmax(array:np.ndarray, func:str):
    """
    Internal function to return the index of the maximum or minimum absolute value of an array.

    Parameters
    ----------
    array: list | numpy.ndarray
        Array from which to compute the absolute maximum or minimum.
    
    func: str ["absmax" | "absmin"]
        Wheter to compute the absolute maximum or the minimum.
    
    Returns
    -------
    value: float
        The absolute maximum or minimum.
    """
    abs_array = np.abs(array)
    
    if func == "absmax":
        max_idx = np.argmax(abs_array)
        return array[max_idx]
    
    elif func == "absmin":
        min_idx = np.argmin(abs_array)
        return array[min_idx]
    

def safe_eval_gpr(expr:GPR, gene_dict:dict, or_func:str):
    """
    Recursive function to parse through gene-protein-reaction (GPR) rules.

    Parameters
    ----------
    expr: cobra.core.gene.GPR or ast.Name or ast.BoolOp
        A GPR expression or an abstract syntax tree (AST) object
    
    gene_dict: dict
        A dictionary of genes to their expression values.

    or_func: str ["absmax" | "max"]
        Function to evaluate OR rules within a GPR rule.
    
    Returns
    -------
    gene_score: float
        The expression value of the resolved GPR rule.
    """
    if isinstance(expr, (Expression, GPR)):
        return safe_eval_gpr(expr.body, gene_dict, or_func)
    
    elif isinstance(expr, Name):
        fgid = re.sub(r"\.\d*", "", expr.id)      # Removes "." notation from genes
        return gene_dict.get(fgid, 0)
    
    elif isinstance(expr, BoolOp):
        op = expr.op
        if isinstance(op, Or):
            if or_func == "max":
                return max([safe_eval_gpr(i, gene_dict, or_func) for i in expr.values])
            elif or_func == "absmax":
                return absolute_minmax([safe_eval_gpr(i, gene_dict, or_func) for i in expr.values], func="absmax")
        elif isinstance(op, And):
            return min([safe_eval_gpr(i, gene_dict, or_func) for i in expr.values])
        else:
            raise TypeError("unsupported operation " + op.__class__.__name__)
    
    # If there is no GPR rule, return 0
    elif expr is None:
        return 0
    
    else:
        raise TypeError("unsupported operation " + repr(expr))
    

def calculate_pvalue(score:float, random_scores:np.ndarray, n_permutations:int):
    """
    Function to calculate the empirical p-value as the probability to observe an equal or more extreme metabolic score using the null distributions generated from the random scores.

    Parameters
    ----------
    score: float
        Actual metabolic score for a given metabolic task (test statistic).
    
    random_scores: numpy.ndarray
        An array of random scores that has a length equal to the number of permutations (null distribution).
    
    Returns
    -------
    pvalue: float
        The empirical p-value.
    """
    pvalue = np.minimum(np.sum(random_scores.astype(float) <= float(score)) / n_permutations,
                        np.sum(random_scores.astype(float) >= float(score)) / n_permutations)
    
    return pvalue


def MTEA_parallel_worker(arguments:tuple) -> list:
    """
    Helper function to parallellise the computation of random metabolic scores to inferr significancy.

    Parameters
    ----------
    arguments: tuple
        A tuple containing different arguments needed to compute a random score array.
    
    Returns
    -------
    random_scores: numpy.ndarray
        An array of random metabolic scores in the same order as the columns of the task structure object.
    """
    # Extracting first last argument to know which framework has the user selected
    args = arguments[-1]
    
    if args.command == "TIDE-essential":
        genes, lfc_vector, task_to_gene, random_seed, args = arguments
        np.random.seed(random_seed)
        np.random.shuffle(lfc_vector)
        random_gene_dict = dict(zip(genes, lfc_vector))
        
        random_scores = [np.mean([random_gene_dict.get(gene, 0.0) for gene in task_to_gene[task]]) for task in task_to_gene]
        
        return np.array(random_scores)
    
    elif args.command == "TIDE":
        genes, lfc_vector, task_structure, gpr_dict, random_seed, args = arguments
        np.random.seed(random_seed)
        np.random.shuffle(lfc_vector)
        random_gene_dict = dict(zip(genes, lfc_vector))
        
        random_scores = [safe_eval_gpr(gpr_dict[rxn], random_gene_dict, args.or_func) \
                        for rxn in task_structure.index]
        
        return np.array(random_scores)


###########################################
# TIDE-essential functions
###########################################

def calculate_TIDEe_scores(gene_dict:dict, gene_essentiality:pd.DataFrame):
    """
    Function to calculate the metabolic score from TIDE using only essential genes.

    Parameters
    ----------
    gene_dict: dict
        A dictionary of genes to their expression values.
    
    gene_essentiality: pandas.DataFrame
        A boolean matrix where rows are genes and columns metabolic tasks. Each cell contains ones or zeros, indicating whether a gene is essential for a metabolic task.
    
    Returns
    -------
    scores: numpy.ndarray
        An array of metabolic scores in the same order as the columns of the gene essentiality matrix.
    """
    task_to_gene = {task: gene_essentiality.index[gene_essentiality[task]].to_list() \
                    for task in gene_essentiality.columns}

    scores = [np.mean([gene_dict.get(gene, 0.0) for gene in task_to_gene[task]]) for task in task_to_gene]

    return np.array(scores)


def calculate_random_TIDEe_scores(gene_dict:dict, gene_essentiality:pd.DataFrame, args):
    """
    Function to compute Nº permutations * random metabolic scores to inferr significancy.

    Parameters
    ----------
    gene_dict: dict
        A dictionary of genes to their expression values. 
    
    gene_essentiality: pandas.DataFrame
        A boolean matrix where rows are genes and columns metabolic tasks. Each cell contains ones or zeros, indicating whether a gene is essential for a metabolic task.
    
    args: argparse.ArgumentParser
        The arguments passed by the user through the command line.
    
    Returns
    -------
    random_scores_df: pandas.DataFrame
        A pandas DataFrame where rows are permutations and columns metabolic tasks. Each cell correspond to a random metabolic score for a specific permutation and task.
    """
    genes = list(gene_dict.keys())
    lfc_vector = np.array(list(gene_dict.values()))
    task_to_gene = {task: gene_essentiality.index[gene_essentiality[task]].to_list() \
                    for task in gene_essentiality.columns}

    if args.n_jobs <= 1:
        random_scores = np.zeros((args.n_permutations, len(gene_essentiality.columns)))
        for i in range(args.n_permutations):
            np.random.shuffle(lfc_vector)
            random_gene_dict = dict(zip(genes, lfc_vector))
            random_scores[i,:] = [np.mean([random_gene_dict.get(gene, 0.0) for gene in task_to_gene[task]]) \
                                  for task in task_to_gene]
            
    else:
        # Argument list for parallel processing
        arguments = [(genes, lfc_vector, task_to_gene, np.random.random_integers(0,1e6), args) \
                     for _ in range(args.n_permutations)]
        
        # Parallel execution
        with Pool(processes=args.n_jobs) as pool:
            map_result = pool.map_async(MTEA_parallel_worker, arguments, chunksize=100)
            random_scores = np.array([array for array in map_result.get()])
            
    return pd.DataFrame(random_scores, columns=gene_essentiality.columns)


def compute_TIDEe(expr_data:pd.DataFrame, gene_essentiality:pd.DataFrame, args):
    """
    Wrapper function to compute the TIDE-essential framework.

    Parameters
    ----------
    expr_data: pandas.DataFrame
        A pandas DataFrame containing gene expression values. Rows should correspond to the different genes, and columns should contain at least a gene column, an expression column, and a p-value column.
    
    gene_essentiality: pandas.DataFrame
        A boolean matrix where rows are genes and columns metabolic tasks. Each cell contains ones or zeros, indicating whether a gene is essential for a metabolic task.
    
    args: argparse.ArgumentParser
        The arguments passed by the user through the command line.

    Returns
    -------
    TIDE_e_results: pandas.DataFrame
        A pandas DataFrame where rows are metabolic tasks and columns correspond to their metabolic score from TIDE-essential, their average random score calculated from the null distribution, and their p-value.
    """
    expr_data.set_index(args.gene_col, inplace=True)
    gene_dict = dict(zip(expr_data.index, expr_data[args.lfc_col]))
    gene_essentiality = gene_essentiality.astype(bool)

    scores = calculate_TIDEe_scores(gene_dict, gene_essentiality)
    random_scores_df = calculate_random_TIDEe_scores(gene_dict, gene_essentiality, args)
    pvalues =  [calculate_pvalue(scores[i], random_scores_df[task], args.n_permutations) \
                for i, task in enumerate(gene_essentiality.columns)]
    
    TIDE_e_results = pd.DataFrame({"task_id": gene_essentiality.columns,
                                   "score": scores,
                                   "random_score": random_scores_df.mean(),
                                   "pvalue": pvalues
                                    })
    
    if args.random_scores_flag:
        TIDE_e_results["random_score_array"] = np.array([";".join(random_scores_df.T.loc[task].astype(str)) \
                                                       for task in TIDE_e_results["task_id"]])
    
    return TIDE_e_results


###########################################
# TIDE functions
###########################################

def calculate_TIDE_scores(gene_dict:dict, task_structure:pd.DataFrame, gpr_dict:dict, args:argparse.ArgumentParser):
    """
    Function to calculate the metabolic score from TIDE.

    Parameters
    ----------
    gene_dict: dict
        A dictionary of genes to their expression values. 
    
    task_structure: pandas.DataFrame
        A boolean matrix where rows are reactions and columns metabolic tasks. Each cell contains ones or zeros, indicating whether a reaction is involved in a metabolic task.
    
    gpr_dict: dict
        A dictionary of reactions to their GPR rules. GPR rules must be cobra.core.gene.GPR or ast.Name or ast.BoolOp objects.
    
    args: argparse.ArgumentParser
        The arguments passed by the user through the command line.

    Returns
    -------
    scores: numpy.ndarray
        An array of metabolic scores in the same order as the columns of the task structure object.
    """
    rxn_projection = {rxn: safe_eval_gpr(gpr_dict[rxn], gene_dict, args.or_func) for rxn in task_structure.index}
    task_to_rxns = {task: task_structure.index[task_structure[task]].to_list() for task in task_structure.columns}
    
    scores = [np.mean([rxn_projection[rxn] for rxn in task_to_rxns[task]]) for task in task_to_rxns]
    
    return np.array(scores)


def calculate_random_TIDE_scores(gene_dict:dict, task_structure:pd.DataFrame, gpr_dict:dict, args:argparse.ArgumentParser):
    """
    Function to compute Nº permutations * random metabolic scores to inferr significancy.

    Parameters
    ----------
    gene_dict: dict
        A dictionary of genes to their expression values. 
    
    task_structure: pandas.DataFrame
        A boolean matrix where rows are reactions and columns metabolic tasks. Each cell contains ones or zeroes, indicating whether a reaction is involved in a metabolic task.
    
    gpr_dict: dict
        A dictionary of reactions to their GPR rules. GPR rules must be cobra.core.gene.GPR or ast.Name or ast.BoolOp objects.
    
    args: argparse.ArgumentParser
        The arguments passed by the user through the command line.
    
    Returns
    -------
    random_scores_df: pandas.DataFrame
        A pandas DataFrame where rows are permutations and columns metabolic tasks. Each cell correspond to a random metabolic score for a specific permutation and task.
    """
    genes = list(gene_dict.keys())
    lfc_vector = np.array(list(gene_dict.values()))

    # Decide whether to parallellise
    if args.n_jobs <= 1:
        random_projection = np.zeros((args.n_permutations, len(task_structure.index)))
        for i in range(args.n_permutations):
            np.random.shuffle(lfc_vector)
            random_gene_dict = dict(zip(genes, lfc_vector))
            random_projection[i,:] = [safe_eval_gpr(gpr_dict[rxn], random_gene_dict, args.or_func) \
                                      for rxn in task_structure.index]

    else:
        # Argument list for parallel processing
        arguments = [(genes, lfc_vector, task_structure, gpr_dict, np.random.random_integers(0,1e6), args) \
                    for _ in range(args.n_permutations)]
        
        # Parallel execution
        with Pool(processes=args.n_jobs) as pool:
            map_result = pool.map_async(MTEA_parallel_worker, arguments, chunksize=100)
            random_projection = np.array([array for array in map_result.get()])
    
    random_projection_df = pd.DataFrame(random_projection, columns=task_structure.index).T

    # Reduce from random reaction projections to random metabolic scores
    # (n_permutations * n_reactions) --> (n_permutations * n_tasks)
    random_scores = np.zeros((args.n_permutations, len(task_structure.columns)))
    for i, task in enumerate(task_structure.columns):
        rxns = [k for k in task_structure.index[task_structure[task]] if k in random_projection_df.index]
        random_scores[:,i] = random_projection_df.loc[rxns].mean().values

    return pd.DataFrame(random_scores, columns=task_structure.columns)


def compute_TIDE(expr_data:pd.DataFrame, task_structure:pd.DataFrame, model:model, args:argparse.ArgumentParser):
    """
    Wrapper function to compute the TIDE framework.

    Parameters
    ----------
    expr_data: pandas.DataFrame
        A pandas DataFrame containing gene expression values. Rows should correspond to the different genes, and columns should contain at least a gene column, an expression column, and a p-value column.
    
    task_structure: pandas.DataFrame
        A boolean matrix where rows are reactions and columns metabolic tasks. Each cell contains ones or zeros, indicating whether a reaction is involved in a metabolic task.
    
    model: cobra.core.model
        A COBRA metabolic model.
    
    args: argparse.ArgumentParser
        The arguments passed by the user through the command line.
    
    Returns
    -------
    TIDE_results: pandas.DataFrame
        A pandas DataFrame where rows are metabolic tasks and columns correspond to their score from TIDE, their average random score calculated from the null distribution, and their p-value.
    """
    expr_data.set_index(args.gene_col, inplace=True)
    gene_dict = dict(zip(expr_data.index, expr_data[args.lfc_col]))
    task_structure = task_structure.astype(bool)
    gpr_dict = {rxn.id: rxn.gpr for rxn in model.reactions if rxn.id in task_structure.index}

    scores = calculate_TIDE_scores(gene_dict, task_structure, gpr_dict, args)
    random_scores_df = calculate_random_TIDE_scores(gene_dict, task_structure, gpr_dict, args)
    pvalues = [calculate_pvalue(scores[i], random_scores_df[task], args.n_permutations) \
               for i, task in enumerate(task_structure.columns)]

    TIDE_results = pd.DataFrame({"task_id": task_structure.columns,
                                 "score": scores,
                                 "random_score": random_scores_df.mean(),
                                 "pvalue": pvalues
                                })
    
    if args.random_scores_flag:
        TIDE_results["random_score_array"] = np.array([";".join(random_scores_df.T.loc[task].astype(str)) \
                                                       for task in TIDE_results["task_id"]])

    return TIDE_results


###########################################
# CellFie functions
###########################################

def calculate_GAL(expr_data:pd.DataFrame, thresh_type:str, local_thresh_type:str, minmaxmean_thresh_type:str, upper_bound:float, lower_bound:float, global_thresh_type:str, global_value:float):
    """
    Function to calculate Gene Activity Levels (GALs) from a gene expression matrix using the CellFie framework. It calculates specific thresholds internally given the inputs from the user.

    Parameters
    ----------
    expr_data: pandas.DataFrame
        A pandas DataFrame containing gene expression values where rows are genes and columns are samples.
    
    thresh_type: str ["local" | "global"]
        Thresholding strategy to use, either locally using specific thresholds for each gene or globally using the same threshold for all genes.
    
    local_thresh_type: str ["minmaxmean" | "mean"]
        Local thresholding strategy to use. Minmaxmean uses the mean expression value and some upper and lower bounds to compute thresholds. Mean uses the mean of each gene as threshold.
    
    minmaxmean_thresh_type: str ["value" | "percentile"]
        Type of upper and lower bounds to apply to the minmaxmean thresholding strategy.
    
    upper_bound: float
        Upper bound for minmaxmean thresholding strategy. If using percentiles, value must be between 0 and 1.
    
    lower_bound: float
        Lower bound for minmaxmean thresholding strategy. If using percentiles, value must be between 0 and 1.
    
    global_thresh_type: str ["value" | "percentile"]
        Global thresholding strategy to use. Value will consider a global value as threshold for all genes, and percentile will consider a global percentile as threshold for all genes.
    
    global_value: float
        Value to use for global thresholding strategy. If using percentile, value must be between 0 and 1.

    Returns
    -------
    gal_df: pandas.DataFrame
        A pandas DataFrame containing GALs.
    
    Notes
    -----
    GALs and thresholds are computed as indicated in the CellFie paper.
    """
    # Transform into numpy and pre-allocate array for results
    expr_data_array = expr_data.to_numpy()
    
    # Delete zero values to compute thresholds
    to_keep = expr_data_array != 0
    linear_data = expr_data_array[to_keep]
    
    gene_scores = np.zeros((len(expr_data.index), len(expr_data.columns)))
    
    # Local approach
    if thresh_type == "local":
        if local_thresh_type == "mean":
            thresholds = np.mean(expr_data_array, axis=1)
        
        elif local_thresh_type == "minmaxmean":
            if minmaxmean_thresh_type == "value":
                local_upper = upper_bound
                local_lower = lower_bound
            elif minmaxmean_thresh_type == "percentile":
                local_upper = np.quantile(linear_data, upper_bound)
                local_lower = np.quantile(linear_data, lower_bound)
            
            expr_mean = np.mean(expr_data_array, axis=1)
            thresholds = np.clip(expr_mean, local_lower, local_upper)

        gene_scores = 5 * np.log(1 + expr_data_array / thresholds[:,np.newaxis])
    
    # Global approach
    elif thresh_type == "global":
        if global_thresh_type == "value":
            global_threshold = global_value
        
        elif global_thresh_type == "percentile":
            global_threshold = np.quantile(linear_data, global_value)
        
        gene_scores = 5 * np.log(1 + expr_data_array / global_threshold)
        
    return pd.DataFrame(gene_scores, index=expr_data.index, columns=expr_data.columns)


def safe_eval_gpr_w_names(expr:GPR, conf_genes:dict):
    """
    Internal function to evaluate a gene-protein rule in an injection-safe manner (hopefully).
    """
    if isinstance(expr, (Expression, GPR)):
        return safe_eval_gpr_w_names(expr.body, conf_genes)
    
    elif isinstance(expr, Name):
        fgid = re.sub(r"\.\d*", "", expr.id)      # Removes "." notation from genes
        return conf_genes.get(fgid, 0), fgid
    
    elif isinstance(expr, BoolOp):
        op = expr.op
        evaluated_values = [safe_eval_gpr_w_names(i, conf_genes) for i in expr.values]
        filtered_values = [(value, gene) for value, gene in evaluated_values if gene in conf_genes]
        # Return default values if no valid genes found
        if len(filtered_values) == 0:
            return 0, "0" 
        if isinstance(op, Or):
            return max(filtered_values, key=lambda x: x[0])
        elif isinstance(op, And):
            return min(filtered_values, key=lambda x: x[0])
        else:
            raise TypeError("unsupported operation " + op.__class__.__name__)
    
    elif expr is None:
        return 0, "0"
    
    else:
        raise TypeError("unsupported operation " + repr(expr))


def calculate_RAL(gal_df:pd.DataFrame, gpr_dict:dict):
    """
    Function to calculate Reaction Activity Levels (RALs) from a Gene Activity Level matrix (GALs) using the CellFie framework.

    Parameters
    ----------
    gal_df: pandas.DataFrame
        A pandas DataFrame containing the GALs.

    gpr_dict: dict
        A dictionary of reactions to their GPR rules. GPR rules must be cobra.core.gene.GPR or ast.Name or ast.BoolOp objects.

    Returns
    -------
    ral_df: pandas.DataFrame
        A pandas DataFrame containing RALs where rows are reactions and columns are samples.
    
    Notes
    -----
    RALs are computed as indicated in the CellFie paper.
    """
    gal_array = gal_df.to_numpy()
    all_reactions = list(gpr_dict.keys())
    all_genes = gal_df.index

    # Pre-allocate ral array
    ral = np.zeros((len(all_reactions),len(gal_df.columns)))

    for i in range(len(gal_df.columns)):
        
        gene_dict = dict(zip(all_genes, gal_array[:,i]))
        rxn_projection = np.array([safe_eval_gpr_w_names(gpr_dict[rxn], gene_dict) for rxn in all_reactions])
        
        # Adjusting RAL by the times a gene appears as the reaction projection within a sample (significancy)
        gene_frequency = np.unique(rxn_projection[:,1], return_counts=True)
        significance_dict = {gene: 1/k for gene, k in zip(gene_frequency[0], gene_frequency[1])}
        significance_array = np.array([significance_dict[gene] for gene in rxn_projection[:,1]])

        ral[:,i] = rxn_projection[:,0].astype(float) * significance_array

    return pd.DataFrame(ral, index=all_reactions, columns=gal_df.columns)


def calculate_CellFie_scores(ral_df:pd.DataFrame, task_structure:pd.DataFrame):
    """
    Function to compute the metabolic scores according to the CellFie framework.

    Parameters
    ----------
    ral_df: pandas.DataFrame
        A pandas DataFrame containing the RALs.
    
    task_structrure: pandas.DataFrame
        A boolean matrix where rows are reactions and columns metabolic tasks. Each cell contains ones or zeros, indicating whether a reaction is involved in a metabolic task.
    
    Returns
    -------
    metabolic_scores_df: pandas.DataFrame
        A pandas DataFrame containing the metabolic scores where rows correspond to metabolic tasks and columns to samples.
    
    binary_scores_df: pandas.DataFrame
        A pandas DataFrame containing the binary metabolic scores after applying a certain threshold of activity. Rows correspond to metabolic tasks and columns to samples.
    
    Notes
    -----
    Metabolic scores and binary metabolic scores are computed as indicated in the CellFie paper.
    """
    task_to_rxns = {task: task_structure.index[task_structure[task]].to_list() for task in task_structure.columns}
    metabolic_scores = np.zeros((len(task_to_rxns.keys()), len(ral_df.columns)))

    for i, sample in enumerate(ral_df.columns):
        rxn_dict = dict(zip(ral_df.index, ral_df[sample]))
        metabolic_scores[:,i] = [np.mean([rxn_dict.get(rxn, 0.0) for rxn in task_to_rxns[task]]) \
                                 for task in task_to_rxns]
    
    metabolic_scores_df = pd.DataFrame(metabolic_scores, index=task_to_rxns.keys(), columns=ral_df.columns)
    binary_scores_df = (metabolic_scores_df >= 5 * np.log(2)).astype(int)

    return metabolic_scores_df, binary_scores_df


def compute_CellFie(expr_data:pd.DataFrame, task_structure:pd.DataFrame, model:model, args:argparse.ArgumentParser):
    """
    Wrapper function to compute the CellFie framework.

    Parameters
    ----------
    expr_data: pandas.DataFrame
        A pandas DataFrame containing gene expression values where rows are genes and columns are samples.

    task_structure: pandas.DataFrame
        A boolean matrix where rows are reactions and columns metabolic tasks. Each cell contains ones or zeros, indicating whether a reaction is involved in a metabolic task.
    
    model: cobra.core.model
        A COBRA metabolic model.
    
    args: argparse.ArgumentParser
        The arguments passed by the user through the command line.
    
    Returns
    -------
    metabolic_scores_df: pandas.DataFrame
        A pandas DataFrame containing the metabolic scores where rows correspond to metabolic tasks and columns to samples.
    
    binary_scores_df: pandas.DataFrame
        A pandas DataFrame containing the binary metabolic scores after applying a certain threshold of activity. Rows correspond to metabolic tasks and columns to samples.
    """
    # expr_data.set_index(args.gene_col, inplace=True)
    task_structure = task_structure.astype(bool)

    # CellFie only uses reactions with valid GPR rules, internally they are considered as 0.0
    gpr_dict = {rxn.id: rxn.gpr for rxn in model.reactions if rxn.gpr != GPR.from_string("")}

    gal_df = calculate_GAL(expr_data, args.thresh_type, args.local_thresh_type, args.minmaxmean_thresh_type,
                        args.upper_bound, args.lower_bound, args.global_thresh_type, args.global_value)
    
    ral_df = calculate_RAL(gal_df, gpr_dict)
    metabolic_scores_df, binary_scores_df = calculate_CellFie_scores(ral_df, task_structure)

    return metabolic_scores_df, binary_scores_df
