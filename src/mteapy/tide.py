import pandas as pd
import numpy as np
from multiprocessing import Pool

import argparse
from cobra.core import model

from mteapy.utils import safe_eval_gpr, calculate_pvalue, MTEA_parallel_worker 


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
        A pandas DataFrame containing gene expression change values. Rows should correspond to the different genes, and columns should contain at least a gene column, an expression column, and a p-value column.
    
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