import re
import pandas as pd
import numpy as np

from ast import Name, And, Or, BoolOp, Expression
from cobra.core.gene import GPR



###########################################
# Checkpoint functions
###########################################

def print_banner():
    """Function to print MTEApy banner"""
    print(r"""
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
               
    """)

###########################################
# Checkpoint functions
###########################################

def check_ensemblid(gene_array: np.ndarray):
    """
    Function to check that all genes are valid EnsemblIDs.
    
    Parameters
    ----------
    gene_array: numpy.ndarray
        An array of gene names/symbols to be checked.
    
    Returns
    -------
    check: bool
        Whether all genes are EnsemblIDs (True) or not (False)
    """
    checked_array = [gene.startswith("ENSG") for gene in gene_array]
    return all(checked_array)

###########################################
# Generic functions
###########################################

def add_task_metadata(results_df:pd.DataFrame, task_metadata_df:pd.DataFrame):
    """
    Helper function to format and add to any analysis results information regarding the metabolic tasks. The information will be added as three new columns: a description of the metabolic task, its metabolic system and subsystem.

    Parameters
    ----------
    results_df: pandas.DataFrame
        A data frame containing the results of an analysis. This data frame must contain a column named 'task_id' with the internal IDs of the metabolic tasks.
    
    task_metadata: pandas.DataFrame
        The data frame containing the metadata of metabolic tasks. Its formatting must be the same as the file stored in the `task_info/` directory in this repository.

    Returns
    -------
    results_df_annotated: pandas.DataFrame
        The original results data frame with the three new columns added.
    """
    task_metadata_df = task_metadata_df.rename(columns={
        "ID": "task_id", 
        "SYSTEM": "metabolic_system", 
        "DESCRIPTION": "task_description", 
        "SUBSYSTEM": "metabolic_subsystem"
    })
    task_metadata_df["metabolic_system"] = [system.title() for system in task_metadata_df["metabolic_system"]]
    task_metadata_df["metabolic_subsystem"] = [system.title() for system in task_metadata_df["metabolic_subsystem"]]
    
    results_df = results_df.merge(
        task_metadata_df[["task_id","task_description","metabolic_system","metabolic_subsystem"]], 
        on="task_id"
    )
    return results_df


def mask_lfc_values(expr_df:pd.DataFrame, lfc_col:str, pvalue_col:str, alpha:float):
    """
    Function that "masks" non-significant log-FC values to 0.

    Parameters
    ----------
    expr_df: pandas.DataFrame
        A pandas DataFrame containing gene expression change values. Rows should correspond to the different genes, and columns should contain at least a gene column, an expression column, and a p-value column.
    
    lfc_col: str
        Name of the column in expr_df with log-FC values.

    pvalue_col
        Name of the column in expr_df with p-value values.

    alpha: float
        Significance threshold.
    
    Returns
    -------
    masked_expr_df: pandas.DataFrame
        Masked pandas DataFrame with gene expression change values.
    """
    masked_expr_df = expr_df.copy()
    masked_expr_df[pvalue_col] = expr_df[pvalue_col].fillna(1)
    masked_expr_df.loc[masked_expr_df[pvalue_col] >= alpha, lfc_col] = 0

    return masked_expr_df



def absmax(array:np.ndarray):
    """
    Function to return the index of the maximum absolute value of an array.

    Parameters
    ----------
    array: list | numpy.ndarray
        Array from which to compute the absolute maximum or minimum.
    
    Returns
    -------
    value: float
        The absolute maximum of the array.
    """
    abs_array = np.abs(array)
    max_idx = np.argmax(abs_array)
    return array[max_idx]

    

def map_gpr(expr:GPR, gene_dict:dict, or_func:str = "absmax"):
    """
    Recursive function to parse through gene-protein-reaction (GPR) rules.

    Parameters
    ----------
    expr: cobra.core.gene.GPR or ast.Name or ast.BoolOp
        A GPR expression or an abstract syntax tree (AST) object
    
    gene_dict: dict
        A dictionary of genes to their expression values.

    or_func: str ["absmax" | "max"]
        Function to evaluate OR rules within a GPR rule (default: absmax).
    
    Returns
    -------
    gene_score: float
        The expression value of the resolved GPR rule.
    """
    if isinstance(expr, (Expression, GPR)):
        return map_gpr(expr.body, gene_dict, or_func)
    
    elif isinstance(expr, Name):
        fgid = re.sub(r"\.\d*", "", expr.id)      # Removes "." notation from genes
        return gene_dict.get(fgid, 0)
    
    elif isinstance(expr, BoolOp):
        op = expr.op
        if isinstance(op, Or):
            if or_func == "max":
                return max([map_gpr(i, gene_dict, or_func) for i in expr.values])
            elif or_func == "absmax":
                return absmax([map_gpr(i, gene_dict, or_func) for i in expr.values])
            else:
                raise TypeError(f"Unsupported OR function ({or_func}). Please, use absmax or max.")
        elif isinstance(op, And):
            return min([map_gpr(i, gene_dict, or_func) for i in expr.values])
        else:
            raise TypeError("unsupported operation " + op.__class__.__name__)
    
    # If there is no GPR rule, return 0
    elif expr is None:
        return 0
    
    else:
        raise TypeError("unsupported operation " + repr(expr))
    

def map_gpr_w_names(expr:GPR, conf_genes:dict):
    """
    Internal function to evaluate a gene-protein rule in an injection-safe manner (hopefully).
    """
    if isinstance(expr, (Expression, GPR)):
        return map_gpr_w_names(expr.body, conf_genes)
    
    elif isinstance(expr, Name):
        fgid = re.sub(r"\.\d*", "", expr.id)      # Removes "." notation from genes
        return conf_genes.get(fgid, 0), fgid
    
    elif isinstance(expr, BoolOp):
        op = expr.op
        evaluated_values = [map_gpr_w_names(i, conf_genes) for i in expr.values]
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


def calculate_pvalue(score:float, random_scores:np.ndarray):
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
    pvalue = np.minimum(np.sum(random_scores.astype(float) <= float(score)) / len(random_scores),
                        np.sum(random_scores.astype(float) >= float(score)) / len(random_scores))
    
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
    framework = arguments[-1]
    
    if framework == "TIDE-essential":
        genes, lfc_vector, task_to_gene, random_seed, _ = arguments
        np.random.seed(random_seed)
        np.random.shuffle(lfc_vector)
        random_gene_dict = dict(zip(genes, lfc_vector))
        
        random_scores = [np.mean([random_gene_dict.get(gene, 0.0) for gene in task_to_gene[task]]) for task in task_to_gene]
        
        return np.array(random_scores)
    
    elif framework == "TIDE":
        genes, lfc_vector, task_structure, gpr_string_dict, random_seed, or_func, _ = arguments
        np.random.seed(random_seed)
        np.random.shuffle(lfc_vector)
        random_gene_dict = dict(zip(genes, lfc_vector))
        
        # Convert string GPR rules back to GPR objects
        gpr_dict = {rxn_id: GPR.from_string(gpr_str) if gpr_str else None 
                   for rxn_id, gpr_str in gpr_string_dict.items()}
        
        random_scores = [map_gpr(gpr_dict[rxn], random_gene_dict, or_func) \
                        for rxn in task_structure.index]
        
        return np.array(random_scores)
    
    else:
        raise TypeError(f"Framework {framework} not available for parallelization.")


# def check_model_compatibility(model, structure, type:str = "reactions"):
#     # TODO: Function to check if a task structure/gene essentiality matrix is compatible with a metabolic model
#     pass