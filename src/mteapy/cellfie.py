import re
import pandas as pd
import numpy as np
from cobra.core.gene import GPR
from ast import Name, And, Or, BoolOp, Expression

from cobra.core import model


###########################################
# CellFie functions
###########################################

def calculate_GAL(expr_data:pd.DataFrame, thresh_type:str = "local", local_thresh_type:str = "minmaxmean", minmaxmean_thresh_type:str = "percentile", upper_bound:float = 0.75, lower_bound:float = 0.25, global_thresh_type:str = None, global_value:float = None):
    """
    Function to calculate Gene Activity Levels (GALs) from a gene expression matrix using the CellFie framework. It calculates specific thresholds internally given the inputs from the user.

    Parameters
    ----------
    expr_data: pandas.DataFrame
        A pandas DataFrame containing gene expression values where rows are genes and columns are samples.
    
    thresh_type: str ["local" | "global"] (default: "local")
        Thresholding strategy to use, either locally using specific thresholds for each gene or globally using the same threshold for all genes.
    
    local_thresh_type: str ["minmaxmean" | "mean"] 
        Local thresholding strategy to use. Minmaxmean uses the mean expression value and some upper and lower bounds to compute thresholds. Mean uses the mean of each gene as threshold (default: "minmaxmean").
    
    minmaxmean_thresh_type: str ["value" | "percentile"] 
        Type of upper and lower bounds to apply to the minmaxmean thresholding strategy (default: "percentile").
    
    upper_bound: float 
        Upper bound for minmaxmean thresholding strategy. If using percentiles, value must be between 0 and 1 (default: 0.75).
    
    lower_bound: float 
        Lower bound for minmaxmean thresholding strategy. If using percentiles, value must be between 0 and 1 (default: 0.25).
    
    global_thresh_type: str ["value" | "percentile"] 
        Global thresholding strategy to use. Value will consider a global value as threshold for all genes, and percentile will consider a global percentile as threshold for all genes (default: None).
    
    global_value: float 
        Value to use for global thresholding strategy. If using percentile, value must be between 0 and 1 (default: None).

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


def compute_CellFie(expr_data:pd.DataFrame, task_structure:pd.DataFrame, model:model, thresh_type:str = "local", local_thresh_type:str = "minmaxmean", minmaxmean_thresh_type:str = "percentile", upper_bound:float = 0.75, lower_bound:float = 0.25, global_thresh_type:str = None, global_value:float = None):
    """
    Wrapper function to compute the CellFie framework.

    Parameters
    ----------
    expr_data: pandas.DataFrame
        A pandas DataFrame containing gene expression values where rows are genes and columns are samples. Genes should be inputed as the index of the DataFrame.

    task_structure: pandas.DataFrame
        A boolean matrix where rows are reactions and columns metabolic tasks. Each cell contains ones or zeros, indicating whether a reaction is involved in a metabolic task.
    
    model: cobra.core.model
        A COBRA metabolic model.
    
    thresh_type: str ["local" | "global"] (default: "local")
        Thresholding strategy to use, either locally using specific thresholds for each gene or globally using the same threshold for all genes.
    
    local_thresh_type: str ["minmaxmean" | "mean"] 
        Local thresholding strategy to use. Minmaxmean uses the mean expression value and some upper and lower bounds to compute thresholds. Mean uses the mean of each gene as threshold (default: "minmaxmean").
    
    minmaxmean_thresh_type: str ["value" | "percentile"] 
        Type of upper and lower bounds to apply to the minmaxmean thresholding strategy (default: "percentile").
    
    upper_bound: float 
        Upper bound for minmaxmean thresholding strategy. If using percentiles, value must be between 0 and 1 (default: 0.75).
    
    lower_bound: float 
        Lower bound for minmaxmean thresholding strategy. If using percentiles, value must be between 0 and 1 (default: 0.25).
    
    global_thresh_type: str ["value" | "percentile"] 
        Global thresholding strategy to use. Value will consider a global value as threshold for all genes, and percentile will consider a global percentile as threshold for all genes (default: None).
    
    global_value: float 
        Value to use for global thresholding strategy. If using percentile, value must be between 0 and 1 (default: None).
    
    Returns
    -------
    metabolic_scores_df: pandas.DataFrame
        A pandas DataFrame containing the metabolic scores where rows correspond to metabolic tasks and columns to samples.
    
    binary_scores_df: pandas.DataFrame
        A pandas DataFrame containing the binary metabolic scores after applying a certain threshold of activity. Rows correspond to metabolic tasks and columns to samples.
    """
    task_structure = task_structure.astype(bool)

    # CellFie only uses reactions with valid GPR rules, internally they are considered as 0.0
    gpr_dict = {rxn.id: rxn.gpr for rxn in model.reactions if rxn.gpr != GPR.from_string("")}

    gal_df = calculate_GAL(
        expr_data, 
        thresh_type, 
        local_thresh_type,
        minmaxmean_thresh_type, 
        upper_bound, lower_bound, 
        global_thresh_type, global_value
    )
    ral_df = calculate_RAL(gal_df, gpr_dict)
    metabolic_scores_df, binary_scores_df = calculate_CellFie_scores(ral_df, task_structure)

    return metabolic_scores_df, binary_scores_df