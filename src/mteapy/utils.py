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
        on="task_id",
        how="left"
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
    Recursive function to parse a GPR rule for CellFie's RAL computation,
    returning (score, gene_used). `score` is None for "no data" -- a gene
    absent from `conf_genes` -- matching the original MATLAB CellFie's `-1`
    sentinel (`findUsedGenesLevels_all.m`/`selectGeneFromGPR_all.m`): a
    missing gene still participates in AND/OR via that sentinel rather than
    being silently dropped as if the rule never mentioned it. AND (a
    complex's subunits) propagates None if *any* subunit is missing -- you
    can't confirm a complex is active without knowing all of it, matching
    `min` naturally sorting a sentinel-below-every-real-value first in the
    original. OR (isoenzymes) recovers past a missing option as long as
    another one has data, matching `max` ignoring it.

    Callers (`calculate_RAL`) must treat a None score as "no data for this
    reaction/sample", to be excluded from downstream aggregation --
    substituting a plain 0 instead would silently pull scores toward zero
    for any reaction with an unmeasured gene, which is not what CellFie's
    published algorithm does.
    """
    if isinstance(expr, (Expression, GPR)):
        return map_gpr_w_names(expr.body, conf_genes)

    elif isinstance(expr, Name):
        fgid = re.sub(r"\.\d*", "", expr.id)      # Removes "." notation from genes
        if fgid in conf_genes:
            return conf_genes[fgid], fgid
        return None, fgid

    elif isinstance(expr, BoolOp):
        op = expr.op
        evaluated_values = [map_gpr_w_names(i, conf_genes) for i in expr.values]
        if isinstance(op, Or):
            real_values = [(value, gene) for value, gene in evaluated_values if value is not None]
            if not real_values:
                return None, evaluated_values[0][1]
            return max(real_values, key=lambda x: x[0])
        elif isinstance(op, And):
            missing = [(value, gene) for value, gene in evaluated_values if value is None]
            if missing:
                return None, missing[0][1]
            return min(evaluated_values, key=lambda x: x[0])
        else:
            raise TypeError("unsupported operation " + op.__class__.__name__)

    elif expr is None:
        return None, "0"

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


def parallel_chunksize(n_permutations:int, n_cpus:int) -> int:
    """A `Pool.map_async` chunksize that actually distributes work across
    `n_cpus` workers, aiming for ~4 chunks per worker so the run stays
    reasonably load-balanced. A fixed chunksize (e.g. 100) silently defeats
    parallelization entirely whenever n_permutations doesn't far exceed it
    -- with 32 permutations and chunksize=100, every permutation lands in
    one chunk sent to a single worker, and n_cpus>1 buys nothing.
    """
    return max(1, n_permutations // max(1, n_cpus * 4))


_WORKER_STATE: dict = {}


def init_MTEA_parallel_worker(state: dict) -> None:
    """`multiprocessing.Pool(initializer=..., initargs=(state,))` target:
    stashes read-only shared state (gene list, GPR/route/complex data) once
    per *worker process*, instead of it being rebuilt into every
    permutation's own argument tuple and re-pickled/sent over IPC
    `n_permutations` times. For framework="TIDE", this is also where the
    (unpicklable-as-is, or at least awkward to trust across processes)
    string-encoded GPR rules get reconstructed into real `GPR` objects --
    once per worker instead of once per permutation.

    `state` must have a "framework" key (one of "TIDE-essential", "TIDE",
    "TIDE-context-aware") plus that framework's own fields -- see
    `MTEA_parallel_worker`'s branches, and the `calculate_random_*`
    functions in `mteapy.tide` that build `state` and open the `Pool`.
    """
    global _WORKER_STATE
    state = dict(state)
    if state["framework"] == "TIDE":
        state["gpr_dict"] = {
            rxn_id: GPR.from_string(gpr_str) if gpr_str else None
            for rxn_id, gpr_str in state.pop("gpr_string_dict").items()
        }
    _WORKER_STATE = state


def MTEA_parallel_worker(random_seed:int) -> list:
    """
    Computes one permutation's random score array, using whichever
    framework's state an earlier `init_MTEA_parallel_worker` call stashed
    in this worker process.

    Parameters
    ----------
    random_seed: int
        This permutation's random seed.

    Returns
    -------
    random_scores: numpy.ndarray
        An array of random metabolic scores in the same order as the columns of the task structure object.
    """
    state = _WORKER_STATE
    framework = state["framework"]
    np.random.seed(random_seed)
    lfc_vector = state["lfc_vector"].copy()  # never mutate the shared array in place
    np.random.shuffle(lfc_vector)
    random_gene_dict = dict(zip(state["genes"], lfc_vector))

    if framework == "TIDE-essential":
        task_to_gene = state["task_to_gene"]
        random_scores = [np.mean([random_gene_dict.get(gene, 0.0) for gene in task_to_gene[task]]) for task in task_to_gene]

        return np.array(random_scores)

    elif framework == "TIDE":
        gpr_dict = state["gpr_dict"]
        or_func = state["or_func"]
        random_scores = [map_gpr(gpr_dict[rxn], random_gene_dict, or_func) for rxn in state["reaction_ids"]]

        return np.array(random_scores)

    elif framework == "TIDE-context-aware":
        from mteapy.context_scoring import score_task  # local import: avoids a module-load-order cycle with mteapy.tide

        tasks_routes = state["tasks_routes"]
        complex_cache = state["complex_cache"]
        or_func = state["or_func"]
        # Fresh per permutation (this worker call = one permutation's
        # gene_dict), shared across every task in it -- see score_task's
        # `reaction_cache` doc for why this matters.
        reaction_cache: dict = {}
        random_scores = [
            score_task(routes, complex_cache, random_gene_dict, aggregation="mean", task_id=task_id, or_func=or_func,
                       reaction_cache=reaction_cache).score
            for task_id, routes in tasks_routes.items()
        ]

        return np.array(random_scores)

    else:
        raise TypeError(f"Framework {framework} not available for parallelization.")


# def check_model_compatibility(model, structure, type:str = "reactions"):
#     # TODO: Function to check if a task structure/gene essentiality matrix is compatible with a metabolic model
#     pass