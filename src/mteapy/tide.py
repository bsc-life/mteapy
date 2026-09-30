import pandas as pd
import numpy as np
from multiprocessing import Pool

from cobra.core import Model

from mteapy.context_scoring import score_task
from mteapy.utils import map_gpr, calculate_pvalue, MTEA_parallel_worker


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


def calculate_random_TIDEe_scores(
        gene_dict:dict, 
        gene_essentiality:pd.DataFrame, 
        n_permutations:int = 1000,
        n_cpus:int = 1,
        random_seed:int = None
    ):
    """
    Function to compute Nº permutations * random metabolic scores to inferr significancy.

    Parameters
    ----------
    gene_dict: dict
        A dictionary of genes to their expression values. 
    
    gene_essentiality: pandas.DataFrame
        A boolean matrix where rows are genes and columns metabolic tasks. Each cell contains ones or zeros, indicating whether a gene is essential for a metabolic task.
    
    n_permutations: int 
        Number of permutations to inferr significancy (default: 1000).
    
    n_cpus: int 
        Number of jobs to parallelize the computation. If 1, the run is not parallellized (default: 1).
    
    Returns
    -------
    random_scores_df: pandas.DataFrame
        A pandas DataFrame where rows are permutations and columns metabolic tasks. Each cell correspond to a random metabolic score for a specific permutation and task.
    """
    genes = list(gene_dict.keys())
    lfc_vector = np.array(list(gene_dict.values()))
    task_to_gene = {task: gene_essentiality.index[gene_essentiality[task]].to_list() \
                    for task in gene_essentiality.columns}

    # Set random seed for reproducibility
    if random_seed is not None:
        np.random.seed(random_seed)

    if n_cpus <= 1:
        random_scores = np.zeros((n_permutations, len(gene_essentiality.columns)))
        for i in range(n_permutations):
            np.random.shuffle(lfc_vector)
            random_gene_dict = dict(zip(genes, lfc_vector))
            random_scores[i,:] = [np.mean([random_gene_dict.get(gene, 0.0) for gene in task_to_gene[task]]) \
                                  for task in task_to_gene]
            
    else:
        # Generate deterministic seeds for parallel processing
        if random_seed is not None:
            seeds = [random_seed + i for i in range(n_permutations)]
        else:
            seeds = [np.random.randint(0, 1_000_000) for _ in range(n_permutations)]
        
        # Argument list for parallel processing
        arguments = [(genes, lfc_vector, task_to_gene, seeds[i], "TIDE-essential") \
                     for i in range(n_permutations)]
        
        # Parallel execution
        with Pool(processes=n_cpus) as pool:
            map_result = pool.map_async(MTEA_parallel_worker, arguments, chunksize=100)
            random_scores = np.array([array for array in map_result.get()])
            
    return pd.DataFrame(random_scores, columns=gene_essentiality.columns)


def compute_TIDEe(
        expr_data:pd.DataFrame, 
        lfc_col:str,
        gene_essentiality:pd.DataFrame,
        n_permutations:int = 1000,
        n_cpus:int = 1,
        random_scores_flag:bool = False,
        random_seed:int = None
    ):
    """
    Wrapper function to compute the TIDE-essential framework.

    Parameters
    ----------
    expr_data: pandas.DataFrame
        A pandas DataFrame containing gene expression change values. Row indices should correspond to the different genes, and columns should contain at least an expression column and a p-value column.

    lfc_col: str
        Column name in expr_data containing log-FC values.
    
    gene_essentiality: pandas.DataFrame
        A boolean matrix where rows are genes and columns metabolic tasks. Each cell contains ones or zeros, indicating whether a gene is essential for a metabolic task.

    n_permutations: int 
        Number of permutations to inferr significancy (default: 1000).
    
    n_cpus: int 
        Number of jobs to parallelize the computation. If 1, the run is not parallellized (default: 1).

    random_scores_flag: bool
        Flag to indicate whether to add to the results DataFrame the null distributions (default: False).
    
    random_seed: int
        Random seed for reproducibility. If None, random seeds will be generated (default: None).
    
    Returns
    -------
    TIDE_e_results: pandas.DataFrame
        A pandas DataFrame where rows are metabolic tasks and columns correspond to their metabolic score from TIDE-essential, their average random score calculated from the null distribution, and their p-value.
    """
    gene_dict = dict(zip(expr_data.index, expr_data[lfc_col]))
    gene_essentiality = gene_essentiality.astype(bool)

    scores = calculate_TIDEe_scores(gene_dict, gene_essentiality)
    random_scores_df = calculate_random_TIDEe_scores(gene_dict, gene_essentiality, n_permutations, n_cpus, random_seed)
    pvalues =  [calculate_pvalue(scores[i], random_scores_df[task]) for i, task in enumerate(gene_essentiality.columns)]
    
    TIDE_e_results = pd.DataFrame({"task_id": gene_essentiality.columns,
                                   "score": scores,
                                   #"random_score": random_scores_df.mean(),
                                   "pvalue": pvalues
                                    })
    
    if random_scores_flag:
        TIDE_e_results["random_score_array"] = np.array([";".join(random_scores_df.T.loc[task].astype(str)) \
                                                       for task in TIDE_e_results["task_id"]])
    
    return TIDE_e_results


###########################################
# TIDE functions
###########################################

def calculate_TIDE_scores(
        gene_dict:dict, 
        task_structure:pd.DataFrame, 
        gpr_dict:dict, 
        or_func:str
    ):
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
    
    or_func: str ["absmax" | "max"]
        Function to evaluate OR rules within a GPR rule.

    Returns
    -------
    scores: numpy.ndarray
        An array of metabolic scores in the same order as the columns of the task structure object.
    """
    rxn_projection = {rxn: map_gpr(gpr_dict[rxn], gene_dict, or_func) for rxn in task_structure.index}
    task_to_rxns = {task: task_structure.index[task_structure[task]].to_list() for task in task_structure.columns}
    
    scores = [np.mean([rxn_projection[rxn] for rxn in task_to_rxns[task]]) for task in task_to_rxns]
    
    return np.array(scores)


def calculate_random_TIDE_scores(
        gene_dict:dict, 
        task_structure:pd.DataFrame, 
        gpr_dict:dict, 
        or_func:str,
        n_permutations:int = 1000,
        n_cpus:int = 1,
        random_seed:int = None
    ):
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
    
    or_func: str ["absmax" | "max"]
        Function to evaluate OR rules within a GPR rule.

    n_permutations: int 
        Number of permutations to inferr significancy (default: 1000).
    
    n_cpus: int 
        Number of jobs to parallelize the computation. If 1, the run is not parallellized (default: 1).
    
    random_seed: int
        Random seed for reproducibility. If None, random seeds will be generated (default: None).
    
    Returns
    -------
    random_scores_df: pandas.DataFrame
        A pandas DataFrame where rows are permutations and columns metabolic tasks. Each cell correspond to a random metabolic score for a specific permutation and task.
    """
    genes = list(gene_dict.keys())
    lfc_vector = np.array(list(gene_dict.values()))

    # Set random seed for reproducibility
    if random_seed is not None:
        np.random.seed(random_seed)

    # Decide whether to parallellise
    if n_cpus <= 1:
        random_projection = np.zeros((n_permutations, len(task_structure.index)))
        for i in range(n_permutations):
            np.random.shuffle(lfc_vector)
            random_gene_dict = dict(zip(genes, lfc_vector))
            random_projection[i,:] = [map_gpr(gpr_dict[rxn], random_gene_dict, or_func) \
                                      for rxn in task_structure.index]

    else:
        # Convert GPR objects to strings for serialization
        gpr_string_dict = {rxn_id: str(gpr) if gpr else None for rxn_id, gpr in gpr_dict.items()}
        
        # Generate deterministic seeds for parallel processing
        if random_seed is not None:
            seeds = [random_seed + i for i in range(n_permutations)]
        else:
            seeds = [np.random.randint(0, 1_000_000) for _ in range(n_permutations)]
        
        # Argument list for parallel processing
        arguments = [(genes, lfc_vector, task_structure, gpr_string_dict, seeds[i], or_func, "TIDE") \
                    for i in range(n_permutations)]
        
        # Parallel execution
        with Pool(processes=n_cpus) as pool:
            map_result = pool.map_async(MTEA_parallel_worker, arguments, chunksize=100)
            random_projection = np.array([array for array in map_result.get()])
    
    random_projection_df = pd.DataFrame(random_projection, columns=task_structure.index).T

    # Reduce from random reaction projections to random metabolic scores
    # (n_permutations * n_reactions) --> (n_permutations * n_tasks)
    random_scores = np.zeros((n_permutations, len(task_structure.columns)))
    for i, task in enumerate(task_structure.columns):
        rxns = [k for k in task_structure.index[task_structure[task]] if k in random_projection_df.index]
        random_scores[:,i] = random_projection_df.loc[rxns].mean().values

    return pd.DataFrame(random_scores, columns=task_structure.columns)


###########################################
# Context-aware mapping strategy
###########################################

def calculate_TIDE_scores_context_aware(gene_dict:dict, tasks_routes:dict, complex_cache:dict, or_func:str):
    """
    Context-aware analogue of `calculate_TIDE_scores`: instead of a single
    fixed per-task reaction set, scores every enumerated alternate route of
    each task (`mteapy.context_scoring.score_task`) and takes the
    best-supported one -- using "mean" aggregation for parity with TIDE's
    own published methodology (`calculate_TIDE_scores` also aggregates with
    `np.mean`; see `score_task`'s own docstring on why "mean", not this
    package's "min" default, matches TIDE).

    Parameters
    ----------
    gene_dict: dict
        A dictionary of genes to their expression values.

    tasks_routes: dict
        `{task_id: {route_id: reaction_id_set}}`, e.g. from
        `mteapy.routes.load_multiroute_tasks`/`load_task_routes`.

    complex_cache: dict
        `{reaction_id: candidate_complexes}`, from
        `mteapy.context_scoring.build_complex_cache`. Must cover every
        reaction appearing in `tasks_routes`.

    or_func: str ["absmax" | "max"]
        Function to evaluate OR rules within a GPR rule.

    Returns
    -------
    scores: numpy.ndarray
        An array of metabolic scores, in the same order as `tasks_routes`.
    """
    scores = [
        score_task(routes, complex_cache, gene_dict, aggregation="mean", task_id=task_id, or_func=or_func).score
        for task_id, routes in tasks_routes.items()
    ]
    return np.array(scores)


def calculate_random_TIDE_scores_context_aware(
        gene_dict:dict,
        tasks_routes:dict,
        complex_cache:dict,
        or_func:str,
        n_permutations:int = 1000,
        n_cpus:int = 1,
        random_seed:int = None
    ):
    """
    Context-aware analogue of `calculate_random_TIDE_scores`: computes Nº
    permutations * random metabolic scores (one per task in `tasks_routes`)
    to infer significancy, using `mteapy.context_scoring.score_task`
    instead of a fixed reaction set. See `calculate_TIDE_scores_context_aware`
    for `tasks_routes`/`complex_cache`; other parameters mirror
    `calculate_random_TIDE_scores`.

    Returns
    -------
    random_scores_df: pandas.DataFrame
        A pandas DataFrame where rows are permutations and columns are the
        task ids of `tasks_routes`.
    """
    genes = list(gene_dict.keys())
    lfc_vector = np.array(list(gene_dict.values()))
    task_ids = list(tasks_routes.keys())

    if random_seed is not None:
        np.random.seed(random_seed)

    if n_cpus <= 1:
        random_scores = np.zeros((n_permutations, len(task_ids)))
        for i in range(n_permutations):
            np.random.shuffle(lfc_vector)
            random_gene_dict = dict(zip(genes, lfc_vector))
            random_scores[i, :] = calculate_TIDE_scores_context_aware(
                random_gene_dict, tasks_routes, complex_cache, or_func
            )
    else:
        if random_seed is not None:
            seeds = [random_seed + i for i in range(n_permutations)]
        else:
            seeds = [np.random.randint(0, 1_000_000) for _ in range(n_permutations)]

        arguments = [
            (genes, lfc_vector, tasks_routes, complex_cache, seeds[i], or_func, "TIDE-context-aware")
            for i in range(n_permutations)
        ]
        with Pool(processes=n_cpus) as pool:
            map_result = pool.map_async(MTEA_parallel_worker, arguments, chunksize=100)
            random_scores = np.array([array for array in map_result.get()])

    return pd.DataFrame(random_scores, columns=task_ids)


def compute_TIDE(
        expr_data:pd.DataFrame,
        lfc_col:str,
        task_structure:pd.DataFrame,
        model:Model,
        or_func:str,
        n_permutations:int = 1000,
        n_cpus:int = 1,
        random_scores_flag:bool = False,
        random_seed:int = None,
        mapping_strategy:str = "classic",
        tasks_routes:dict = None,
        complex_cache:dict = None,
    ):
    """
    Wrapper function to compute the TIDE framework.

    Parameters
    ----------
    expr_data: pandas.DataFrame
        A pandas DataFrame containing gene expression values. Rows indices should correspond to the different genes, and columns should contain at least an expression column and a p-value column.
    
    lfc_col: str
        Column name in expr_data containing log-FC values.
    
    task_structure: pandas.DataFrame
        A boolean matrix where rows are reactions and columns metabolic tasks. Each cell contains ones or zeros, indicating whether a reaction is involved in a metabolic task.
    
    model: cobra.core.Model
        A COBRA metabolic model.
    
    or_func: str ["absmax" | "max"]
        Function to evaluate OR rules within a GPR rule.
    
    n_permutations: int 
        Number of permutations to inferr significancy (default: 1000).
    
    n_cpus: int 
        Number of jobs to parallelize the computation. If 1, the run is not parallellized (default: 1).
    
    random_scores_flag: bool
        Flag to indicate whether to add to the results DataFrame the null distributions (default: False).
    
    random_seed: int
        Random seed for reproducibility. If None, random seeds will be generated (default: None).

    mapping_strategy: str ["classic" | "context-aware"]
        "classic" (default) scores each task's single, fixed reaction set from
        `task_structure` (the traditional TIDE approach). "context-aware"
        instead scores every enumerated alternate route of each task and
        takes the best-supported one (`mteapy.context_scoring.score_task`),
        which needs `tasks_routes` and `complex_cache` instead of
        `task_structure`.

    tasks_routes: dict
        `{task_id: {route_id: reaction_id_set}}`. Required when
        mapping_strategy="context-aware"; ignored otherwise.

    complex_cache: dict
        `{reaction_id: candidate_complexes}`, from
        `mteapy.context_scoring.build_complex_cache`. Required when
        mapping_strategy="context-aware"; ignored otherwise.

    Returns
    -------
    TIDE_results: pandas.DataFrame
        A pandas DataFrame where rows are metabolic tasks and columns correspond to their score from TIDE, their average random score calculated from the null distribution, and their p-value.
    """
    gene_dict = dict(zip(expr_data.index, expr_data[lfc_col]))

    if mapping_strategy == "context-aware":
        if tasks_routes is None or complex_cache is None:
            raise ValueError("mapping_strategy='context-aware' requires both tasks_routes and complex_cache.")
        task_ids = list(tasks_routes.keys())
        scores = calculate_TIDE_scores_context_aware(gene_dict, tasks_routes, complex_cache, or_func)
        random_scores_df = calculate_random_TIDE_scores_context_aware(
            gene_dict, tasks_routes, complex_cache, or_func, n_permutations, n_cpus, random_seed
        )
    elif mapping_strategy == "classic":
        task_structure = task_structure.astype(bool)
        gpr_dict = {rxn.id: rxn.gpr for rxn in model.reactions if rxn.id in task_structure.index}
        task_ids = list(task_structure.columns)
        scores = calculate_TIDE_scores(gene_dict, task_structure, gpr_dict, or_func)
        random_scores_df = calculate_random_TIDE_scores(
            gene_dict, task_structure, gpr_dict, or_func, n_permutations, n_cpus, random_seed
        )
    else:
        raise ValueError(f"Unsupported mapping_strategy {mapping_strategy!r}. Please, use 'classic' or 'context-aware'.")

    pvalues = [calculate_pvalue(scores[i], random_scores_df[task]) for i, task in enumerate(task_ids)]

    TIDE_results = pd.DataFrame({"task_id": task_ids,
                                 "score": scores,
                                #"random_score": random_scores_df.mean(),
                                 "pvalue": pvalues
                                })
    
    if random_scores_flag:
        TIDE_results["random_score_array"] = np.array([";".join(random_scores_df.T.loc[task].astype(str)) \
                                                       for task in TIDE_results["task_id"]])

    return TIDE_results
