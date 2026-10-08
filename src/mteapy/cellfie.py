import pandas as pd
import numpy as np

from cobra.core.gene import GPR
from cobra.core import Model

from mteapy.context_scoring import score_tasks_matrix
from mteapy.utils import map_gpr_w_names


###########################################
# CellFie functions
###########################################

def calculate_GAL(
        expr_data:pd.DataFrame,
        thresh_type:str = "local",
        local_thresh_type:str = "minmaxmean",
        minmaxmean_thresh_type:str = "percentile",
        upper_bound:float = 0.75, lower_bound:float = 0.25,
        global_thresh_type:str = None, global_value:float = None,
        log_transformed:bool = False,
    ):
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

    log_transformed: bool
        Whether `expr_data` is already log-transformed (default: False --
        raw linear expression, e.g. TPM/FPKM). The original MATLAB CellFie
        (`CellFie.m`) computes percentile-based thresholds by taking the
        percentile in log10 space and exponentiating back
        (`10^percentile(log10(x))`), not a percentile taken directly on x --
        those give different values for continuous data, since percentile
        *rank* is preserved under a monotonic transform but linear
        interpolation between order statistics is not. When
        `log_transformed=False` (default), this replicates that exact
        MATLAB computation, on the assumption `expr_data` is raw, matching
        what the rest of CellFie's formula (`5*log(1+value/threshold)`)
        itself assumes throughout. Setting `log_transformed=True` is a
        deliberate opt-out for a caller who knows their data is already
        log-transformed and wants percentiles computed directly on it as
        given, rather than treating it as MATLAB's algorithm does -- it
        does not change the final `5*log(1+value/threshold)` step, so a
        caller choosing this should be sure that's still the right formula
        for their data.

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

    def _percentile_threshold(p):
        if log_transformed:
            return np.quantile(linear_data, p)
        return 10 ** np.quantile(np.log10(linear_data), p)

    if thresh_type not in ("local", "global"):
        raise ValueError(f"Unsupported thresh_type ({thresh_type}). Please, use 'local' or 'global'.")

    # Local approach
    if thresh_type == "local":
        if local_thresh_type == "mean":
            thresholds = np.mean(expr_data_array, axis=1)

        elif local_thresh_type == "minmaxmean":
            if minmaxmean_thresh_type == "value":
                local_upper = upper_bound
                local_lower = lower_bound
            elif minmaxmean_thresh_type == "percentile":
                local_upper = _percentile_threshold(upper_bound)
                local_lower = _percentile_threshold(lower_bound)

            expr_mean = np.mean(expr_data_array, axis=1)
            thresholds = np.clip(expr_mean, local_lower, local_upper)

        gene_scores = 5 * np.log(1 + expr_data_array / thresholds[:,np.newaxis])

    # Global approach
    elif thresh_type == "global":
        if global_thresh_type == "value":
            global_threshold = global_value

        elif global_thresh_type == "percentile":
            global_threshold = _percentile_threshold(global_value)

        gene_scores = 5 * np.log(1 + expr_data_array / global_threshold)

    return pd.DataFrame(gene_scores, index=expr_data.index, columns=expr_data.columns)


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
        A cell is NaN where the reaction's GPR has no gene with data for
        that sample (matching the original MATLAB CellFie's `-1` "no data"
        sentinel -- see `mteapy.utils.map_gpr_w_names`) -- callers must
        exclude these from aggregation rather than treat them as a real
        zero score.

    Notes
    -----
    RALs are computed as indicated in the CellFie paper.
    """
    all_reactions = list(gpr_dict.keys())
    all_genes = gal_df.index

    ral = np.full((len(all_reactions), len(gal_df.columns)), np.nan)

    for i, sample in enumerate(gal_df.columns):
        gene_dict = dict(zip(all_genes, gal_df[sample].to_numpy()))
        rxn_projection = [map_gpr_w_names(gpr_dict[rxn], gene_dict) for rxn in all_reactions]

        # Significance (1 / how often a gene wins a reaction across the
        # whole model, for this sample) is computed only over reactions
        # that actually have data -- a "no data" reaction never has a
        # winning gene, so it shouldn't dilute another gene's count.
        genes_with_data = [gene for score, gene in rxn_projection if score is not None]
        if genes_with_data:
            genes_unique, counts = np.unique(genes_with_data, return_counts=True)
            significance_dict = {gene: 1 / k for gene, k in zip(genes_unique, counts)}
        else:
            significance_dict = {}

        for j, (score, gene) in enumerate(rxn_projection):
            if score is not None:
                ral[j, i] = score * significance_dict[gene]
            # else leave as NaN -- "no data" for this reaction/sample.

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
        A cell is NaN if the task has no reaction with any data for that
        sample -- this package's equivalent of the original MATLAB
        CellFie's `-1` "no data" sentinel for a whole task/sample -- rather
        than a real, computed score of zero.

    binary_scores_df: pandas.DataFrame
        A pandas DataFrame containing the binary metabolic scores after applying a certain threshold of activity. Rows correspond to metabolic tasks and columns to samples.
        NaN wherever `metabolic_scores_df` is NaN, for the same reason.

    Notes
    -----
    Metabolic scores and binary metabolic scores are computed as indicated in the CellFie paper.
    A reaction absent from `ral_df` (no GPR at all) or NaN in it (GPR
    present but no measured genes for that sample -- see `calculate_RAL`)
    is excluded from a task's average rather than treated as a zero score,
    matching the original MATLAB CellFie's behavior of dropping any
    reaction that has no data for every sample from that task's
    essential-reaction list, instead of scoring it in as a confirmed zero.
    """
    task_to_rxns = {task: task_structure.index[task_structure[task]].to_list() for task in task_structure.columns}
    metabolic_scores = np.full((len(task_to_rxns), len(ral_df.columns)), np.nan)

    for i, sample in enumerate(ral_df.columns):
        rxn_dict = dict(zip(ral_df.index, ral_df[sample]))
        for k, task in enumerate(task_to_rxns):
            values = [
                rxn_dict[rxn] for rxn in task_to_rxns[task]
                if rxn in rxn_dict and not np.isnan(rxn_dict[rxn])
            ]
            if values:
                metabolic_scores[k, i] = np.mean(values)

    metabolic_scores_df = pd.DataFrame(
        metabolic_scores,
        index=pd.Index(task_to_rxns.keys(), name="task_id"),
        columns=ral_df.columns
    )
    is_active = metabolic_scores_df >= 5 * np.log(2)
    binary_scores_df = is_active.astype(float)
    binary_scores_df[metabolic_scores_df.isna()] = np.nan

    return metabolic_scores_df, binary_scores_df


###########################################
# Context-aware mapping strategy
###########################################

def calculate_CellFie_scores_context_aware(gal_df:pd.DataFrame, tasks_routes:dict, model:Model):
    """
    Context-aware analogue of `calculate_CellFie_scores`: instead of a
    single fixed per-task reaction set, scores every enumerated alternate
    route of each task (`mteapy.context_scoring.score_task`, via
    `score_tasks_matrix`) and takes the best-supported one for each sample
    -- using "mean" aggregation and `or_func="max"` to match CellFie's own
    non-negative-expression convention (`calculate_CellFie_scores` also
    aggregates with `np.mean`).

    Parameters
    ----------
    gal_df: pandas.DataFrame
        A pandas DataFrame containing the GALs (genes x samples) -- the
        same input `calculate_RAL` takes in the classic pipeline; this
        skips straight from gene-level GALs to task scores, since a route's
        reaction set isn't fixed and so can't be pre-reduced to a flat RAL
        matrix the way `calculate_RAL` does.

    tasks_routes: dict
        `{task_id: {route_id: reaction_id_set}}`, e.g. from
        `mteapy.taskdb.load_task_list_routes`.

    model: cobra.core.Model
        The COBRA model the routes' reaction ids come from, used to look up
        each reaction's GPR for `build_complex_cache`.

    Returns
    -------
    metabolic_scores_df: pandas.DataFrame
        A pandas DataFrame containing the metabolic scores where rows correspond to metabolic tasks and columns to samples.

    binary_scores_df: pandas.DataFrame
        A pandas DataFrame containing the binary metabolic scores after applying a certain threshold of activity. Rows correspond to metabolic tasks and columns to samples.
    """
    metabolic_scores_df, _ = score_tasks_matrix(tasks_routes, model, gal_df, aggregation="mean", or_func="max")
    metabolic_scores_df.index.name = "task_id"
    binary_scores_df = (metabolic_scores_df >= 5 * np.log(2)).astype(int)

    return metabolic_scores_df, binary_scores_df


def compute_CellFie(
        expr_data:pd.DataFrame,
        task_structure:pd.DataFrame,
        model:Model,
        thresh_type:str = "local",
        local_thresh_type:str = "minmaxmean",
        minmaxmean_thresh_type:str = "percentile",
        upper_bound:float = 0.75, lower_bound:float = 0.25,
        global_thresh_type:str = None, global_value:float = None,
        log_transformed:bool = False,
        mapping_strategy:str = "classic",
        tasks_routes:dict = None,
    ):
    """
    Wrapper function to compute the CellFie framework.

    Parameters
    ----------
    expr_data: pandas.DataFrame
        A pandas DataFrame containing gene expression values where rows are genes and columns are samples. Genes should be inputed as the index of the DataFrame.

    task_structure: pandas.DataFrame
        A boolean matrix where rows are reactions and columns metabolic tasks. Each cell contains ones or zeros, indicating whether a reaction is involved in a metabolic task.
    
    model: cobra.core.Model
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

    mapping_strategy: str ["classic" | "context-aware"]
        "classic" (default) scores each task's single, fixed reaction set from
        `task_structure` (the traditional CellFie approach). "context-aware"
        instead scores every enumerated alternate route of each task and
        takes the best-supported one for each sample
        (`mteapy.context_scoring.score_tasks_matrix`), which needs
        `tasks_routes` instead of `task_structure`.

    tasks_routes: dict
        `{task_id: {route_id: reaction_id_set}}`. Required when
        mapping_strategy="context-aware"; ignored otherwise.

    log_transformed: bool
        Passed through to `calculate_GAL` -- see its own docstring. Default
        False assumes `expr_data` is raw (linear) expression, matching what
        the original CellFie algorithm assumes throughout.

    Returns
    -------
    metabolic_scores_df: pandas.DataFrame
        A pandas DataFrame containing the metabolic scores where rows correspond to metabolic tasks and columns to samples.

    binary_scores_df: pandas.DataFrame
        A pandas DataFrame containing the binary metabolic scores after applying a certain threshold of activity. Rows correspond to metabolic tasks and columns to samples.
    """
    gal_df = calculate_GAL(
        expr_data,
        thresh_type,
        local_thresh_type,
        minmaxmean_thresh_type,
        upper_bound, lower_bound,
        global_thresh_type, global_value,
        log_transformed,
    )

    if mapping_strategy == "context-aware":
        if tasks_routes is None:
            raise ValueError("mapping_strategy='context-aware' requires tasks_routes.")
        metabolic_scores_df, binary_scores_df = calculate_CellFie_scores_context_aware(gal_df, tasks_routes, model)
    elif mapping_strategy == "classic":
        task_structure = task_structure.astype(bool)
        # CellFie only uses reactions with valid GPR rules, internally they are considered as 0.0
        gpr_dict = {rxn.id: rxn.gpr for rxn in model.reactions if rxn.gpr != GPR.from_string("")}
        ral_df = calculate_RAL(gal_df, gpr_dict)
        metabolic_scores_df, binary_scores_df = calculate_CellFie_scores(ral_df, task_structure)
    else:
        raise ValueError(f"Unsupported mapping_strategy {mapping_strategy!r}. Please, use 'classic' or 'context-aware'.")

    return metabolic_scores_df, binary_scores_df