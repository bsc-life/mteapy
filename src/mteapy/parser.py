import argparse
from rich_argparse import ArgumentDefaultsRichHelpFormatter, RichHelpFormatter

import importlib.metadata

from mteapy import methods


RichHelpFormatter.group_name_formatter = str.upper
FORMATTER = ArgumentDefaultsRichHelpFormatter
FORMATTER.group_name_formatter = str.upper
FORMATTER.styles["argparse.prog"] = "bold"

VERSION = importlib.metadata.version("mteapy")


def _add_mapping_strategy_args(parser):
    """Shared by TIDE and CellFie: which reaction set represents a task.

    "classic" (default) is each method's traditional single, fixed
    reaction set. "context-aware" instead enumerates alternate routes for
    each task and scores every one of them (mteapy.context_scoring),
    taking the best-supported route per sample -- an orthogonal axis to
    which scoring convention (TIDE vs CellFie) is used, not a third
    framework of its own; see mteapy.tide.compute_TIDE/mteapy.cellfie
    .compute_CellFie's own `mapping_strategy` docs for what each method
    actually does with it.
    """
    parser.add_argument("--mapping-strategy", action="store", type=str, dest="mapping_strategy",
                         choices=["classic", "context-aware"], default="classic",
                         help="'classic': score each task's single, fixed reaction set. 'context-aware': "
                              "enumerate alternate routes per task and score the best-supported one "
                              "(requires routes already built with `run-mtea tasks enumerate-routes`).")
    parser.add_argument("--model", action="store", type=str, dest="model", default=None,
                         help="Which registered model to use: a key from `run-mtea models list` or a model folder "
                              "(default: the bundled Human-GEM). Supplies the routes database and model file. "
                              "Only used with --mapping-strategy context-aware.")
    parser.add_argument("--routes-db", action="store", type=str, dest="routes_db", default=None,
                         help="Override --model's routes database with this file. "
                              "Only used with --mapping-strategy context-aware.")
    parser.add_argument("--task-list", action="store", type=str, dest="task_list", default="HumanGEM-Full",
                         help="Which stored task list's enumerated routes to use (see `run-mtea tasks import "
                              "--task-list`). Only used with --mapping-strategy context-aware.")
    parser.add_argument("--routes-model-file", action="store", type=str, dest="routes_model_file", default=None,
                         help="Override --model's model file with this one (it must match the routes database's "
                              "reaction ids). Only used with --mapping-strategy context-aware -- "
                              "this is independent of, and does not need to match, the model implied by "
                              "--mapping-strategy classic's bundled task_structure_matrix.tsv.")


def mtea_parser():

    ###########################################
    # General parser
    ###########################################

    parser = argparse.ArgumentParser(
        description="""
        Command line tool to perform Metabolic Task Enrichment Analysis (MTEA) using several methods.\n
        All analysis uses the Human1 metabolic model (Robinson et al., 2020).
        Uses metabolic task list from Richelle et. al 2021 and curated to Human1 (see -t for more info on metabolic tasks).
        """,
        formatter_class=RichHelpFormatter
    )

    parser.add_argument("-v", "--version", action="version", version=f"[argparse.prog]%(prog)s[/] {VERSION}")

    parser.add_argument("-c", "--cite", action="store_true", dest="citation_flag", help="prints information regarding citation of methods")

    subparser = parser.add_subparsers(title="commands", required=False, dest="command")

    ###########################################
    # "tasks" group -- task/route database preparation
    ###########################################

    models_parser = subparser.add_parser("models", help="list the registered models and check their files", formatter_class=FORMATTER)
    models_subparser = models_parser.add_subparsers(title="models commands", required=False, dest="models_command")
    models_subparser.add_parser("list", help="list every model found (bundled, ~/.mteapy/models, or $MTEAPY_MODELS) "
                                             "with its task lists and whether its files check out",
                                formatter_class=FORMATTER)

    tasks_parser = subparser.add_parser("tasks", help="prepare/manage the task-route database", formatter_class=FORMATTER)
    tasks_subparser = tasks_parser.add_subparsers(title="tasks commands", required=False, dest="tasks_command")

    import_parser = tasks_subparser.add_parser(
        "import",
        help="store a RAVEN-style task list file, validated against a model, as a named task list in the "
             "task/route database (creating the database if needed)",
        formatter_class=FORMATTER,
    )
    import_parser.add_argument("task_file", action="store", help="Path to a RAVEN-style metabolic task list file.")
    import_parser.add_argument("--task-list", action="store", type=str, dest="task_list", required=True,
                                help="Name to store the list under (e.g. 'HumanGEM-Full'). Task lists are immutable: "
                                     "to change a task, import the edited file under a new name.")
    import_parser.add_argument("--description", action="store", type=str, dest="description", default=None)
    import_parser.add_argument("--model", action="store", type=str, dest="model", default=None,
                                help="Registered model (key from `run-mtea models list`, or a model folder) whose "
                                     "model file and database to use (default: the bundled Human-GEM).")
    import_parser.add_argument("--model-file", action="store", type=str, dest="model_file", default=None,
                                help="Override the model file the tasks are validated against.")
    import_parser.add_argument("--db", action="store", type=str, dest="db", default=None,
                                help="Override the task/route database to write into (created if missing). Giving "
                                     "both --db and --model-file builds a new database without any registered model.")
    import_parser.add_argument("--model-name", action="store", type=str, dest="model_name", default=None,
                                help="Name to register the model under; required the first time a model is used with a database.")
    import_parser.add_argument("--model-version", action="store", type=str, dest="model_version", default=None)

    enum_parser = tasks_subparser.add_parser(
        "enumerate-routes",
        help="enumerate alternate reaction routes (with fluxes) for every valid task of an imported task list "
             "(dataset-agnostic; run once per model+task-list, reused by every sample later scored with "
             "--mapping-strategy context-aware)",
        formatter_class=FORMATTER,
    )
    enum_parser.add_argument("--task-list", action="store", type=str, dest="task_list", required=True,
                              help="Name of a task list already stored with `run-mtea tasks import`.")
    enum_parser.add_argument("--model", action="store", type=str, dest="model", default=None,
                              help="Registered model (key from `run-mtea models list`, or a model folder) whose "
                                   "model file and database to use (default: the bundled Human-GEM).")
    enum_parser.add_argument("--model-file", action="store", type=str, dest="model_file", default=None,
                              help="Override the SBML model file to enumerate against; must be the model the "
                                   "database was built with.")
    enum_parser.add_argument("--db", action="store", type=str, dest="db", default=None,
                              help="Override the task/route database to write into.")
    enum_parser.add_argument("--solver", action="store", type=str, dest="solver", default=None,
                              help="COBRApy solver name (e.g. 'gurobi', 'cplex', 'glpk'). Default: whatever cobra picks.")
    enum_parser.add_argument("--max-routes", action="store", type=int, dest="max_routes", default=10,
                              help="Maximum number of alternate routes to search for per task.")
    enum_parser.add_argument("--mode", action="store", type=str, dest="mode", choices=["reset", "resume"], default=None,
                              help="'reset': permanently delete existing routes for this task list first, "
                                   "after confirmation (see --yes). 'resume': skip already-exhaustive tasks, "
                                   "continue capped ones from their stored routes. Default: full re-run of every task.")
    enum_parser.add_argument("--yes", action="store_true", dest="skip_confirmation",
                              help="Don't prompt for confirmation before --mode reset permanently deletes existing routes.")
    enum_parser.add_argument("--only", action="store", type=str, dest="only", default=None,
                              help="Comma-separated task ids to process, ignoring --limit/--start-at (for targeted repairs).")
    enum_parser.add_argument("--limit", action="store", type=int, dest="limit", default=None,
                              help="Only process the first N tasks (for timing/partial runs).")
    enum_parser.add_argument("--start-at", action="store", type=int, dest="start_at", default=0,
                              help="Skip tasks with a numeric id < this.")

    ###########################################
    # "analyze" group -- sample expression analysis
    ###########################################

    analyze_parser = subparser.add_parser("analyze", help="score sample gene expression data against metabolic tasks", formatter_class=FORMATTER)
    analyze_subparser = analyze_parser.add_subparsers(title="analyze commands", required=False, dest="analyze_command")

    ###########################################
    # TIDE parser
    ###########################################

    TIDE_parser = analyze_subparser.add_parser("TIDE", help="performs TIDE analysis (Doughberty et al., 2021) and complementary TIDE-essential", formatter_class=FORMATTER)

    TIDE_parser.add_argument("dea_file", action="store", help="Filename for a differential expression analysis results file. It should contain at least three columns: genic (string), log-FC (numeric) and significance (numeric, e.g.: p-value, adjusted p-value, FDR). Genes must be stored as EnsemblIDs.")

    TIDE_parser.add_argument("-d", "--delim", action="store", type=str, dest="sep", default="\t", help="Field delimiter for inputed file.")

    TIDE_parser.add_argument("-o", "--out", action="store", type=str, dest="out_filename", default="tide_results.tsv", help="Name (and location) to store the analysis’ results. They will be stored in a tab-sepparated file, so filenames should contain the .tsv or .txt extensions.")

    TIDE_parser.add_argument("--gene_col", action="store", type=str, dest="gene_col", default="gene_id", help="Name of the column in the inputed file containing gene names/symbols. Genes must be stored as EnsemblIDs.")

    TIDE_parser.add_argument("--lfc_col", action="store", type=str, dest="lfc_col", default="log2FoldChange", help="Name of the column in the inputed file containing log-FC values.")

    TIDE_parser.add_argument("--pvalue_col", action="store", type=str, dest="pvalue_col", default="pvalue", help="Name of the column in the inputed file containing significance values. Only required if the flag --mask_lfc_values is True.")

    TIDE_parser.add_argument("-a", "--alpha", action="store", type=float, dest="alpha", default=0.05, help="Significance threshold to mask log-FC. Only required if the flag --mask_lfc_values is True.")

    TIDE_parser.add_argument("-n", "--n_permutations", action="store", type=int, dest="n_permutations", default=1000, help="Number of permutations to infer p-values for the metabolic scores. The resolution of the computed p-values will depend on this number.")

    TIDE_parser.add_argument("--or_func", action="store", type=str, dest="or_func", choices=["max", "absmax"], default="absmax", help="Name of the function that will be used to resolve OR relationships in gene-protein-reaction (GPR) rules. Possible values are absmax, which will return the absolute maximum value, and max, which will return the maximum value.")

    TIDE_parser.add_argument("--n_cpus", action="store", type=int, dest="n_cpus", default=1,help="Number of CPUs for parallel execution.")

    TIDE_parser.add_argument("--mask_lfc_values", action="store_true", dest="filter_lfc", help="Flag to indicate whether to mask log-FC values to 0 according to their significance. That is, if a log-FC value is non-significant (determined by the user), they will be masked to 0.")

    TIDE_parser.add_argument("--random_scores", action="store_true", dest="random_scores_flag", help="Flag to indicate whether to return the null distribution of random scores used to inferr significance with the results file.")

    TIDE_parser.add_argument("--random_seed", action="store", type=int, dest="random_seed", default=None, help="Random seed for reproducibility. If not provided, a random seed will be generated for each run.")

    _add_mapping_strategy_args(TIDE_parser)

    TIDE_parser.add_argument("--permutation-strategy", action="store", type=str, dest="permutation_strategy",
                              choices=["argmax", "fixed-route"], default="argmax",
                              help="Only used with --mapping-strategy context-aware (CellFie has no permutation "
                                   "test to apply this to). 'argmax' (default, more rigorous): re-runs the full "
                                   "best-of-routes search for every permutation, exactly mirroring the real "
                                   "score's procedure. 'fixed-route' (faster): fixes each task's real-data "
                                   "winning route once and only scores that single reaction set under every "
                                   "permutation -- cheaper, but the resulting null distribution runs slightly "
                                   "less conservative than 'argmax', since it never re-earns \"best of K "
                                   "routes\" under permutation.")

    ###########################################
    # CellFie parser
    ###########################################

    CellFie_parser = analyze_subparser.add_parser("CellFie", help="performs CellFie analysis (Richelle et al., 2021)", formatter_class=FORMATTER)

    CellFie_parser.add_argument("expr_file", action="store", help="Filename for a normalized gene expression file (e.g., TPM). It should contain at least one column with gene names/symbols. Genes must be stored as EnsemblIDs.")

    CellFie_parser.add_argument("-d", "--delim", action="store", type=str, dest="sep", default="\t", help="	Field delimiter for inputed file.")

    CellFie_parser.add_argument("-o", "--out", action="store", type=str, dest="out_dir", default="cellfie_results", help="Directory to store the analysis' results. The result file(s) will be stored in the specified directory in a tab-sepparated format (.tsv).")

    CellFie_parser.add_argument("--gene_col", action="store", type=str, dest="gene_col", default="geneID", help="Name of the column in the inputed file containing gene names/symbols. Genes must be stored as EnsemblIDs.")

    methods.add_arguments(CellFie_parser, methods.CELLFIE)

    CellFie_parser.add_argument("--binary_scores", action="store_true", dest="binary_scores_flag", help="Flag to indicate whether to also return the binary metabolic score matrix as a second result file. See the original publication for more details.")

    _add_mapping_strategy_args(CellFie_parser)

    ###########################################
    # TAS parser
    ###########################################

    TAS_parser = analyze_subparser.add_parser("TAS", help="performs Task Activity Score (TAS) analysis: plain expression-through-GPR projection, no threshold transform", formatter_class=FORMATTER)

    TAS_parser.add_argument("expr_file", action="store", help="Filename for a normalized gene expression file (e.g., TPM). It should contain at least one column with gene names/symbols. Genes must be stored as EnsemblIDs.")

    TAS_parser.add_argument("-d", "--delim", action="store", type=str, dest="sep", default="\t", help="Field delimiter for inputed file.")

    TAS_parser.add_argument("-o", "--out", action="store", type=str, dest="out_dir", default="tas_results", help="Directory to store the analysis' results. The result file(s) will be stored in the specified directory in a tab-sepparated format (.tsv).")

    TAS_parser.add_argument("--gene_col", action="store", type=str, dest="gene_col", default="geneID", help="Name of the column in the inputed file containing gene names/symbols. Genes must be stored as EnsemblIDs.")

    methods.add_arguments(TAS_parser, methods.TAS)

    _add_mapping_strategy_args(TAS_parser)

    return parser
