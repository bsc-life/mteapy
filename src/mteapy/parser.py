import argparse
from rich_argparse import ArgumentDefaultsRichHelpFormatter, RichHelpFormatter

import importlib.metadata


RichHelpFormatter.group_name_formatter = str.upper
FORMATTER = ArgumentDefaultsRichHelpFormatter
FORMATTER.group_name_formatter = str.upper
FORMATTER.styles["argparse.prog"] = "bold"

VERSION = importlib.metadata.version("mteapy")


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
    # TIDE-essential parser
    ###########################################
    
    TIDE_essential_parser = subparser.add_parser("TIDE-essential", help="performs TIDE analysis using only essential genes", formatter_class=FORMATTER)
    
    TIDE_essential_parser.add_argument("dea_file", action="store", help="Filename for a differential expression analysis results file. It should contain at least three columns: genic (string), log-FC (numeric) and significance (numeric, e.g.: p-value, adjusted p-value, FDR). Genes must be stored as EnsemblIDs.")
    
    TIDE_essential_parser.add_argument("-d", "--delim", action="store", type=str, dest="sep", default="\t", help="Field delimiter for inputed file.")
    
    TIDE_essential_parser.add_argument("-o", "--out", action="store", type=str, dest="out_filename",  default="tide_e_results.tsv", help="Name (and location) to store the analysis' results. They will be stored in a tab-sepparated file, so filenames should contain the .tsv or .txt extensions.")
    
    TIDE_essential_parser.add_argument("--gene_col", action="store", type=str, dest="gene_col", default="geneID", help="Name of the column in the inputed file containing gene names/symbols. Genes must be stored as EnsemblIDs.")

    TIDE_essential_parser.add_argument("--lfc_col", action="store", type=str, dest="lfc_col", default="log2FoldChange", help="Name of the column in the inputed file containing log-FC values.")

    TIDE_essential_parser.add_argument("--pvalue_col", action="store", type=str, dest="pvalue_col", default="pvalue", help="Name of the column in the inputed file containing significance values. Only required if the flag --mask_lfc_values is True.")
    
    TIDE_essential_parser.add_argument("-a", "--alpha", action="store", type=float, dest="alpha", default=0.05, help="Significance threshold to mask log-FC. Only required if the flag --mask_lfc_values is True.")

    TIDE_essential_parser.add_argument("-n", "--n_permutations", action="store", type=int, dest="n_permutations", default=1000, help="Number of permutations to infer p-values for the metabolic scores. The resolution of the computed p-values will depend on this number.")
    
    TIDE_essential_parser.add_argument("--n_cpus", action="store", type=int, dest="n_cpus", default=1, help="Number of CPUs for parallel execution.")
    
    TIDE_essential_parser.add_argument("--mask_lfc_values", action="store_true", dest="filter_lfc", help="Flag to indicate whether to mask log-FC values to 0 according to their significance. That is, if a log-FC value is non-significant (determined by the user), they will be masked to 0.")
    
    TIDE_essential_parser.add_argument("--random_scores", action="store_true", dest="random_scores_flag", help="Flag to indicate whether to return the null distribution of random scores used to inferr significance with the results file.")


    ###########################################
    # TIDE parser
    ###########################################

    TIDE_parser = subparser.add_parser("TIDE", help="performs TIDE analysis (Doughberty et al., 2021)", formatter_class=FORMATTER)
    
    TIDE_parser.add_argument("dea_file", action="store", help="Filename for a differential expression analysis results file. It should contain at least three columns: genic (string), log-FC (numeric) and significance (numeric, e.g.: p-value, adjusted p-value, FDR). Genes must be stored as EnsemblIDs.") 
    
    TIDE_parser.add_argument("-s", "--secretory", action="store_true", dest="secretory_flag", help="whether to also use the secretory tasks expansion")
    
    TIDE_parser.add_argument("-d", "--delim", action="store", type=str, dest="sep", default="\t", help="Field delimiter for inputed file.")
    
    TIDE_parser.add_argument("-o", "--out", action="store", type=str, dest="out_filename", default="tide_results.tsv", help="Name (and location) to store the analysis’ results. They will be stored in a tab-sepparated file, so filenames should contain the .tsv or .txt extensions.")

    TIDE_parser.add_argument("--gene_col", action="store", type=str, dest="gene_col", default="geneID", help="Name of the column in the inputed file containing gene names/symbols. Genes must be stored as EnsemblIDs.")
    
    TIDE_parser.add_argument("--lfc_col", action="store", type=str, dest="lfc_col", default="log2FoldChange", help="Name of the column in the inputed file containing log-FC values.")

    TIDE_parser.add_argument("--pvalue_col", action="store", type=str, dest="pvalue_col", default="pvalue", help="Name of the column in the inputed file containing significance values. Only required if the flag --mask_lfc_values is True.")
    
    TIDE_parser.add_argument("-a", "--alpha", action="store", type=float, dest="alpha", default=0.05, help="Significance threshold to mask log-FC. Only required if the flag --mask_lfc_values is True.")

    TIDE_parser.add_argument("-n", "--n_permutations", action="store", type=int, dest="n_permutations", default=1000, help="Number of permutations to infer p-values for the metabolic scores. The resolution of the computed p-values will depend on this number.")
    
    TIDE_parser.add_argument("--or_func", action="store", type=str, dest="or_func", choices=["max", "absmax"], default="absmax", help="Name of the function that will be used to resolve OR relationships in gene-protein-reaction (GPR) rules. Possible values are absmax, which will return the absolute maximum value, and max, which will return the maximum value.")

    TIDE_parser.add_argument("--n_cpus", action="store", type=int, dest="n_cpus", default=1,help="Number of CPUs for parallel execution.")

    TIDE_parser.add_argument("--mask_lfc_values", action="store_true", dest="filter_lfc", help="Flag to indicate whether to mask log-FC values to 0 according to their significance. That is, if a log-FC value is non-significant (determined by the user), they will be masked to 0.")
    
    TIDE_parser.add_argument("--random_scores", action="store_true", dest="random_scores_flag", help="Flag to indicate whether to return the null distribution of random scores used to inferr significance with the results file.")
    
    
    ###########################################
    # CellFie parser
    ###########################################
    
    CellFie_parser = subparser.add_parser("CellFie", help="performs CellFie analysis (Richelle et al., 2021)", formatter_class=FORMATTER)
    
    CellFie_parser.add_argument("expr_file", action="store", help="Filename for a normalized gene expression file (e.g., TPM). It should contain at least one column with gene names/symbols. Genes must be stored as EnsemblIDs.")

    # CellFie_parser.add_argument("-s", "--secretory", action="store_true", dest="secretory_flag", help="whether to use the secretory tasks expansion")
    
    CellFie_parser.add_argument("-d", "--delim", action="store", type=str, dest="sep", default="\t", help="	Field delimiter for inputed file.")
    
    CellFie_parser.add_argument("-o", "--out", action="store", type=str, dest="out_dir", default="Directory to store the analysis' results. The result file(s) will be stored in the specified directory in a tab-sepparated format (.tsv).")

    CellFie_parser.add_argument("--gene_col", action="store", type=str, dest="gene_col", default="geneID", help="Name of the column in the inputed file containing gene names/symbols. Genes must be stored as EnsemblIDs.")
    
    CellFie_parser.add_argument("--threshold_type", action="store", type=str, dest="thresh_type", default="local", choices=["local","global"], help="Determines the threshold approach to be used. A global approach used the same threshold for all genes whereas a local approach uses a different threshold for each gene when computing the gene activity levels.")
    
    CellFie_parser.add_argument("--global_threshold_type", action="store", type=str, dest="global_thresh_type", default="percentile", choices=["value","percentile"], help="Whether to use a value or a percentile of the distribution of all genes as global treshold for all genes.")

    CellFie_parser.add_argument("--global_value", action="store", type=float, dest="global_value", default=0.75, help="Value to use as global threshold according to the global_threshold_type option selected. Note that percentile values must be between 0 and 1.")
    
    CellFie_parser.add_argument("--local_threshold_type", action="store", type=str, dest="local_thresh_type", default="minmaxmean", choices=["minmaxmean","mean"], help="Determines the threshold type to be used in a local approach. minmaxmean: the threshold for each gene is determined by the mean of expression values across all conditions/samples but must be higher or equal than a lower bound and lower or equal to an upper bound. mean: the threshold of a gene is determined as its mean expression across all conditions/samples.")

    CellFie_parser.add_argument("-s", "--secretory", action="store_true", dest="secretory_flag", help="whether to also use the secretory tasks expansion")
    
    CellFie_parser.add_argument("--minmaxmean_threshold_type", action="store", type=str, dest="minmaxmean_thresh_type", default="percentile", choices=["percentile","value"], help="Whether to use value or percentile of the distribution of all genes as upper and lower bounds.")
    
    CellFie_parser.add_argument("--upper_bound", action="store", type=float, dest="upper_bound", default=0.75, help="Upper bound value to be used according to the minmaxmean_threshold_type. Note that percentile values must be between 0 and 1")
    
    CellFie_parser.add_argument("--lower_bound", action="store", type=float, dest="lower_bound", default=0.25, help="Lower bound value to be used according to the minmaxmean_threshold_type. Note that percentile values must be between 0 and 1.") 

    CellFie_parser.add_argument("--binary_scores", action="store_true", dest="binary_scores_flag", help="Flag to indicate whether to also return the binary metabolic score matrix as a second result file. See the original publication for more details.")
    
    return parser