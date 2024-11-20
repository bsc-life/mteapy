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
        All analysis are based on the Human1 metabolic model (Robinson et al., 2020). 
        Uses metabolic task list from Richelle et. al 2021 and curated to Human1 (see -t for more info on metabolic tasks).
        """,
        formatter_class=RichHelpFormatter
    )

    parser.add_argument("-v", "--version", action="version", version=f"[argparse.prog]%(prog)s[/] {VERSION}")
    
    parser.add_argument("-c", "--cite", action="store_true", dest="citation_flag", help="prints information regarding citation of methods")

    # parser.add_argument("-t", "--task_metadata", action="store_true", dest="task_metadata", help="saves metabolic task metadata into current directory")

    # parser.add_argument("-s", "--task_metadata_sec", action="store_true", dest="task_metadata_sec", help="saves secretory metabolic task metadata into current directory")

    subparser = parser.add_subparsers(title="commands", required=False, dest="command")
    

    ###########################################
    # TIDE-essential parser
    ###########################################
    
    TIDE_essential_parser = subparser.add_parser("TIDE-essential", help="performs TIDE analysis using only essential genes", formatter_class=FORMATTER)
    
    TIDE_essential_parser.add_argument("dea_file", action="store", help="file containing DEA results")
    
    TIDE_essential_parser.add_argument("-d", "--delim", action="store", type=str, dest="sep", default="\t", help="files field separator")
    
    TIDE_essential_parser.add_argument("-o", "--out", action="store", type=str, dest="out_filename",  default="tide_e_results.tsv", help="out file name to store permutation results")
    
    TIDE_essential_parser.add_argument("-n", "--n_permutations", action="store", type=int, dest="n_permutations", default=1000, help="number of permutations to generate")
    
    TIDE_essential_parser.add_argument("--n_cpus", action="store", type=int, dest="n_cpus", default=1, help="number of CPUs for parallel processing")
    
    TIDE_essential_parser.add_argument("--lfc_col", action="store", type=str, dest="lfc_col", default="log2FoldChange", help="column name containing LFC values")
    
    TIDE_essential_parser.add_argument("--gene_col", action="store", type=str, dest="gene_col", default="geneID", help="column name containing gene symbols")
    
    TIDE_essential_parser.add_argument("--mask_lfc_values", action="store_true", dest="filter_lfc", help="flag to indicate whether to mask to 0 non-significant LFC values using their p-values")
    
    TIDE_essential_parser.add_argument("--random_scores", action="store_true", dest="random_scores_flag", help="flag to indicate whether to return the random score distributions with the results")
    
    TIDE_essential_parser.add_argument("--pvalue_col", action="store", type=str, dest="pvalue_col", default="pvalue", help="column name containing p-values (only used when masking LFC values)")
    
    TIDE_essential_parser.add_argument("-a", "--alpha", action="store", type=float, dest="alpha", default=0.05, help="significance threshold for p-values (only used when masking LFC values)")


    ###########################################
    # TIDE parser
    ###########################################

    TIDE_parser = subparser.add_parser("TIDE", help="performs TIDE analysis (Doughberty et al., 2021)", formatter_class=FORMATTER)
    
    TIDE_parser.add_argument("dea_file", action="store", help="file containing DEA results") 
    
    TIDE_parser.add_argument("-s", "--secretory", action="store_true", dest="secretory_flag", help="whether to also use the secretory tasks expansion")
    
    TIDE_parser.add_argument("-d", "--delim", action="store", type=str, dest="sep", default="\t", help="file field separator")
    
    TIDE_parser.add_argument("-o", "--out", action="store", type=str, dest="out_filename", default="tide_results.tsv", help="out file name to store permutation results")

    TIDE_parser.add_argument("-n", "--n_permutations", action="store", type=int, dest="n_permutations", default=1000, help="number of permutations")
    
    TIDE_parser.add_argument("--n_cpus", action="store", type=int, dest="n_cpus", default=1,help="number of CPUs for parallel processing")
    
    TIDE_parser.add_argument("--lfc_col", action="store", type=str, dest="lfc_col", default="log2FoldChange", help="column name containing LFC values")
    
    TIDE_parser.add_argument("--gene_col", action="store", type=str, dest="gene_col", default="geneID", help="column name containing gene symbols")
    
    TIDE_parser.add_argument("--or_func", action="store", type=str, dest="or_func", choices=["max", "absmax"], default="absmax", help="function to evaluate OR relationships in GPR rules")

    TIDE_parser.add_argument("--mask_lfc_values", action="store_true", dest="filter_lfc", help="flag to indicate whether to mask to 0 non-significant LFC values using their p-values")
    
    TIDE_parser.add_argument("--random_scores", action="store_true", dest="random_scores_flag", help="flag to indicate whether to return the random score distributions with the results")
    
    TIDE_parser.add_argument("--pvalue_col", action="store", type=str, dest="pvalue_col", default="pvalue", help="column name containing p-values (only used when masking LFC values)")
    
    TIDE_parser.add_argument("-a", "--alpha", action="store", type=float, dest="alpha", default=0.05, help="significance threshold for p-values (only used when masking LFC values)")


    ###########################################
    # CellFie parser
    ###########################################
    
    CellFie_parser = subparser.add_parser("CellFie", help="performs CellFie analysis (Richelle et al., 2021)", formatter_class=FORMATTER)
    
    CellFie_parser.add_argument("expr_file", action="store", help="file containing gene expression data")

    CellFie_parser.add_argument("-s", "--secretory", action="store_true", dest="secretory_flag", help="whether to use the secretory tasks expansion")
    
    CellFie_parser.add_argument("--gene_col", action="store", type=str, dest="gene_col", default="geneID", help="column in expression file containing genes")

    CellFie_parser.add_argument("-d", "--delim", action="store", type=str, dest="sep", default="\t", help="file field separator")
    
    CellFie_parser.add_argument("-o", "--out", action="store", type=str, dest="out_dir", default="CellFie_results", help="out directory name to store result files")
    
    CellFie_parser.add_argument("--binary_scores", action="store_true", dest="binary_scores_flag", help="flag to indicate whether also return the binary metabolic score matrix")
    
    CellFie_parser.add_argument("--threshold_type", action="store", type=str, dest="thresh_type", default="local", choices=["local","global"], help="determines the threshold approach to be used (a global approach used the same threshold for all genes whereas a local approach uses a different threshold for each gene when computing the gene activity levels)")
    
    CellFie_parser.add_argument("--global_threshold_type", action="store", type=str, dest="global_thresh_type", default="percentile", choices=["value","percentile"], help="option to determine whether to use a value or a percentile of the distribution of all genes as global treshold")

    CellFie_parser.add_argument("--global_value", action="store", type=float, dest="global_value", default=0.75, help="value to use as global threshold according to the global_threshold_type option (percentile values must be between 0 and 1)")
    
    CellFie_parser.add_argument("--local_threshold_type", action="store", type=str, dest="local_thresh_type", default="minmaxmean", choices=["minmaxmean","mean"], help="determines the threshold type to be used in a local approach ('minmaxmean': the threshold for each gene is determined by the mean of expression values across all conditions/samples BUT must be higher or equal than a lower bound and lower or equal to an upper bound; 'mean': the threshold of a gene is determined as its mean expression across all conditions/samples)")
    
    CellFie_parser.add_argument("--minmaxmean_threshold_type", action="store", type=str, dest="minmaxmean_thresh_type", default="percentile", choices=["percentile","value"], help="option to determine whether to use a value or a percentile of the distribution of all genes as upper and lower bounds")
    
    CellFie_parser.add_argument("--upper_bound", action="store", type=float, dest="upper_bound", default=0.75, help="upper bound value to be used according to the minmaxmean_threshold_type (percentile values must be between 0 and 1)")
    
    CellFie_parser.add_argument("--lower_bound", action="store", type=float, dest="lower_bound", default=0.25, help="lower bound value to be used according to the minmaxmean_threshold_type (percentile values must be between 0 and 1)") 

    
    return parser