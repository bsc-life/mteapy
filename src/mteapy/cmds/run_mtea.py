#!/usr/bin/env python3

import os
import pandas as pd
from time import time
from cobra.io import read_sbml_model

from mteapy.parser import mtea_parser
from mteapy.colors import bcolors as bc

from mteapy.utils import mask_lfc_values
from mteapy.tide import compute_TIDEe, compute_TIDE
from mteapy.cellfie import compute_CellFie


def main() -> None:

    curdir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(curdir, "../../../logo_ascii.txt"), "r") as f:
        logo = f.read()
        print(logo)
    
    parser = mtea_parser()
    args = parser.parse_args()
    start = time()
    print()
    
    if args.citation_flag:
        print(f"{bc.CYAN}Cite TIDE:{bc.ENDC}\thttps://doi.org/10.1016/j.celrep.2021.108836")
        print(f"{bc.CYAN}Cite CellFie:{bc.ENDC}\thttps://doi.org/10.1016/j.crmeth.2021.100040\n")
        exit(1)
   
    if args.task_metadata:
        task_metadata = pd.read_csv(os.path.join(curdir, "../../data/task_metadata.tsv"), sep="\t", index_col=0)
        task_metadata.to_csv("./task_metadata.tsv", sep="\t")
        print("Downloaded metabolic task metadata at ./task_metadata.tsv\n")
        exit(1)

    if args.task_metadata_sec:
        task_metadata_sec = pd.read_csv(os.path.join(curdir, "../../data/task_metadata_sec.tsv"), sep="\t", index_col=0)
        task_metadata_sec.to_csv("./task_metadata_sec.tsv", sep="\t")
        print("Downloaded secretory metabolic task metadata at ./task_metadata_sec.tsv\n")
        exit(1)

        
    ###########################################
    # Performs TIDE with essential genes (TIDE-essential)
    ###########################################
    
    if args.command == "TIDE-essential":
        
        print("Starting Tasks Inferred from Differential Expression (TIDE) analysis using  essential genes.")
        
        # File and dir status check
        if not os.path.isfile(args.dea_file):
            print(f"{bc.FAIL}ERROR: File {args.dea_file} could not be found.{bc.FAIL}\n")
            exit(1)
        
        # Out file handling
        if os.path.isfile(args.out_filename):
            print(f"File {args.out_filename} already exists, will be overwritten.")
        elif not os.path.isdir(os.path.dirname(args.out_filename)) and os.path.dirname(args.out_filename) != "":
            print(f"Creating out directory {os.path.dirname(args.out_filename)}")
            os.mkdir(os.path.dirname(args.out_filename))
        if os.path.isdir(args.out_filename):
            args.out_filename += "tide-essential-results.tsv"
            print(f"Results will be saved into {args.out_filename}")
        else:
            print(f"Results will be saved into {args.out_filename}.")
        
        # Reading in data and gene essentiality
        expr_data_df = pd.read_csv(args.dea_file, delimiter=args.sep)
        gene_essentiality = pd.read_csv(os.path.join(curdir, "../../data/HumanGEM_essential_genes_matrix.tsv"),\
                                        delimiter="\t", index_col=0)

        # Column names check
        if args.lfc_col not in list(expr_data_df): 
            print(f"{bc.FAIL}ERROR: LFC column name '{args.lfc_col}' not found in results file.{bc.ENDC}\n")
            exit(1)
        if args.gene_col not in list(expr_data_df):
            print(f"{bc.FAIL}ERROR: Gene column name '{args.gene_col}' not found in results file.{bc.ENDC}\n")
            exit(1)
        if any(expr_data_df[args.gene_col].duplicated()):
            print(f"{bc.FAIL}ERROR: gene column contains duplicated entries.{bc.ENDC}\n")
            exit(1)
            
        # Filtering LFC
        if args.filter_lfc:
            if args.pvalue_col not in list(expr_data_df):
                print(f"{bc.FAIL}ERROR: P-value column name {args.pvalue_col} not in results file {args.dea_file}.{bc.ENDC}\n")
                exit(1)
            
            print(f"LFC values will be masked using their p-values (alpha = {args.alpha})", end = " ")   
            expr_data_df = mask_lfc_values(expr_data_df, args.lfc_col, args.pvalue_col, args.alpha)
            print("- OK")
            
        # Main execution
        print("Starting analysis:")
        print(f"\tNº jobs      = {args.n_jobs}")
        print(f"\tPermutations = {args.n_permutations}")

        print("Saving results", end=" ")
        TIDE_e_results = compute_TIDEe(
            expr_data_df.set_index(args.gene_col), 
            args.lfc_col,
            gene_essentiality, 
            args.n_permutations,
            args.n_jobs,
            args.random_scores_flag
        )
        TIDE_e_results.sort_values(by="pvalue").to_csv(args.out_filename, index=False, sep="\t")
        print("- OK.")
        
    
    ###########################################
    # Performs TIDE
    ###########################################
    
    elif args.command == "TIDE":
        
        print("Starting Tasks Inferred from Differential Expression (TIDE) analysis.")

        # File and dir status check
        if not os.path.isfile(args.dea_file):
            print(f"{bc.FAIL}ERROR: File {args.dea_file} could not be found.{bc.FAIL}\n")
            exit(1)

        # Out file handling
        if os.path.isfile(args.out_filename):
            print(f"File {args.out_filename} already exists, will be overwritten.")
        elif not os.path.isdir(os.path.dirname(args.out_filename)) and os.path.dirname(args.out_filename) != "":
            print(f"Creating out directory {os.path.dirname(args.out_filename)}")
            os.mkdir(os.path.dirname(args.out_filename))
        if os.path.isdir(args.out_filename):
            args.out_filename += "tide-results.tsv"
            print(f"Results will be saved into {args.out_filename}")
        else:
            print(f"Results will be saved into {args.out_filename}.")
    
        # Reading in data, model and task structure
        # TODO: allow user to input their own metabolic model and task structures
        expr_data_df = pd.read_csv(args.dea_file, delimiter=args.sep)
        task_structure = pd.read_csv(os.path.join(curdir, "../../data/task_structure_matrix.tsv"), \
                                     sep="\t", index_col=0)
        
        print("Loading metabolic model", end=" ")
        if args.secretory_flag:
            model = read_sbml_model(os.path.join(curdir, "../../data/HumanGEM_secretory.xml.gz"))
            task_structure_sec = pd.read_csv(os.path.join(curdir, "../../data/task_structure_matrix_sec.tsv"), \
                                             sep="\t", index_col=0)
            task_structure = pd.concat([task_structure, task_structure_sec]).fillna(0)
        else:
            model = read_sbml_model(os.path.join(curdir, "../../data/HumanGEM.xml.gz"))
        print("- OK.")
                
        # Column names check
        if args.lfc_col not in list(expr_data_df): 
            print(f"{bc.FAIL}ERROR: Log2-FC column '{args.lfc_col}' not found in input file.{bc.ENDC}\n")
            exit(1)
        if args.gene_col not in list(expr_data_df):
            print(f"{bc.FAIL}ERROR: Gene column '{args.gene_col}' not found in input file.{bc.ENDC}\n")
            exit(1)
        if any(expr_data_df[args.gene_col].duplicated()):
            print(f"{bc.FAIL}ERROR: gene column contains duplicated entries.{bc.ENDC}\n")
            exit(1)

        # Filtering LFC
        if args.filter_lfc:
            if args.pvalue_col not in list(expr_data_df):
                print(f"{bc.FAIL}ERROR: P-value column name {args.pvalue_col} not in results file {args.dea_file}.{bc.ENDC}\n")
                exit(1)
            
            print(f"LFC values will be masked using their p-values (alpha = {args.alpha})", end = " ")   
            expr_data_df = mask_lfc_values(expr_data_df, args.lfc_col, args.pvalue_col, args.alpha)
            print("- OK.")

        # Main execution
        print("Starting analysis:")
        print(f"\tNº jobs      = {args.n_jobs}")
        print(f"\tPermutations = {args.n_permutations}")
        print(f"\tOR function  = {args.or_func}")
        if args.secretory_flag:
            print(f"\tModules      = metabolic & secretory")
        else:
            print(f"\tModules      = metabolic")
        
        print("Saving results", end = " ")
        TIDE_results = compute_TIDE(
            expr_data_df.set_index(args.gene_col), 
            args.lfc_col,
            task_structure, 
            model, 
            args.or_func,
            args.n_permutations,
            args.n_jobs,
            args.random_scores_flag
        )
        TIDE_results.sort_values(by="pvalue").to_csv(args.out_filename, index=False, sep="\t")
        print("- OK.")
        
        
    ###########################################
    # Performs CellFie
    ###########################################
    
    elif args.command == "CellFie":
        
        print("Starting CellFie analysis.")

        # File and dir status check
        if not os.path.isfile(args.expr_file):
            print(f"{bc.FAIL}ERROR: File {args.expr_file} could not be found.{bc.ENDC}\n")
            exit(1)
        
        # Out directory handling
        if not os.path.isdir(args.out_dir):
            os.mkdir(args.out_dir)
        print(f"Results will be saved into {args.out_dir}")
        
        # Reading in data
        expr_data_df = pd.read_csv(args.expr_file, delimiter=args.sep)
        task_structure = pd.read_csv(os.path.join(curdir, "../../data/task_structure_matrix.tsv"), \
                                     sep="\t", index_col=0)
        
        print("Loading metabolic model", end=" ")
        if args.secretory_flag:
            model = read_sbml_model(os.path.join(curdir, "../../data/HumanGEM_secretory.xml.gz"))
            task_structure_sec = pd.read_csv(os.path.join(curdir, "../../data/task_structure_matrix_sec.tsv"), \
                                             sep="\t", index_col=0)
            task_structure = pd.concat([task_structure, task_structure_sec]).fillna(0)
        else:
            model = read_sbml_model(os.path.join(curdir, "../../data/HumanGEM.xml.gz"))
        print("- OK.")
        
        # Column names check
        if args.gene_col not in expr_data_df.columns:
            print(f"{bc.FAIL}ERROR: Gene column '{args.gene_col}' not found in input file.{bc.ENDC}\n")
            exit(1)
        if any(expr_data_df[args.gene_col].duplicated()):
            print(f"{bc.FAIL}ERROR: gene column contains duplicated entries.{bc.ENDC}\n")
            exit(1)

        # Filtering genes not in model
        model_genes = [str(gene) for gene in model.genes]
        print(f"Nº genes in metabolic model: {len(model_genes)}).")
        
        expr_data_df.set_index(args.gene_col, inplace=True)
        initial_genes = expr_data_df.index.values 
        expr_data_df = expr_data_df.loc[[gene for gene in model_genes if gene in expr_data_df.index]]
        print(f"A total of {len(initial_genes)-len(expr_data_df.index)} genes were not found in the model and were removed.")

        # Thresholds strategy summary
        print("Starting analysis:")
        print(f"\tNº samples detected   = {len(expr_data_df.columns)}")
        print(f"\tThreshold type        = {args.thresh_type}")
        if args.thresh_type == "global":
            print(f"\tGlobal approach       = {args.global_thresh_type}")
            print(f"\tThreshold value       = {args.global_value}")
        elif args.thresh_type == "local":
            print(f"\tLocal approach        = {args.local_thresh_type}")
            if args.local_thresh_type == "minmaxmean":
                print(f"\tUpper bound           = {args.upper_bound}")
                print(f"\tLower bound           = {args.lower_bound}")
        if args.secretory_flag:
            print(f"\tModules               = metabolic & secretory")
        else:
            print(f"\tModules               = metabolic")
        
        # Main execution
        metabolic_scores_df, binary_scores_df = compute_CellFie(
            expr_data_df, 
            task_structure, 
            model, 
            args.thresh_type, 
            args.local_thresh_type, 
            args.minmaxmean_thresh_type, 
            args.upper_bound, args.lower_bound, 
            args.global_thresh_type, args.global_value
        )
        
        # Saving results
        print(f"Saving results", end=" ")
        metabolic_scores_df.to_csv(f"{args.out_dir}/cellfie-scores.tsv", sep="\t")
        if args.binary_scores_flag:
            binary_scores_df.to_csv(f"{args.out_dir}/cellfie-binary-scores.tsv", sep="\t")
        print("- OK.")
        
    
    else:
        print("Command line tool to perform Metabolic Task Enrichment Analysis (MTEA) using several methods.", end=" ")
        print("All analysis are based on the Human1 metabolic model (Robinson et al., 2020).", end=" ")
        print("Uses metabolic task list from Richelle et. al 2021 curated to Human1 (use --task_metadata to download list of tasks).\n")
        print("Usage: run-mtea [-h] [-c] [-t] [-s] { TIDE-essential, TIDE, CellFie } [command options]")
        print(f"{bc.FAIL}Select one of the available commands (see --help for more information).{bc.ENDC}\n")
        exit(1)
    
    elapsed = time() - start
    print(f"\n{bc.OKGREEN}Elapsed time:\t{elapsed:.3f} seconds.{bc.ENDC}")
    print()


# MAIN
if __name__ == "__main__":
    main()