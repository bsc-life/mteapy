#!/usr/bin/env python3

import os
import pandas as pd
from time import time
from cobra.io import read_sbml_model

from mteapy.parser import mtea_parser
from mteapy.colors import bcolors as bc

from mteapy.utils import print_banner, mask_lfc_values, add_task_metadata, check_ensemblid
from mteapy.tide import compute_TIDEe, compute_TIDE
from mteapy.cellfie import compute_CellFie


def main() -> None:

    curdir = os.path.dirname(os.path.realpath(__file__))
    parser = mtea_parser()
    args = parser.parse_args()
    start = time()
    print_banner()
    print()
    
    if args.citation_flag:
        print(f"{bc.CYAN}Cite MTEApy:{bc.ENDC}\thttps://doi.org/10.1101/2025.05.08.652850")
        print()
        print(f"{bc.CYAN}Cite TIDE:{bc.ENDC}\thttps://doi.org/10.1016/j.celrep.2021.108836")
        print(f"{bc.CYAN}Cite TIDE-e:{bc.ENDC}\thttps://doi.org/10.1101/2025.05.08.652850")
        print(f"{bc.CYAN}Cite CellFie:{bc.ENDC}\thttps://doi.org/10.1016/j.crmeth.2021.100040\n")
        exit(1)
        
    
    ###########################################
    # Performs TIDE (and TIDE-essential)
    ###########################################
    
    if args.command == "TIDE":
        
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
        task_metadata = pd.read_csv(os.path.join(curdir, "../data/task_metadata.tsv"), sep="\t")
        task_structure = pd.read_csv(os.path.join(curdir, "../data/task_structure_matrix.tsv"), \
                                     sep="\t", index_col=0)
        gene_essentiality = pd.read_csv(os.path.join(curdir, "../data/HumanGEM_essential_genes_matrix.tsv"),\
                                        delimiter="\t", index_col=0)
        
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
        if not check_ensemblid(expr_data_df[args.gene_col].values):
            print(f"{bc.FAIL}ERROR: one or more genes in the column {args.gene_col} are not EnsemblIDs.{bc.ENDC}\n")
        
        # Loading in metabolic model
        print("Loading metabolic model", end=" ")
        model = read_sbml_model(os.path.join(curdir, "../data/HumanGEM.xml.gz"))
        print("- OK.")        

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
        print(f"\tNº jobs      = {args.n_cpus}")
        print(f"\tPermutations = {args.n_permutations}")
        print(f"\tOR function  = {args.or_func}")
        
        if args.random_seed is not None:
            print(f"\tRandom seed  = {args.random_seed}")


        TIDE_results = compute_TIDE(
            expr_data_df.set_index(args.gene_col), 
            args.lfc_col,
            task_structure, 
            model, 
            args.or_func,
            args.n_permutations,
            args.n_cpus,
            args.random_scores_flag,
            args.random_seed
        )
        TIDE_e_results = compute_TIDEe(
            expr_data_df.set_index(args.gene_col), 
            args.lfc_col,
            gene_essentiality, 
            args.n_permutations,
            args.n_cpus,
            args.random_scores_flag,
            args.random_seed
        )
        print("Saving results", end = " ")
        TIDE_results.rename(columns={"score": "TIDE_score", "pvalue": "TIDE_pvalue"}, inplace=True)
        TIDE_e_results.rename(columns={"score": "TIDEe_score", "pvalue": "TIDEe_pvalue"}, inplace=True)
        
        # Merging TIDE and TIDE-essential results using left join
        all_results = pd.merge(TIDE_results, TIDE_e_results, on="task_id", how="left")
        all_results = add_task_metadata(all_results, task_metadata)
        all_results.to_csv(args.out_filename, index=False, sep="\t", na_rep="NA")
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
        task_metadata = pd.read_csv(os.path.join(curdir, "../data/task_metadata.tsv"), sep="\t")
        task_structure = pd.read_csv(os.path.join(curdir, "../data/task_structure_matrix.tsv"), \
                                     sep="\t", index_col=0)
        
        # Column names check
        if args.gene_col not in expr_data_df.columns:
            print(f"{bc.FAIL}ERROR: Gene column '{args.gene_col}' not found in input file.{bc.ENDC}\n")
            exit(1)
        if any(expr_data_df[args.gene_col].duplicated()):
            print(f"{bc.FAIL}ERROR: gene column contains duplicated entries.{bc.ENDC}\n")
            exit(1)
        if not check_ensemblid(expr_data_df[args.gene_col].values):
            print(f"{bc.FAIL}ERROR: one or more genes in the column {args.gene_col} are not EnsemblIDs.{bc.ENDC}\n")
            exit(1)
        
        # Loading in metabolic model
        print("Loading metabolic model", end=" ")
        model = read_sbml_model(os.path.join(curdir, "../data/HumanGEM.xml.gz"))
        print("- OK.")

        # Filtering genes not in model
        model_genes = [str(gene) for gene in model.genes]
        print(f"Nº genes in metabolic model: {len(model_genes)}.")
        
        expr_data_df.set_index(args.gene_col, inplace=True)
        initial_genes = expr_data_df.index.values 
        print(f"Nº genes in input file: {len(initial_genes)}.")
        
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
        
        # Adding metadata
        metabolic_scores_df = add_task_metadata(metabolic_scores_df, task_metadata)
        binary_scores_df = add_task_metadata(binary_scores_df, task_metadata)

        # Saving results
        print(f"Saving results", end=" ")
        metabolic_scores_df.to_csv(f"{args.out_dir}/cellfie_scores.tsv", sep="\t", index=False)
        if args.binary_scores_flag:
            binary_scores_df.to_csv(f"{args.out_dir}/cellfie_binary_scores.tsv", sep="\t", index=False)
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