import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches


def true_abundances_func(pures):

    if pures: 

        true_abundance = {
            "01_100": {"BA.1":0, "Wuhan": 1},
            "02_100": {"BA.1":0, "Wuhan": 1},
            "08_0": {"BA.1":1.0, "Wuhan": 0.0},
            "09_0": {"BA.1":1.0, "Wuhan": 0.0}
            }
    else: 

        true_abundance = {
            "03_50": {"BA.1":0.5, "Wuhan": 0.5},
            "04_75": {"BA.1":0.25, "Wuhan": 0.75},
            "05_90": {"BA.1":0.1, "Wuhan": 0.9},
            "06_95": {"BA.1":0.05, "Wuhan": 0.95},
            "07_98": {"BA.1":0.02, "Wuhan": 0.98}
        }
    
    return true_abundance

def init_rel_abundances(samples, strains, genomic_regions):
    rel_abundances = {}
    for sample in samples:
        rel_abundances[sample] = {}
        for strain in strains:
            rel_abundances[sample][strain] = {}
            for gr in genomic_regions:
                rel_abundances[sample][strain][gr] =0
    return rel_abundances

def calculate_rel_ab_error(true_abundances, rel_abundances):
    rel_ab_errors = {}
    for sample in true_abundances.keys():
        rel_ab_errors[sample] = {}
        for strain in true_abundances[sample].keys():
            rel_ab_errors[sample][strain] = {}
            for gr in rel_abundances[sample][strain].keys():
                if true_abundances[sample][strain] == 0:
                    if abs(true_abundances[sample][strain] - rel_abundances[sample][strain][gr]) == 0:
                        rel_ab_errors[sample][strain][gr] = 0
                    else:
                        rel_ab_errors[sample][strain][gr] = 1
                    
                    continue
                    
                rel_ab_errors[sample][strain][gr] = abs(true_abundances[sample][strain] - rel_abundances[sample][strain][gr])/ (true_abundances[sample][strain])
                                                    
    return rel_ab_errors

def main():
    directory = "Experiments/lumc_subsample/"
    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                    (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                    (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                    (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                    (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                    (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                    (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]
    
    genomic_regions = [f"{gr[0]}_{gr[1]}" for gr in genomic_regions]
    
    samples =  ["03_50", "04_75", "05_90", "06_95", "07_98"] #["01_100", "03_50", "04_75", "05_90","07_98", "08_0", "09_0"]
    pure_samples = ["01_100", "02_100", "08_0", "09_0"]
    rel_abundances = init_rel_abundances(samples, ["BA.1", "Wuhan"], genomic_regions)
    t_abs = true_abundances_func(pures=False)

    for sample in samples:
        sample_dir = os.path.join(directory, sample)
        input_file = os.path.join(sample_dir, "consensus_lumc_comparison_post_processed.tsv")
        df = pd.read_csv(input_file, sep="\t", header=0)

        #sort such that "Both" is last in the list per genomic region
        order = ["BA.1", "Wuhan", "Both"]
        # sort the dataframe based on the order of the strains
        # If any entry is "BA1" replace with "BA.1"
        df["Min_consensus"] = df["Min_consensus"].replace("BA1", "BA.1")


        df["Min_consensus"] = pd.Categorical(df["Min_consensus"], categories=order, ordered=True)

        df = df.sort_values(by=["Genomic_region", "Min_consensus"])
       
        for index, row in df.iterrows():
            gr = row["Genomic_region"]
            strain_assigned_to = row["Min_consensus"]
            rel_ab = row["Relative_abundance"]
            

            
            if strain_assigned_to != "Both":   
                rel_abundances[sample][strain_assigned_to][gr] += rel_ab
            else:
                wuhan_true = t_abs[sample]["Wuhan"]
                ba1_true = t_abs[sample]["BA.1"]
                current_wuhan = rel_abundances[sample]["Wuhan"][gr]
                current_ba1 = rel_abundances[sample]["BA.1"][gr]
                leftover_wuhan = wuhan_true - current_wuhan
                leftover_ba1 = ba1_true - current_ba1

                abs_diff_wuhan = abs(leftover_wuhan - rel_ab)
                abs_diff_ba1 = abs(leftover_ba1 - rel_ab)


                if abs_diff_wuhan < abs_diff_ba1:
                    rel_abundances[sample]["Wuhan"][gr] += rel_ab
                else:
                    rel_abundances[sample]["BA.1"][gr] += rel_ab
    
    
    rel_ab_errors = calculate_rel_ab_error(t_abs, rel_abundances)
    
    
    
    

if __name__ == '__main__':
    main()
