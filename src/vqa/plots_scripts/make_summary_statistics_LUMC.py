import argparse
import json
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

def is_the_coverage_sufficient(reads, gr, low_limit=100):
    reads_region = reads[reads["Genomic_regions"] == gr]
    reads_coverage = len(reads_region)
    if reads_coverage >= low_limit:
        return True
    else:
        return False

def regions_with_sufficie_coverage_in_sample(reads, low_limit=100):
    genomic_regions = reads["Genomic_regions"].unique()
    regions_with_sufficient_coverage = []
    for gr in genomic_regions:
        if is_the_coverage_sufficient(reads, gr, low_limit):
            regions_with_sufficient_coverage.append(gr)
    return regions_with_sufficient_coverage

def get_grs_with_sufficient_coverage(samples, low_limit=100):
    sufficient_coverage_per_sample = {}

    for sample in samples:
        reads_file = "Experiments/lumc_subsample/" + sample + '/communities.tsv'
        reads = pd.read_csv(reads_file, sep="\t", header=0)
        sufficient_coverage_per_sample[sample] = regions_with_sufficie_coverage_in_sample(reads, low_limit=100)
    
    return sufficient_coverage_per_sample

def get_edit_distances(input_dir, strains, genomic_regions, subdirs, sufficient_coverage_across_samples):

    edit_distances = {}

    for subdir in subdirs:
        file = input_dir + '/' + subdir + '/consensus_lumc_comparison_post_processed.tsv'
        consensus_eval_file = pd.read_csv(file, sep='\t', header=0)
        edit_distances[subdir] = {}
        for strain in strains:
            edit_distances[subdir][strain] = {}

            for region in genomic_regions:
                if region not in sufficient_coverage_across_samples[subdir]:
                    continue
                edit_distances[subdir][strain][region] = np.nan

        for index, row in consensus_eval_file.iterrows():
            min_strain = row['Min_consensus']
            genomic_region = row['Genomic_region']
            edit_distance = row['Min_edit_distance']

            if genomic_region not in sufficient_coverage_across_samples[subdir]:
                continue
            if min_strain == "Both": 
                min_strain = "Wuhan"
            if min_strain == "BA1": 
                min_strain = "BA.1"
            edit_distances[subdir][min_strain][genomic_region] = edit_distance
   
    # get average per subdirectory
    
    ed_avg_per_subdir = {}

    for subdir in edit_distances.keys():
        values_exclude_nan = [x for x in list(edit_distances[subdir][strains[0]].values()) if (math.isnan(x) == False)]

        if len(values_exclude_nan) != 0:
            region_avg_seq_1 = sum(values_exclude_nan)/len(values_exclude_nan)
        else: 
            region_avg_seq_1 = np.nan
        
        values_exclude_nan = [x for x in list(edit_distances[subdir][strains[1]].values()) if (math.isnan(x) == False)]
        if len(values_exclude_nan) != 0:
            region_avg_seq_2 = sum(values_exclude_nan)/len(values_exclude_nan)
        else: 
            region_avg_seq_2 = np.nan
        
        if math.isnan(region_avg_seq_1) and math.isnan(region_avg_seq_2):  
            ed_avg_per_subdir[subdir] = np.nan
        elif math.isnan(region_avg_seq_1):
            ed_avg_per_subdir[subdir] = region_avg_seq_2
        elif math.isnan(region_avg_seq_2): 
            ed_avg_per_subdir[subdir] = region_avg_seq_1
        else: 
            ed_avg_per_subdir[subdir] = (region_avg_seq_1 + region_avg_seq_2)/2

    return ed_avg_per_subdir

def get_true_abundances():

    true_abundance = {
        "01_100": {"BA.1":0, "Wuhan": 1},
        "02_100": {"BA.1":0, "Wuhan": 1},
        "03_50": {"BA.1":0.5, "Wuhan": 0.5},
        "04_75": {"BA.1":0.25, "Wuhan": 0.75},
        "05_90": {"BA.1":0.1, "Wuhan": 0.9},
        "06_95": {"BA.1":0.05, "Wuhan": 0.95},
        "07_98": {"BA.1":0.02, "Wuhan": 0.98},
        "08_0": {"BA.1":1.0, "Wuhan": 0.0},
        "09_0": {"BA.1":1.0, "Wuhan": 0.0}
    }

    return true_abundance

def init_rel_abundances(samples, strains, genomic_regions, sufficient_coverage_across_samples):
    rel_abundances = {}
    for sample in samples:
        rel_abundances[sample] = {}
        for strain in strains:
            rel_abundances[sample][strain] = {}
            for gr in genomic_regions:
                if gr in sufficient_coverage_across_samples[sample]:
                    rel_abundances[sample][strain][gr] =0
    return rel_abundances

def get_relative_abundances(directory, samples, sufficient_coverage_across_samples, genomic_regions, strains, t_abs):
    rel_abundances = init_rel_abundances(samples, strains, genomic_regions, sufficient_coverage_across_samples)
    haplotype_assignment = {}

    for sample in samples:
            sample_dir = os.path.join(directory, sample)
            input_file = os.path.join(sample_dir, "consensus_lumc_comparison_post_processed.tsv")
            df = pd.read_csv(input_file, sep="\t", header=0)

            haplotype_assignment[sample] = {}

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

                if gr not in sufficient_coverage_across_samples[sample]:
                    continue

                if strain_assigned_to != "Both":   
                    rel_abundances[sample][strain_assigned_to][gr] += rel_ab
                    haplotype_assignment[sample][gr] = strain_assigned_to
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
                        haplotype_assignment[sample][gr] = "Wuhan"
                    else:
                        rel_abundances[sample]["BA.1"][gr] += rel_ab
                        haplotype_assignment[sample][gr] = "BA.1"
    
    return rel_abundances, haplotype_assignment

def get_rel_ab_error(true_abundances, rel_abundances, directory):

    rel_ab_errors = {}
    for sample in true_abundances.keys():
        sample_dir = os.path.join(directory, sample)
        input_file = os.path.join(sample_dir, "consensus_lumc_comparison_post_processed.tsv")
        df = pd.read_csv(input_file, sep="\t", header=0)
        df_same = list(df[df["Edit_distance_between_consensus"] == 0]['Genomic_region'])
        rel_ab_errors[sample] = {}
            
        for strain in true_abundances[sample].keys():
            
            if strain not in rel_ab_errors[sample].keys(): 
                rel_ab_errors[sample][strain] = {}

            for gr in rel_abundances[sample][strain].keys():
                
                rel_ab_errors[sample][strain][gr] = abs(true_abundances[sample][strain] - rel_abundances[sample][strain][gr])/ ((true_abundances[sample][strain]) + 0.0000001)
                
                if gr in df_same:
                    min_distance = np.argmin([abs(1 - rel_abundances[sample][strain][gr]), abs(0 - rel_abundances[sample][strain][gr])]) 
                    if min_distance == 0:
                        rel_ab_errors[sample][strain][gr] = abs(1 - rel_abundances[sample][strain][gr]) / ((1) + 0.0000001)
                    else:
                        rel_ab_errors[sample][strain][gr] = abs(0 - rel_abundances[sample][strain][gr]) / ((0) + 0.0000001)

    # get average per subdirectory
    rel_ab_errors_avg_per_subdir_and_strain = {}

    for subdir in rel_ab_errors.keys():
        rel_ab_errors_avg_per_subdir_and_strain[subdir] = {}
        for strain in rel_ab_errors[subdir].keys(): 
            rel_ab_errors_avg_per_subdir_and_strain[subdir][strain] = sum(rel_ab_errors[subdir][strain].values())/len(rel_ab_errors[subdir][strain].values())
    
    # get average per subdirectory
    rel_ab_errors_avg_per_subdir = {}

    for subdir in rel_ab_errors_avg_per_subdir_and_strain.keys():
        rel_ab_errors_avg_per_subdir[subdir] = sum(rel_ab_errors_avg_per_subdir_and_strain[subdir].values())/len(rel_ab_errors_avg_per_subdir_and_strain[subdir].values())
    
    return rel_ab_errors_avg_per_subdir

def get_estimated_haplotypes_numbers(input_dir, subdirs, regions_with_sufficient_coverage):
    n_haplotypes = {}

    for subdir in subdirs:
        file = input_dir + '/' + subdir + '/consensus_lumc_comparison_post_processed.tsv'
        n_haplotypes[subdir] = {}
       # load the file if it exists
        if os.path.exists(file):
            consensus_tsv = pd.read_csv(file, sep='\t', header=0)
        
        for gr in consensus_tsv['Genomic_region'].unique():
            if gr not in regions_with_sufficient_coverage[subdir]:
                continue
            n_haplotypes[subdir][gr] = len(consensus_tsv[consensus_tsv['Genomic_region'] == gr])

    # get average number of haplotypes per subdirectory
    n_haplotypes_avg_per_subdir = {}

    for subdir in n_haplotypes.keys():
        n_haplotypes_avg_per_subdir[subdir] = sum(n_haplotypes[subdir].values())/len(n_haplotypes[subdir].values())
        
    return n_haplotypes_avg_per_subdir

def get_estimated_haplotypes_numbers(input_dir, subdirs, regions_with_sufficient_coverage):
    n_haplotypes = {}

    for subdir in subdirs:
        file = input_dir + '/' + subdir + '/consensus_lumc_comparison_post_processed.tsv'
        n_haplotypes[subdir] = {}
       # load the file if it exists
        if os.path.exists(file):
            consensus_tsv = pd.read_csv(file, sep='\t', header=0)
        
        for gr in consensus_tsv['Genomic_region'].unique():
            if gr not in regions_with_sufficient_coverage[subdir]:
                continue
            n_haplotypes[subdir][gr] = len(consensus_tsv[consensus_tsv['Genomic_region'] == gr])

    # get average number of haplotypes per subdirectory
    n_haplotypes_avg_per_subdir = {}

    for subdir in n_haplotypes.keys():
        n_haplotypes_avg_per_subdir[subdir] = sum(n_haplotypes[subdir].values())/len(n_haplotypes[subdir].values())
        
    return n_haplotypes_avg_per_subdir

def get_recall_and_duplication_ratio(input_dir, genomic_regions, subdirs, strains, haploype_assignment, sufficient_coverage_across_samples):
    recall = {} # number of haplotypes that were correctly identified
    duplication_ratio = {} # number of haplotypes that were identified more than once
    
    for subdir in subdirs:
        file = input_dir + '/' + subdir + '/consensus_lumc_comparison_post_processed.tsv'
        recall[subdir] = {}
        duplication_ratio[subdir] = {}

        # load the file if it exists
        if os.path.exists(file):
            consensus_tsv = pd.read_csv(file, sep='\t', header=0)
 
        for strain in strains:
            recall[subdir][strain] = {}
            for gr in genomic_regions:
                if gr not in sufficient_coverage_across_samples[subdir]: 
                    continue
                recall[subdir][strain][gr] = False
                duplication_ratio[subdir][gr] = 0
        
        for index, row in consensus_tsv.iterrows():
            sequence_id = row['Min_consensus']
            genomic_region = row['Genomic_region']

            if genomic_region not in sufficient_coverage_across_samples[subdir]: 
                continue

            if sequence_id == "Both":
                if haploype_assignment[subdir][genomic_region] == "Wuhan":
                    sequence_id = "Wuhan"
                else:
                    sequence_id = "BA.1"

            if sequence_id == "BA1": 
                sequence_id = "BA.1"
            
            recall[subdir][sequence_id][genomic_region] = True
            
            diff_between_true_haps = row['Edit_distance_between_consensus']
            
            if (diff_between_true_haps == 0):
                recall[subdir]["Wuhan"][genomic_region] = True
                recall[subdir]["BA.1"][genomic_region] = True


            duplication_ratio[subdir][genomic_region] += 1
        
    # get the recall and duplication ratio per subdirectory
    recall_avg_per_subdir_and_seqs = {}
    duplication_ratio_avg_per_subdir = {}
    recall_per_subdir = {}

    for subdir in recall.keys():
        recall_avg_per_subdir_and_seqs[subdir] = {}
        
        true_values = [1 if value == True else 0 for value in recall[subdir][strains[0]].values()]
        recall_seq_1 = sum(true_values)/len(true_values)

        true_values = [1 if value == True else 0 for value in recall[subdir][strains[1]].values()]
        recall_seq_2 = sum(true_values)/len(true_values)

        recall_avg_per_subdir_and_seqs[subdir][strains[0]] = recall_seq_1
        recall_avg_per_subdir_and_seqs[subdir][strains[1]] = recall_seq_2
        recall_per_subdir[subdir] = (recall_seq_1 + recall_seq_2)/2

        avges_subdir = []
        for gr in duplication_ratio[subdir].keys():
            if gr not in sufficient_coverage_across_samples[subdir]:
                continue
            if subdir in ["08_0", "09_0", "01_100", "02_100"]:
                expected_n_haplotypes = 1
                temp_avg = duplication_ratio[subdir][gr]/expected_n_haplotypes
            else:
                # take first row of consensus_tsv matching the genomic region and check Edit_distance_between_consensu
                if duplication_ratio[subdir][gr] == 0:
                    avges_subdir.append(0)
                    continue

                file = input_dir + '/' + subdir + '/consensus_lumc_comparison_post_processed.tsv'

                consensus_tsv = pd.read_csv(file, sep='\t', header=0)
                diff_between_true_haps = list(consensus_tsv[consensus_tsv["Genomic_region"] == gr]["Edit_distance_between_consensus"])[0]
                expected_n_haplotypes = 2 if diff_between_true_haps > 0 else 1
                temp_avg = duplication_ratio[subdir][gr]/expected_n_haplotypes
            avges_subdir.append(temp_avg)
        duplication_ratio_avg_per_subdir[subdir] = sum(avges_subdir)/len(avges_subdir)
        
    return recall_per_subdir, duplication_ratio_avg_per_subdir, recall_avg_per_subdir_and_seqs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    args = parser.parse_args()

    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                    (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                    (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                        (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                    (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                        (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                        (26766, 27872), (27808, 28985), (28699, 29768)]

    genomic_regions = [ str(gr[0]) + "_" + str(gr[1]) for gr in genomic_regions]

    strains = ["Wuhan", "BA.1"]

    subdirs = ["01_100", "02_100" ,"03_50", "04_75", "05_90",'06_95',"07_98", "08_0", "09_0"]
    
    sufficient_coverage_across_samples = get_grs_with_sufficient_coverage(subdirs, low_limit=100)

    ed_avges = get_edit_distances(args.input_dir, strains, genomic_regions, subdirs, sufficient_coverage_across_samples)


    true_abundances = get_true_abundances()

    rel_abundances, haplotype_assignment = get_relative_abundances(args.input_dir, subdirs, sufficient_coverage_across_samples, genomic_regions, strains, true_abundances)


    relative_abundance_error = get_rel_ab_error(true_abundances, rel_abundances, args.input_dir)

    n_haplotypes = get_estimated_haplotypes_numbers(args.input_dir, subdirs, sufficient_coverage_across_samples)
    
    recall, duplication_ratio, recall_per_subdir_and_seqs = get_recall_and_duplication_ratio(args.input_dir, genomic_regions, subdirs, strains, haplotype_assignment, sufficient_coverage_across_samples)

    # make dataframe where rows are the datasets and columns are the metrics
    df = pd.DataFrame(columns=['Dataset', 'Relative_abundance_error', 'Edit_distance','Number_of_haplotypes', 'Recall', 'Recall_{}'.format(strains[0]), 'Recall_{}'.format(strains[1]), 'Duplication_ratio'])
    for subdir in ed_avges.keys():
        # round to 3 decimal places
        df = df._append({'Dataset': subdir, 'Relative_abundance_error': round(relative_abundance_error[subdir], 3), 'Edit_distance': round(ed_avges[subdir], 3), 'Number_of_haplotypes': round(n_haplotypes[subdir], 3), 'Recall': round(recall[subdir], 3), 'Recall_{}'.format(strains[0]): round(recall_per_subdir_and_seqs[subdir][strains[0]], 3), 'Recall_{}'.format(strains[1]): round(recall_per_subdir_and_seqs[subdir][strains[1]], 3), 'Duplication_ratio': round(duplication_ratio[subdir], 3)}, ignore_index=True)

    # save the dataframe to a file
    df.to_csv(args.output_dir + '/summary_statistics_extensive_LUMC.csv', sep=',', index=False)
    

if __name__ == "__main__":
    sys.exit(main())