import argparse
import json
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

## TO BE ADJUSTED FOR LUMC


def get_edit_distances(input_dir, seq_ids, genomic_regions):

    edit_distances = {}

    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
         # make sure is a directory
        if not os.path.isdir(input_dir + '/' + subdir):
            continue
        file = input_dir + '/' + subdir + '/consensus_evaluation.tsv'
        consensus_eval_file = pd.read_csv(file, sep='\t', header=0)
        edit_distances[subdir] = {}
        for seq_id in seq_ids:
            edit_distances[subdir][seq_id] = {}

            for region in genomic_regions:
                edit_distances[subdir][seq_id][region] = np.nan

        for index, row in consensus_eval_file.iterrows():
            sequence_id = row['Sequence_id']
            genomic_region = row['Genomic_region']
            edit_distance = row['Edit_distance']

            edit_distances[subdir][sequence_id][genomic_region] = edit_distance
   
    # get average per subdirectory
    
    ed_avg_per_subdir = {}

    for subdir in edit_distances.keys():
        values_exclude_nan = [x for x in list(edit_distances[subdir][seq_ids[0]].values()) if (math.isnan(x) == False)]

        if len(values_exclude_nan) != 0:
            region_avg_seq_1 = sum(values_exclude_nan)/len(values_exclude_nan)
        else: 
            region_avg_seq_1 = np.nan
        
        values_exclude_nan = [x for x in list(edit_distances[subdir][seq_ids[1]].values()) if (math.isnan(x) == False)]
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

def get_accuracy(input_dir): 
    accuracies = {}
    # open the file if it exists
    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
        # make sure is a directory
        if not os.path.isdir(input_dir + '/' + subdir):
            continue
        file = input_dir + '/' + subdir + '/predictions.tsv'
        predictions_tsv = pd.read_csv(file, sep='\t', header=0)
        accuracies[subdir] =  sum(predictions_tsv['Predicted_label'] == predictions_tsv['True_label'])/len(predictions_tsv)
    return accuracies

def get_true_abundances(input_dir, seq_ids, genomic_regions):
    true_abundances = {}

    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
        # make sure is a directory
        if not os.path.isdir(input_dir + '/' + subdir):
            continue
        file = input_dir + '/' + subdir + '/mixture.json'
        with open(file) as f:
            data = json.load(f)
            for key in data:
                if key not in true_abundances.keys():
                    true_abundances[key] = {}
                true_abundances[key][subdir] = data[key]

    return true_abundances

def get_relative_abundance_error(input_dir, seq_ids, genomic_regions, true_abundances):
    
    predicted_abundance = {}

    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
        # make sure is a directory
        if not os.path.isdir(input_dir + '/' + subdir):
            continue
        file = input_dir+ '/' + subdir + '/' + 'consensus.tsv'
        predicted_abundance[subdir] = {}
        for seq_id in seq_ids:
            predicted_abundance[subdir][seq_id] = {}
            for gr in genomic_regions:
                predicted_abundance[subdir][seq_id][gr] = 0

        # open the file if it exists
        if os.path.exists(file):
            df = pd.read_csv(file, sep='\t', header=0)
       
            for index, row in df.iterrows():
                sequence_id = row['Sequence_id']
                genomic_region = row['Genomic_region']
                abundance = row['Relative_abundance']

                predicted_abundance[subdir][sequence_id][genomic_region] = abundance
    
    relative_abundance_error = {}

    
    for subdir in predicted_abundance.keys():
        relative_abundance_error[subdir] = {}

        for seq_id in predicted_abundance[subdir].keys():
            relative_abundance_error[subdir][seq_id] = {}

            for gr in predicted_abundance[subdir][seq_id].keys():
                relative_abundance_error[subdir][seq_id][gr] = abs(predicted_abundance[subdir][seq_id][gr] - true_abundances[seq_id][subdir])/true_abundances[seq_id][subdir] 
           

    # get average relative abundance error per subdirectory
    avg_rel_ab_error_per_subdir = {}

    for subdir in relative_abundance_error.keys():
         
        avg_rel_ab_error_seq_1 = sum(relative_abundance_error[subdir][seq_ids[0]].values())/len(relative_abundance_error[subdir][seq_ids[0]].values())
        avg_rel_ab_error_seq_2 = sum(relative_abundance_error[subdir][seq_ids[1]].values())/len(relative_abundance_error[subdir][seq_ids[1]].values())

        avg_rel_ab_error_per_subdir[subdir] = (avg_rel_ab_error_seq_1 + avg_rel_ab_error_seq_2)/2

    return avg_rel_ab_error_per_subdir

def get_estimated_haplotypes_numbers(input_dir, genomic_regions):
    n_haplotypes = {}

    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
        # make sure is a directory
        if not os.path.isdir(input_dir + '/' + subdir):
            continue
        file = input_dir + '/' + subdir + '/consensus.tsv'
        n_haplotypes[subdir] = {}
       # load the file if it exists
        if os.path.exists(file):
            consensus_tsv = pd.read_csv(file, sep='\t', header=0)
        
        for gr in consensus_tsv['Genomic_region'].unique():
            n_haplotypes[subdir][gr] = len(consensus_tsv[consensus_tsv['Genomic_region'] == gr])

    # get average number of haplotypes per subdirectory
    n_haplotypes_avg_per_subdir = {}

    for subdir in n_haplotypes.keys():
        n_haplotypes_avg_per_subdir[subdir] = sum(n_haplotypes[subdir].values())/len(n_haplotypes[subdir].values())
        
    return n_haplotypes_avg_per_subdir

def get_recall_and_duplication_ratio(input_dir, genomic_regions, sequence_ids):
    recall = {} # number of haplotypes that were correctly identified
    duplication_ratio = {} # number of haplotypes that were identified more than once
    
    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
        # make sure is a directory
        if not os.path.isdir(input_dir + '/' + subdir):
            continue
        file = input_dir + '/' + subdir + '/consensus.tsv'
        recall[subdir] = {}
        duplication_ratio[subdir] = {}

        # load the file if it exists
        if os.path.exists(file):
            consensus_tsv = pd.read_csv(file, sep='\t', header=0)

        for sequence_id in sequence_ids:
            recall[subdir][sequence_id] = {}
            for gr in genomic_regions:
                recall[subdir][sequence_id][gr] = False
                duplication_ratio[subdir][gr] = 0
        
        for index, row in consensus_tsv.iterrows():
            sequence_id = row['Sequence_id']
            genomic_region = row['Genomic_region']

            recall[subdir][sequence_id][genomic_region] = True
            duplication_ratio[subdir][genomic_region] += 1
        
    # get the average recall and duplication ratio per subdirectory
    recall_avg_per_subdir = {}
    recall_avg_per_subdir_and_seqs = {}
    duplication_ratio_avg_per_subdir = {}

    for subdir in recall.keys():
        recall_avg_per_subdir[subdir] = {}
        recall_avg_per_subdir_and_seqs[subdir] = {}
        
        true_values = [1 if value == True else 0 for value in recall[subdir][sequence_ids[0]].values()]
        recall_seq_1 = sum(true_values)/len(true_values)

        true_values = [1 if value == True else 0 for value in recall[subdir][sequence_ids[1]].values()]
        recall_seq_2 = sum(true_values)/len(true_values)

        recall_avg_per_subdir[subdir] = (recall_seq_1 + recall_seq_2)/2
        recall_avg_per_subdir_and_seqs[subdir][sequence_ids[0]] = recall_seq_1
        recall_avg_per_subdir_and_seqs[subdir][sequence_ids[1]] = recall_seq_2
        duplication_ratio_avg_per_subdir[subdir] = sum([ value/2 for value in duplication_ratio[subdir].values()])/len(duplication_ratio[subdir].values())

    return recall_avg_per_subdir, duplication_ratio_avg_per_subdir, recall_avg_per_subdir_and_seqs

def avg_homogeinity_and_completeness(input_dir, genomic_regions, sequence_ids):
    homogeinity = {}
    completeness = {}

    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
        # make sure is a directory
        if not os.path.isdir(input_dir + '/' + subdir):
            continue
        
        file = input_dir + '/' + subdir + '/community_evaluation_results.tsv'
        community_eval_tsv = pd.read_csv(file, sep='\t', header=0)
        
        homogeinity[subdir] =0
        completeness[subdir] = 0

        # load the file if it exists
        if os.path.exists(file):
            for gr in community_eval_tsv['Genomic_regions'].unique():
                homogeinity[subdir] += community_eval_tsv[community_eval_tsv['Genomic_regions'] == gr]['Homogeneity'].values[0]
                completeness[subdir] += community_eval_tsv[community_eval_tsv['Genomic_regions'] == gr]['Completeness'].values[0]
            
        homogeinity[subdir] = homogeinity[subdir]/len(community_eval_tsv['Genomic_regions'].unique())
        completeness[subdir] = completeness[subdir]/len(community_eval_tsv['Genomic_regions'].unique())

    return homogeinity, completeness

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True)
    parser.add_argument('--gt_input_dir', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--sequence_ids', type=str, required=False, help="Options: HIV-1, HCV-1b_95, HCV-1b_90")
    args = parser.parse_args()

    
    if args.sequence_ids == "HIV-1":
        genomic_regions = ["30_1028", "947_1907", "1823_2775", "2687_3621", "3516_4474", "4374_5354", "5272_6237", "6126_7120", "7029_8011", "7931_8862"]
        sequence_ids = ["OQ551959.1", "OQ551961.1"]
    elif args.sequence_ids == "HCV-1b_95":
        genomic_regions = ["72_1065", "985_1946", "1842_2800", "2703_3698", "3495_4459", "4314_5279", "5215_6167", "6068_7008", "6930_7899", "7740_8681", "8300_9280"]
        sequence_ids = ["DQ480517.1", "DQ480523.1"] 
    elif args.sequence_ids == "HCV-1b_90":
        genomic_regions = ["72_1065", "985_1946", "1842_2800", "2703_3698", "3495_4459", "4314_5279", "5215_6167", "6068_7008", "6930_7899", "7740_8681", "8300_9280"]
        sequence_ids = ["D50482.1", "D50483.1"]
    elif args.sequence_ids == "sars-cov-2": 
        genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                       (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                        (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                         (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                        (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                         (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                          (26766, 27872), (27808, 28985), (28699, 29768)]

        genomic_regions = [ str(gr[0]) + "_" + str(gr[1]) for gr in genomic_regions]
        sequence_ids = ["PP357841.1", "OR887434.1"]


   
    
    true_abundances = get_true_abundances(args.gt_input_dir, sequence_ids, genomic_regions)
    relative_abundance_error = get_relative_abundance_error(args.input_dir, sequence_ids, genomic_regions, true_abundances)
    ed_avges = get_edit_distances(args.input_dir, sequence_ids, genomic_regions)
    accuracies = get_accuracy(args.input_dir)
    n_haplotypes = get_estimated_haplotypes_numbers(args.input_dir, genomic_regions)
    recall, duplication_ratio, recall_per_subdir_and_seqs = get_recall_and_duplication_ratio(args.input_dir, genomic_regions, sequence_ids)
    homogeinity, completeness = avg_homogeinity_and_completeness(args.input_dir, genomic_regions, sequence_ids)

    # make dataframe where rows are the datasets and columns are the metrics
    df = pd.DataFrame(columns=['Dataset', 'Relative_abundance_error', 'Edit_distance', 'Accuracy', 'Number_of_haplotypes', 'Recall', 'Recall_seq1', 'Recall_seq2', 'Duplication_ratio', 'Homogeinity', 'Completeness'])
    
    for subdir in ed_avges.keys():
        # round to 3 decimal places
        df = df._append({'Dataset': subdir, 'Relative_abundance_error': round(relative_abundance_error[subdir], 3), 'Edit_distance': round(ed_avges[subdir], 2), 'Accuracy': round(accuracies[subdir], 2), 'Number_of_haplotypes': round(n_haplotypes[subdir], 2), 'Recall': round(recall[subdir], 2), 'Recall_seq1': round(recall_per_subdir_and_seqs[subdir][sequence_ids[0]], 2), 'Recall_seq2': round(recall_per_subdir_and_seqs[subdir][sequence_ids[1]], 2), 'Duplication_ratio': round(duplication_ratio[subdir], 2), 'Homogeinity': round(homogeinity[subdir], 2), 'Completeness': round(completeness[subdir], 2)}, ignore_index=True)

    # save the dataframe to a file
    df.to_csv(args.output_dir + '/summary_statistics_extensive.csv', sep=',', index=False)
    

if __name__ == "__main__":
    sys.exit(main())