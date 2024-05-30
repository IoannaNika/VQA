import argparse
import json
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def get_edit_distances(input_dir, seq_ids, genomic_regions):

    edit_distances = {}

    for seq_id in seq_ids:
        edit_distances[seq_id] = {}

        for subdir in os.listdir(input_dir):
            if subdir[:2] != "ab":
                continue

            file = input_dir + '/' + subdir + '/consensus_evaluation.tsv'

            edit_distances[seq_id][subdir] = {}

            for region in genomic_regions:
                edit_distances[seq_id][subdir][region] = np.nan

    # open the file if it exists
    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
        
        file = input_dir + '/' + subdir + '/consensus_evaluation.tsv'
        
        consensus_eval_file = pd.read_csv(file, sep='\t', header=0)

        for index, row in consensus_eval_file.iterrows():
            sequence_id = row['Sequence_id']
            genomic_region = row['Genomic_region']
            edit_distance = row['Edit_distance']

            edit_distances[sequence_id][subdir][genomic_region] = edit_distance

    return edit_distances

def get_avg_edit_distances(input_dir, genomic_regions, sequence_ids):

    # get edit distances
    edit_distances = get_edit_distances(input_dir, sequence_ids, genomic_regions)

    # get average over genomic regions
    ed_avges = {}
    na_num = {}

    for seq_id in sequence_ids:
        ed_avges[seq_id] = {}
        for subdir in ["ab_003_997", "ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03", "ab_997_003"]:
            if subdir not in na_num.keys():
                na_num[subdir] = 0
            ed_avges[seq_id][subdir] = 0
            count = 0
            count_na = 0
            for region in genomic_regions:
                if not np.isnan(edit_distances[seq_id][subdir][region]):
                    ed_avges[seq_id][subdir] += edit_distances[seq_id][subdir][region]
                    count += 1
                else:
                   na_num[subdir] += 1

            ed_avges[seq_id][subdir] = ed_avges[seq_id][subdir] / count
    
    return ed_avges

def get_accuracy(input_dir): 
    accuracies = {}
    # open the file if it exists
    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
        file = input_dir + '/' + subdir + '/predictions.tsv'
        predictions_tsv = pd.read_csv(file, sep='\t', header=0)
        try:
            accuracies[subdir] =  sum(predictions_tsv['Predicted_label'] == predictions_tsv['True_label'])/len(predictions_tsv)
        except: 
            print(predictions_tsv)
    return accuracies

def average_predicted_abundance(input_dir, seq_ids, genomic_regions):
    
    predicted_abundance = {}

    for seq_id in seq_ids:
        predicted_abundance[seq_id] = {}

    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
        file = input_dir+ '/' + subdir + '/' + 'consensus.tsv'

        # predicted_abundance[seq_id][subdir] = {}

        for seq_id in seq_ids:
            predicted_abundance[seq_id][subdir] = {}
        
        for gr in genomic_regions:
            for seq_id in seq_ids:
                predicted_abundance[seq_id][subdir][gr] = 0

        # open the file if it exists
        if os.path.exists(file):
            df = pd.read_csv(file, sep='\t', header=0)
       
            for index, row in df.iterrows():
                sequence_id = row['Sequence_id']
                genomic_region = row['Genomic_region']
                abundance = row['Relative_abundance']

                predicted_abundance[sequence_id][subdir][genomic_region] = abundance
                    
    return predicted_abundance

def calculate_relative_abundance_error(true_abundances, predicted_abundance):
    relative_abundance_error = {}
    for seq_id in true_abundances.keys():
        relative_abundance_error[seq_id] = {}
        for k in predicted_abundance[seq_id].keys():
            relative_abundance_error[seq_id][k] = {}
            for gr in predicted_abundance[seq_id][k].keys():
                relative_abundance_error[seq_id][k][gr] = 0

    for seq_id in true_abundances.keys():
        for k in predicted_abundance[seq_id].keys():
            for gr in predicted_abundance[seq_id][k].keys():
                relative_abundance_error[seq_id][k][gr] = abs(predicted_abundance[seq_id][k][gr] - true_abundances[seq_id][k])/true_abundances[seq_id][k]

    return relative_abundance_error

def get_n_haplotypes(input_dir, genomic_regions):
    n_haplotypes = {}
    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
        file = input_dir + '/' + subdir + '/community_evaluation_results.tsv'
        
        for region in genomic_regions:
            if os.path.exists(file):
                predictions_tsv = pd.read_csv(file, sep='\t', header=0)
                n_haplotypes[subdir] = sum(predictions_tsv["Predicted_number_of_communities"]) / len(predictions_tsv["Predicted_number_of_communities"])

    return n_haplotypes

def get_n_haplotypes(input_dir, genomic_regions):
    n_haplotypes = {}
    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
        file = input_dir + '/' + subdir + '/consensus.tsv'
        
        for region in genomic_regions:
            if os.path.exists(file):
                consensus_tsv = pd.read_csv(file, sep='\t', header=0)
                n_haplotypes[subdir] = len(consensus_tsv) / len(consensus_tsv["Genomic_region"].unique())

    return n_haplotypes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True)
    parser.add_argument('--gt_input_dir', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    args = parser.parse_args()


    genomic_regions = ["30_1028", "947_1907", "1823_2775", "2687_3621", "3516_4474", "4374_5354", "5272_6237", "6126_7120", "7029_8011", "7931_8862"]
    sequence_ids = ["OQ551959.1", "OQ551961.1"]

    ed_avges_per_seq =  get_avg_edit_distances(args.input_dir, genomic_regions, sequence_ids)

    ed_avges = {}

    for seq_id in ed_avges_per_seq.keys():
        for dataset in ed_avges_per_seq[seq_id].keys():
            if dataset not in ed_avges.keys():
                ed_avges[dataset] = ed_avges_per_seq[seq_id][dataset]
            else:
                ed_avges[dataset] = (ed_avges[dataset] + ed_avges_per_seq[seq_id][dataset])/2
            
    
    true_abundances = {}

    # get subdirectories
    for subdir in os.listdir(args.gt_input_dir):
       mixture = subdir + '/' + 'mixture.json'

       if os.path.exists(args.gt_input_dir + '/' + mixture):
           with open(args.gt_input_dir + '/' + mixture) as f:
                data = json.load(f)
                for key in data:
                    if key not in true_abundances.keys():
                        true_abundances[key] = {}
                    true_abundances[key][subdir] = data[key]    
    
    
    avg_pred_abundance =  average_predicted_abundance(args.input_dir, true_abundances.keys(), genomic_regions)

    relative_abundance_error = calculate_relative_abundance_error(true_abundances, avg_pred_abundance)

    avg_rel_ab_error = {}

    for seq_id in relative_abundance_error.keys(): 
        for dataset in relative_abundance_error[seq_id].keys():
            if dataset not in avg_rel_ab_error.keys():
                avg_rel_ab_error[dataset] = sum(relative_abundance_error[seq_id][dataset].values())/len(relative_abundance_error[seq_id][dataset].values())
            else:
                avg_rel_ab_error[dataset] = (avg_rel_ab_error[dataset] + sum(relative_abundance_error[seq_id][dataset].values())/len(relative_abundance_error[seq_id][dataset].values()))/2

    accuracies = get_accuracy(args.input_dir)

    n_haplotypes = get_n_haplotypes(args.input_dir, genomic_regions)

    # make dataframe where rows are the datasets and columns are the metrics
    df = pd.DataFrame(columns=['Dataset', 'Average_accuracy', 'Average_edit_distance', 'Average_relative_abundance_error', 'Average_number_of_haplotypes', 'True_number_of_haplotypes'])

    for dataset in ed_avges.keys():
        df = df._append({'Dataset': dataset, 'Average_accuracy': accuracies[dataset], 'Average_edit_distance': ed_avges[dataset], 'Average_relative_abundance_error': avg_rel_ab_error[dataset], 'Average_number_of_haplotypes': n_haplotypes[dataset], 'True_number_of_haplotypes': 2}, ignore_index=True)

    # save the dataframe to a file
    df.to_csv(args.output_dir + '/summary_statistics.csv', sep=',', index=False)
    

if __name__ == "__main__":
    sys.exit(main())