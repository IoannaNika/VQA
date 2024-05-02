import argparse
import json
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict


def average_predicted_abundance(input_dir, seq_ids, genomic_regions):
    
    avg_predicted_abundance = {}

    for seq_id in seq_ids:
        avg_predicted_abundance[seq_id] = {}

    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
        file = input_dir+ '/' + subdir + '/' + 'consensus.tsv'

        # avg_predicted_abundance[seq_id][subdir] = {}

        for seq_id in seq_ids:
            avg_predicted_abundance[seq_id][subdir] = {}
        
        for gr in genomic_regions:
            for seq_id in seq_ids:
                avg_predicted_abundance[seq_id][subdir][gr] = 0

        # open the file if it exists
        if os.path.exists(file):
            df = pd.read_csv(file, sep='\t', header=0)
       
            for index, row in df.iterrows():
                sequence_id = row['Sequence_id']
                genomic_region = row['Genomic_region']
                abundance = row['Relative_abundance']

                avg_predicted_abundance[sequence_id][subdir][genomic_region] = abundance
                    
            # calculate the average predicted abundance
            for seq_id in seq_ids:
                avg = 0
                for k in avg_predicted_abundance[seq_id][subdir].keys(): 
                    avg += avg_predicted_abundance[seq_id][subdir][k]
                
                average_abundance = avg / len(avg_predicted_abundance[seq_id][subdir].keys())
                avg_predicted_abundance[seq_id][subdir]['average_abundance'] = average_abundance

    return avg_predicted_abundance

def calculate_relative_abundance_error(true_abundances, avg_predicted_abundance):
    relative_abundance_error = {}
    for seq_id in true_abundances.keys():
        relative_abundance_error[seq_id] = {}
        for k in avg_predicted_abundance[seq_id].keys():
            if k == "average_abundance": 
                continue
            relative_abundance_error[seq_id][k] = {}
            for gr in avg_predicted_abundance[seq_id][k].keys():
                relative_abundance_error[seq_id][k][gr] = 0

    for seq_id in true_abundances.keys():
        for k in avg_predicted_abundance[seq_id].keys():
            if k == "average_abundance": 
                continue
            for gr in avg_predicted_abundance[seq_id][k].keys():
                relative_abundance_error[seq_id][k][gr] = abs(avg_predicted_abundance[seq_id][k][gr] - true_abundances[seq_id][k])/true_abundances[seq_id][k]

    return relative_abundance_error

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True)
    parser.add_argument('--gt_dir', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    args = parser.parse_args()


    gt_input_dir = args.gt_dir
    true_abundances = {}

    # get subdirectories
    for subdir in os.listdir(gt_input_dir):
       mixture = subdir + '/' + 'mixture.json'

       if os.path.exists(gt_input_dir + '/' + mixture):
           with open(gt_input_dir + '/' + mixture) as f:
                data = json.load(f)
                for key in data:
                    if key not in true_abundances.keys():
                        true_abundances[key] = {}
                    true_abundances[key][subdir] = data[key]    

    # plot a line that passes from (3, 97), (30, 70), (50, 50), (70, 30), (97, 3)
    x = [3, 30, 50, 70, 97]
    y = [3, 30, 50, 70, 97]

    # genomic_regions
    genomic_regions = ["72_1065", "985_1946", "1842_2800", "2703_3698", "3495_4459", "4314_5279", "5215_6167", "6068_7008", "6930_7899", "7740_8681", "8300_9280"]
    # plot the predicted abundances 
    avg_pred_abundance =  average_predicted_abundance(args.input_dir, true_abundances.keys(), genomic_regions)

    relative_abundance_error = calculate_relative_abundance_error(true_abundances, avg_pred_abundance)


    # plot the relative abundance error average per genomic region and per sequence and plot error bars

    bar_1 = []
    seq_ids = list(true_abundances.keys())
    var_1 = []


    for subdirs in ["ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03"]:
        bar_1.append(np.mean([relative_abundance_error[seq_ids[0]][subdirs][gr] for gr in genomic_regions]))
        var_1.append(np.var([relative_abundance_error[seq_ids[0]][subdirs][gr] for gr in genomic_regions]))
    


        
    bar_2 = []
    var_2 = []
    for  subdirs in["ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03"]:
        bar_2.append(np.mean([relative_abundance_error[seq_ids[1]][subdirs][gr] for gr in genomic_regions]))
        var_2.append(np.var([relative_abundance_error[seq_ids[1]][subdirs][gr] for gr in genomic_regions]))
    

   # plot the relative abundance error average per genomic region and per sequence and plot error bars 
    
    xticks = []
    for subdir in ["ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03"]:
        subdir= subdir.split('_')
        if subdir[1][0] == '0':
            subdir[1] ="" + subdir[1][1:]
        if subdir[2][0] == '0':
            subdir[2] ="" + subdir[2][1:]
        subdir = subdir[1] + "% (seq1)\n" + "and\n" + subdir[2] + "% (seq2)" 
        xticks.append(subdir)
    x = np.arange(len(bar_1))

    plt.bar(x, bar_1, yerr=var_1, capsize=3, width=0.2, label='Seq1: ' + seq_ids[0], color ="#5D3A9B")
    plt.bar(x + 0.3, bar_2, yerr=var_2, capsize=3, width=0.2, label='Seq2: ' + seq_ids[1], color ="#E66100")
  # x-ticks are the subdirectories
    plt.xticks(x, xticks, fontsize = 11)
    plt.xlabel("Abundance distribution", fontsize=13) #, fontweight='bold')
    plt.ylabel("Average relative abundance error", fontsize=13) #, fontweight='bold')
    plt.legend(fontsize=11)
    plt.ylim([0, 1.1])
    plt.tight_layout()

    plt.savefig(args.output_dir + '/' + 'relative_abundance_error_2.pdf')
    plt.close()


if __name__ == "__main__":
    sys.exit(main())