import argparse
import json
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
from matplotlib.lines import Line2D
import seaborn as sns


def average_predicted_abundance(input_dir, seq_ids, genomic_regions):
    
    avg_predicted_abundance = {}

    for seq_id in seq_ids:
        avg_predicted_abundance[seq_id] = {}

    for subdir in os.listdir(input_dir):
        if subdir[:2] != "ab":
            continue
        file = input_dir+ '/' + subdir + '/' + 'consensus.tsv'

        avg_predicted_abundance[seq_id][subdir] = {}

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

                avg_predicted_abundance[sequence_id][subdir][genomic_region] += abundance
                    
            # calculate the average predicted abundance
            for seq_id in seq_ids:
                avg = 0
                for gr in avg_predicted_abundance[seq_id][subdir].keys(): 
                    avg += avg_predicted_abundance[seq_id][subdir][gr]
                
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

    # genomic_regions
    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                       (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                        (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                         (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                        (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                         (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                          (26766, 27872), (27808, 28985), (28699, 29768)]

    genomic_regions = [ str(gr[0]) + "_" + str(gr[1]) for gr in genomic_regions]

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

    

    # PLOT BOXPLOTS per sequence next to each oter  PER SUBDIR
    # PLOT BOXPLOTS per sequence next to each oter  PER SUBDIR
    positions_seq_1 = [-0.1, 0.9, 1.9, 2.9, 3.9]
    positions_seq_2 = [0.1, 1.1, 2.1, 3.1, 4.1]

    data_seq_1 = []
    data_seq_2 = []

    fig, ax = plt.subplots(figsize=(10, 5))

    print(relative_abundance_error["OR887434.1"]["ab_97_03"]["21562_22590"])

    for subdir in ["ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03"]:
        data_seq_1.append([relative_abundance_error[seq_ids[0]][subdir][gr] for gr in genomic_regions])
        data_seq_2.append([relative_abundance_error[seq_ids[1]][subdir][gr] for gr in genomic_regions])

    bp1 = plt.boxplot(data_seq_1, positions = positions_seq_1, widths = 0.2, patch_artist=True, boxprops=dict(facecolor="#5D3A9B"),  flierprops={'marker': ".", 'markersize': 5, 'markerfacecolor': '#5D3A9B', 'markeredgecolor': '#5D3A9B', "alpha":0.5})
    bp2 = plt.boxplot(data_seq_2, positions = positions_seq_2, widths = 0.2, patch_artist=True, boxprops=dict(facecolor="#E66100"),  flierprops={'marker': ".", 'markersize': 5, 'markerfacecolor': '#E66100', 'markeredgecolor': '#E66100', "alpha":0.5})

    for element in ['means', 'medians']:
        plt.setp(bp1[element], color="#30D5C8")
        plt.setp(bp2[element], color="#30D5C8")
    # x-ticks are the subdirectories
    plt.xticks(x, xticks, fontsize = 12)
    plt.xlabel("Abundance distribution", fontsize = 13) #, fontweight='bold')
    plt.ylabel("Relative abundance error", fontsize = 13) #, fontweight='bold')
    # plot legend for the two sequences, with the colors of the boxplots with line2D etc
    custom_lines = [Line2D([0], [0], color="#5D3A9B", lw=4),
                    Line2D([0], [0], color="#E66100", lw=4)]
    
    # plt.legend(custom_lines, [seq_ids[0], seq_ids[1]], fontsize = 12)
    plt.legend(custom_lines, ["Seq 1: " + seq_ids[0], "Seq 2: " +  seq_ids[1]], loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2, fontsize = 12)

    
    # plt.ylim([0, 1.1])
    plt.tight_layout()

    plt.savefig(args.output_dir + '/' + 'relative_abundance_error_sars2_ont.pdf')
    plt.close()


if __name__ == "__main__":
    sys.exit(main())