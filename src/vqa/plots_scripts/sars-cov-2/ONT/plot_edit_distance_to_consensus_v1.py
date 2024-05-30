import argparse
import json
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import seaborn as sns

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

            
            if isinstance(edit_distances[sequence_id][subdir][genomic_region], list) == False:
                edit_distances[sequence_id][subdir][genomic_region] = [edit_distance]
            else:
                edit_distances[sequence_id][subdir][genomic_region].append(edit_distance)

                
    return edit_distances

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

    
    sequence_ids = ["PP357841.1", "OR887434.1"]

    # get edit distances
    edit_distances = get_edit_distances(args.input_dir, sequence_ids, genomic_regions)

    # get average over genomic regions
    na_num = {}

    boxplots_per_subdir ={}
    boxplots_per_subdir[sequence_ids[0]] = {}
    boxplots_per_subdir[sequence_ids[1]] = {}
    

    for seq_id in sequence_ids:
        for  subdir in ['ab_03_97', 'ab_30_70', 'ab_50_50', 'ab_70_30', 'ab_97_03']:
            if subdir not in na_num.keys():
                na_num[subdir] = 0

            if subdir not in boxplots_per_subdir[seq_id].keys():
                boxplots_per_subdir[seq_id][subdir] = []

            for region in genomic_regions:
                if isinstance(edit_distances[seq_id][subdir][region], list) == False:
                   na_num[subdir] += 1
                   if subdir == "ab_97_03":
                    print(seq_id, region)
                else: 
                    boxplots_per_subdir[seq_id][subdir].extend(edit_distances[seq_id][subdir][region])

    colors = {sequence_ids[0]:"#DC267F", sequence_ids[1]:"#FE6100"}


    # plot boxplots per subdir and per sequence next to each other
    fig, ax = plt.subplots(figsize=(10, 5))
    positions_seq1 = np.arange(5) - 0.2
    positions_seq2 = np.arange(5) + 0.2
    for seq_id in sequence_ids:
        data = [boxplots_per_subdir[seq_id][subdir] for subdir in ['ab_03_97', 'ab_30_70', 'ab_50_50', 'ab_70_30', 'ab_97_03']]
        bxplot = ax.boxplot(data, positions=positions_seq1 if seq_id == sequence_ids[0] else positions_seq2, widths=0.3, patch_artist=True)
        # set the color
        for element in ['boxes']:
            plt.setp(bxplot[element], color=colors[seq_id], alpha=0.5)
        
        # set mean line to yellow
        plt.setp(bxplot['means'], color='#F9F250')
            

    xticks = []
    for subdir in ["ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03"]:
        subdir= subdir.split('_')
        if subdir[1][0] == '0':
            subdir[1] ="" + subdir[1][1:]
        if subdir[2][0] == '0':
            subdir[2] ="" + subdir[2][1:]
        subdir = subdir[1] + "% (seq1)\n" + "and\n" + subdir[2] + "% (seq2)" 
        xticks.append(subdir)  
    

    ax.set_ylabel('Edit distance to true haplotype\n(over genomic regions)', fontsize=15) #, fontweight = "bold")
    ax.set_xlabel('Relative abundance for the SARS-CoV-2 sequences', fontsize=15) #, fontweight = "bold")

    # ax.set_xticks(xticks)
    plt.legend([plt.Line2D([0], [0], color=colors[sequence_ids[0]], lw=4),
            plt.Line2D([0], [0], color=colors[sequence_ids[1]], lw=4)],
            [sequence_ids[0], sequence_ids[1]], fontsize=14)
    # plt.ylim([0, 500])

    plt.xticks(np.arange(5), xticks)
    
    plt.xticks(fontsize=13, fontweight ="normal")
    plt.yticks(fontsize=13, fontweight ="normal")


    plt.tight_layout()

    plt.savefig(args.output_dir + '/edit_distance_to_consensus_sars.pdf')

    

if __name__ == "__main__":
    sys.exit(main())