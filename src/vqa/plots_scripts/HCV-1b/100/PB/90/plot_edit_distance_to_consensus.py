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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    args = parser.parse_args()

    genomic_regions = ["72_1065", "985_1946", "1842_2800", "2703_3698", "3495_4459", "4314_5279", "5215_6167", "6068_7008", "6930_7899", "7740_8681", "8300_9280"]
    sequence_ids = ["D50482.1", "D50483.1"]

    # get edit distances
    edit_distances = get_edit_distances(args.input_dir, sequence_ids, genomic_regions)

    # get average over genomic regions
    ed_avges = {}
    na_num = {}

    for seq_id in sequence_ids:
        ed_avges[seq_id] = {}
        for  subdir in ['ab_03_97', 'ab_30_70', 'ab_50_50', 'ab_70_30', 'ab_97_03']:
            if subdir not in na_num.keys():
                na_num[subdir] = 0
            ed_avges[seq_id][subdir] = 0
            count = 0

            for region in genomic_regions:
                if not np.isnan(edit_distances[seq_id][subdir][region]):
                    ed_avges[seq_id][subdir] += edit_distances[seq_id][subdir][region]
                    count += 1
                else:
                   na_num[subdir] += 1
                   print(subdir, seq_id, region)
            try: 
                ed_avges[seq_id][subdir] = ed_avges[seq_id][subdir] / count
            except: 
                ed_avges[seq_id][subdir] = np.nan

    bar_seq_1 = []
    for subdir in ['ab_03_97', 'ab_30_70', 'ab_50_50', 'ab_70_30', 'ab_97_03']:
        bar_seq_1.append(ed_avges[sequence_ids[0]][subdir])
    
    
    bar_seq_2 = []
    for subdir in ['ab_03_97', 'ab_30_70', 'ab_50_50', 'ab_70_30', 'ab_97_03']:
            bar_seq_2.append(ed_avges[sequence_ids[1]][subdir])

    
    bar_3_nas = list(na_num.values())
    

    xticks = ['ab_03_97', 'ab_30_70', 'ab_50_50', 'ab_70_30', 'ab_97_03'] #list(edit_distances[sequence_ids[0]].keys())
    xticks = [x[3:].split("_") for x in xticks]
    xticks_new =[]

    for x in xticks: 
        if x[0][0] == "0": 
            xticks_new.append( x[0][1] + "% (seq1)" + "\nand\n" + x[1] + "% (seq2)")
        else: 
            if x[1][0] == "0":
                xticks_new.append(x[0]  + "% (seq1)" + "\nand\n" + x[1][1] + "% (seq2)")
            else: 
                xticks_new.append(x[0]  + "% (seq1)" + "\nand\n" + x[1] + "% (seq2)")

    xticks = xticks_new
    # xticks = [ if x[0][0] == "0" then x[0][1] else  x[0] + "%" + " and " + x[1] + "%" for x in xticks]
    

    x = np.arange(len(bar_seq_1))
    width = 0.2

    fig, ax = plt.subplots()
    ax.bar(x + 0.2, bar_seq_1, width, label="Seq 1: " + sequence_ids[0], color = "#DDCC77")
    ax.errorbar(x + 0.2, bar_seq_1, yerr= np.var(bar_seq_1), fmt='o', color='black', ecolor='black', elinewidth=0.2, capsize=2)
    ax.bar(x + 0.4, bar_seq_2, width, label="Seq 2: " +sequence_ids[1], color = "#332288")
    ax.errorbar(x + 0.4, bar_seq_2, yerr= np.var(bar_seq_2), fmt='o', color='black', ecolor='black', elinewidth=0.2, capsize=2)
    ax.bar(x + 0.6, bar_3_nas, width, label="Number of haplotypes not found", color = "#AA4499")
    # ax.axhline(y=11, color='black', linestyle='--', label='Number of genomic regions')


    ax.set_ylabel('Edit distance to consensus\n(average over genomic regions)', fontsize=15) #, fontweight = "bold")
    ax.set_xlabel('Relative abundance for the HCV-1b sequences', fontsize=15) #, fontweight = "bold")

    # ax.set_xticks(xticks)
    plt.legend(fontsize=12)
    # plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    xr = [i  + 0.4 for i in range(0,len(xticks))]

    plt.xticks(xr, xticks)
    
    plt.xticks(fontsize=12, fontweight ="normal")
    plt.yticks(fontsize=12, fontweight ="normal")
    plt.ylim([0, 5])

    plt.tight_layout()

    plt.savefig(args.output_dir + '/edit_distance_to_consensus_2.pdf')

    

if __name__ == "__main__":
    sys.exit(main())