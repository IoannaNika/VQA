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

                avg_predicted_abundance[sequence_id][subdir][genomic_region] = abundance
                    
            # calculate the average predicted abundance
            for seq_id in seq_ids:
                avg = 0
                for k in avg_predicted_abundance[seq_id][subdir].keys(): 
                    avg += avg_predicted_abundance[seq_id][subdir][k]
                
                average_abundance = avg / len(avg_predicted_abundance[seq_id][subdir].keys())
                avg_predicted_abundance[seq_id][subdir]['average_abundance'] = average_abundance


    return avg_predicted_abundance


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
                print(data)
                for key in data:
                    true_abundances[key] = data[key]    

    # plot a line that passes from (3, 97), (30, 70), (50, 50), (70, 30), (97, 3)
    x = [3, 30, 50, 70, 97]
    y = [3, 30, 50, 70, 97]

   # plot the line
    plt.plot(x, y, color='black', linestyle='dashed', label="True abundance")
    # genomic_regions
    genomic_regions = ["30_1028", "947_1907", "1823_2775", "2687_3621", "3516_4474", "4374_5354", "5272_6237", "6126_7120", "7029_8011", "7931_8862"]
    # plot the predicted abundances 
    avg_pred_abundance =  average_predicted_abundance(args.input_dir, true_abundances.keys(), genomic_regions)

    # plot the predicted abundances
    seq_id_1 = list(true_abundances.keys())[0]
    seq_id_2 = list(true_abundances.keys())[1]
    seq_id_1_avgs = []
    seq_id_2_avgs = []

    y_err_1 = {}
    var_1 = []
    y_err_2 = {}
    var_2 = []
    for k in["ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03"]:
        seq_id_1_avgs.append(avg_pred_abundance[seq_id_1][k]['average_abundance']*100)
        y_err_1[k] =  list(avg_pred_abundance[seq_id_1][k].values())
        y_err_1[k] = [x*100 for x in y_err_1[k]]

    for k in ["ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03"]:
        var_1.append(np.var(y_err_1[k]))

    for k in ["ab_97_03",  "ab_70_30", "ab_50_50", "ab_30_70", "ab_03_97"]:
        seq_id_2_avgs.append(avg_pred_abundance[seq_id_2][k]['average_abundance']*100)
        y_err_2[k] =  list(avg_pred_abundance[seq_id_2][k].values())
        y_err_2[k] = [x*100 for x in y_err_2[k]]

    
    for k in ["ab_97_03",  "ab_70_30", "ab_50_50", "ab_30_70", "ab_03_97"]:
        var_2.append(np.var(y_err_2[k]))
                     

    plt.scatter(x, seq_id_1_avgs, label= seq_id_1 , color = 'blue')
    plt.errorbar(x, seq_id_1_avgs, yerr=var_1, fmt='o', color='blue')

    plt.scatter(x, seq_id_2_avgs, label=seq_id_2, color = 'orange')
    plt.errorbar(x, seq_id_2_avgs, yerr=var_2, fmt='o', color='orange')

    plt.xlabel('True Abundance (%)')
    plt.yticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    plt.ylabel('Predicted Abundance (%)\n(average over genomic regions)')
    plt.title('True vs predicted abundance\nHIV-1 sequences simulated at 100x coverage')

    plt.legend()

    # save plot 
    plt.savefig(args.output_dir + '/true_vs_estimated_abundance.pdf')

        
        







    


        



if __name__ == "__main__":
    sys.exit(main())