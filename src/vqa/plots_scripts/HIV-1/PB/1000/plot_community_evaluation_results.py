import argparse
import json
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    args = parser.parse_args()

    results = {}

    # get all community evaluation results in the input directory
    for sub_directory in os.listdir(args.input_dir):
        input_file = os.path.join(args.input_dir, sub_directory, "community_evaluation_results.tsv")
        consensus_file = os.path.join(args.input_dir, sub_directory, "consensus.tsv")

        if os.path.exists(consensus_file):
            consensus_df = pd.read_csv(consensus_file, sep='\t')

        if os.path.exists(input_file):
            df = pd.read_csv(input_file, sep='\t')

            homogeneity = df['Homogeneity'].values
            completeness = df['Completeness'].values
            n_pred_communities = df['Predicted_number_of_communities'].values
            n_true_communities = df['True_number_of_communities'].values
            n_pred_edges = df['Predicted_number_of_edges'].values
            n_true_edges = df['True_number_of_edges'].values
            regions = df['Genomic_regions'].values

            results[sub_directory] = {}
            for region in regions: 
                results[sub_directory][region] = {}
                results[sub_directory][region]['homogeneity'] = homogeneity[regions == region][0]
                results[sub_directory][region]['completeness'] = completeness[regions == region][0]
                consensus_gr_df = consensus_df[consensus_df['Genomic_region'] ==region]
                n_pred_communities = len(consensus_gr_df[consensus_gr_df['Nreads'] >= 3])
                results[sub_directory][region]['n_pred_communities'] = n_pred_communities #n_pred_communities[regions == region][0]                results[sub_directory][region]['n_true_communities'] = n_true_communities[regions == region][0]
                results[sub_directory][region]['n_true_communities'] = n_true_communities[regions == region][0]
                results[sub_directory][region]['n_pred_edges'] = n_pred_edges[regions == region][0]
                results[sub_directory][region]['n_true_edges'] = n_true_edges[regions == region][0]
            
    avg_results = {}
    # plot homogeneity and completeness for subdirectories and genomic regions in one plot as bar plot groups
    for  subdir in results.keys():
        avg_results[subdir] = {'homogeneity': 0, 'completeness': 0, 'n_pred_communities': 0, 'n_true_communities': 0, 'n_pred_edges': 0, 'n_true_edges': 0}
        for gr in results[subdir].keys():
            homogeneity = results[subdir][gr]['homogeneity']
            completeness = results[subdir][gr]['completeness']
            n_pred_communities = results[subdir][gr]['n_pred_communities']
            n_true_communities = results[subdir][gr]['n_true_communities']
            n_pred_edges = results[subdir][gr]['n_pred_edges']
            n_true_edges = results[subdir][gr]['n_true_edges']

            avg_results[subdir]['homogeneity'] += homogeneity
            avg_results[subdir]['completeness'] += completeness
            avg_results[subdir]['n_pred_communities'] += n_pred_communities
            avg_results[subdir]['n_true_communities'] += n_true_communities
            avg_results[subdir]['n_pred_edges'] += n_pred_edges
            avg_results[subdir]['n_true_edges'] += n_true_edges
        
        avg_results[subdir]['homogeneity'] = avg_results[subdir]['homogeneity']/len(results[subdir].keys())
        avg_results[subdir]['completeness'] = avg_results[subdir]['completeness']/len(results[subdir].keys())
        avg_results[subdir]['n_pred_communities'] = avg_results[subdir]['n_pred_communities']/len(results[subdir].keys())
        avg_results[subdir]['n_true_communities'] = avg_results[subdir]['n_true_communities']/len(results[subdir].keys())
        avg_results[subdir]['ratio_communities'] = avg_results[subdir]['n_pred_communities']/avg_results[subdir]['n_true_communities']
        avg_results[subdir]['n_pred_edges'] = avg_results[subdir]['n_pred_edges']/len(results[subdir].keys())
        avg_results[subdir]['n_true_edges'] = avg_results[subdir]['n_true_edges']/len(results[subdir].keys())
        avg_results[subdir]['ratio_edges'] = avg_results[subdir]['n_pred_edges']/avg_results[subdir]['n_true_edges']

    # plot homogeneity and completeness for subdirectories and genomic regions in one plot as bar plot groups
    fig = plt.figure(figsize=(20, 10))


    # bar with all homogeneity values
    bars_homogeneity = [avg_results[subdir]['homogeneity'] for subdir in ["ab_003_997", "ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03", "ab_997_003"] ]
    #  variance of homogeneity values
    variance_homogeneity = []
    for subdir in ["ab_003_997", "ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03", "ab_997_003"]:
            variance = np.var([results[subdir][gr]['homogeneity'] for gr in results[subdir].keys()])
            variance_homogeneity.append(variance)
    
    # bar with all completeness values
    bars_completeness = [avg_results[subdir]['completeness'] for subdir in ["ab_003_997", "ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03", "ab_997_003"] ]
    # variance of completeness values
    variance_completeness = []
    for subdir in ["ab_003_997", "ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03", "ab_997_003"]:
            variance = np.var([results[subdir][gr]['completeness'] for gr in results[subdir].keys()])
            variance_completeness.append(variance)

    # bar with all ratio of predicted communities
    bars_ratio_communities = [avg_results[subdir]['ratio_communities'] for subdir in ["ab_003_997", "ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03", "ab_997_003"] ]

    # variance of ratio of predicted communities
    variance_ratio_communities = []
    for subdir in ["ab_003_997", "ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03", "ab_997_003"]:
            variance = np.var([results[subdir][gr]['n_pred_communities']/results[subdir][gr]['n_true_communities'] for gr in results[subdir].keys()])
            variance_ratio_communities.append(variance) 

    # bar with all ratio of predicted edges
    bars_ratio_edges = [avg_results[subdir]['ratio_edges'] for subdir in avg_results.keys()]

    # variance of ratio of predicted edges
    variance_ratio_edges = []
    for subdir in avg_results.keys():
            variance = np.var([results[subdir][gr]['n_pred_edges']/results[subdir][gr]['n_true_edges'] for gr in results[subdir].keys()])
            variance_ratio_edges.append(variance)
    plt.axhline(y=1, color='black', linestyle='--')
    plt.bar(np.arange(len(avg_results.keys())), bars_homogeneity, color='#882255', width=0.2, label='Homogeneity')
    plt.errorbar(np.arange(len(avg_results.keys())), bars_homogeneity, yerr=variance_homogeneity, fmt='o', color='black', capsize=5)
    plt.bar(np.arange(len(avg_results.keys())) + 0.2, bars_completeness, color='#44AA99', width=0.2, label='Completeness')
    plt.errorbar(np.arange(len(avg_results.keys())) + 0.2, bars_completeness, yerr=variance_completeness, fmt='o', color='black', capsize=5)
    plt.bar(np.arange(len(avg_results.keys())) + 0.4, bars_ratio_communities, color='#DDCC77', width=0.2, label='Ratio of predicted communities to true communities')
    plt.errorbar(np.arange(len(avg_results.keys())) + 0.4, bars_ratio_communities, yerr=variance_ratio_communities, fmt='o', color='black', capsize=5)
    # plt.bar(np.arange(len(avg_results.keys())) + 0.6, bars_ratio_edges, color='#88CCEE', width=0.2, label='Ratio of predicted edges\nto true edges')
    # plt.errorbar(np.arange(len(avg_results.keys())) + 0.6, bars_ratio_edges, yerr=variance_ratio_edges, fmt='o', color='black', capsize=5)

    # limit results for y-axis above 0.5
    plt.ylim(0.5, 1.25)
    # make x-ticks more readable
    plt.yticks(fontsize=30)

    xticks = []
    for subdir in ["ab_003_997", "ab_03_97", "ab_30_70", "ab_50_50", "ab_70_30", "ab_97_03", "ab_997_003"] :
        subdir= subdir.split('_')
        if subdir[1][0] == '0':
            subdir[1] ="" + subdir[1][1:]
        if subdir[2][0] == '0':
            subdir[2] ="" + subdir[2][1:]
        subdir = subdir[1] + "% " + "and " + subdir[2] + "%"
        xticks.append(subdir)

        
    plt.xticks(np.arange(len(avg_results.keys())) + 0.3, xticks, fontsize=30)
    # put legend outside of the plot in the middle of the right side
    plt.legend(fontsize=30)
    plt.ylabel('Average scores\nover genomic regions', fontsize=35) #, fontweight ="bold")
    plt.xlabel('Relative abundances for the two simulated HIV-1 sequences', fontsize=35) #, fontweight ="bold")
    # plt.title('HIV-1: Community evaluation results (100x), testing for relative abundance', fontsize=20, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "community_evaluation_results_v2.pdf"))
    plt.close(fig)
    





if __name__ == "__main__":
    sys.exit(main())

      
    



    
            





 
       


            

            

        




if __name__ == "__main__":
    sys.exit(main())