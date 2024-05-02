import argparse
import sys
import os
from Bio import SeqIO
import editdistance
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', dest = 'directory', required=True, type=str, help="Folder with data, will be used as output folder too")
    args = parser.parse_args()

    haplotypes = [1, 5, 10, 15, 20]
    results = dict()
    for h in haplotypes:
        meta_results = os.path.join(args.directory, str(h), "meta_results_2.tsv")
        consensus_results = os.path.join(args.directory, str(h), "consensus.tsv")

        # count consensus haplotypes per genomic region
        consensus = pd.read_csv(consensus_results, sep="\t", header=0)
        
        # load meta results
        meta_results = pd.read_csv(meta_results, sep="\t", header=0)

        results[h] = dict()

        for gr in meta_results["genomic_region"].unique():
            results[h][gr] = dict()
            gr_results = meta_results[meta_results["genomic_region"] == gr]
            adjusted_abundances = gr_results["abs_relative_ab_error"].values
            editdistances = gr_results["edit_distance"].values
            n_consensus = len(gr_results)
            extra_haplotypes = n_consensus - h
            print(h, gr, extra_haplotypes, h, n_consensus)

            results[h][gr]["extra_haplotypes"] = extra_haplotypes
            results[h][gr]["mean_abs_relative_ab_error"] = sum(adjusted_abundances)/n_consensus
            results[h][gr]["mean_edit_distance"] = sum(editdistances)/n_consensus
            results[h][gr]["false_positives"]  = 0 

            for i in range(len(consensus[consensus["Genomic_region"] == gr])):
                consensus_id  = str(consensus[consensus["Genomic_region"] == gr].iloc[i]["Community"]) + "_" + consensus[consensus["Genomic_region"] == gr].iloc[i]["Genomic_region"] + "_" + str(consensus[consensus["Genomic_region"] == gr].iloc[i]["Relative_abundance"])
                if consensus_id not in gr_results["consensus_id"].values:
                    results[h][gr]["false_positives"] += 1


   

    # plot results
    # each haplotype number will be a group of bars in the plot
    extra_haplotypes = dict()
    mean_abs_relative_ab_error = dict()
    mean_edit_distance = dict()
    false_positives = dict()

    # color per haplotype
    colors = dict()

    # colordblind palette
    cbf = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
    for h in haplotypes:
        # pick a color from the palette randomly
        colors[h] = np.random.choice(cbf)
        # remove the color from the palette
        cbf.remove(colors[h])

    for h in haplotypes:
        extra_haplotypes[h] = [results[h][gr]["extra_haplotypes"] for gr in results[h]]
        mean_abs_relative_ab_error[h] = [results[h][gr]["mean_abs_relative_ab_error"] for gr in results[h]]
        mean_edit_distance[h] = [results[h][gr]["mean_edit_distance"] for gr in results[h]]
        false_positives[h] = [results[h][gr]["false_positives"] for gr in results[h]]
    
    fig, ax = plt.subplots(2, 2, figsize=(20, 18))
    bplots1 =[]
    bplots2 =[]
    bplots3 =[]
    bplots4 =[]

    for h in haplotypes:
        # each bar for each number of haplotypes should be next to each other        # plot mean and std per number of haplotypes, one bar per metric
        # ax[0].bar(h, sum(extra_haplotypes[h])/len(extra_haplotypes[h]), label=str(h))
        # overlay a box plot
        extra_haplotypes[h] = [abs(hup) for hup in extra_haplotypes[h]]

        bplot1 = ax[0, 0].boxplot(extra_haplotypes[h], positions=[h], widths=0.6, patch_artist=True, vert=True)
        bplots1.append(bplot1)

        # ax[1].bar(h, sum(mean_abs_relative_ab_error[h])/len(mean_abs_relative_ab_error[h]), label=str(h))
        bplot2 = ax[0,1].boxplot(mean_abs_relative_ab_error[h],  positions=[h], widths=0.6, patch_artist=True, vert=True)
        bplots2.append(bplot2)
        # ax[2].bar(h, sum(mean_edit_distance[h])/len(mean_edit_distance[h]), label=str(h))
        bplot3 = ax[1,0].boxplot(mean_edit_distance[h],  positions=[h], widths=0.6, patch_artist=True, vert=True)
        bplots3.append(bplot3)

        bplot4 = ax[1,1].boxplot(false_positives[h],  positions=[h], widths=0.6, patch_artist=True)
        bplots4.append(bplot4)

    for bplots in [bplots1, bplots2, bplots3, bplots4]:
        for bars,color in zip(bplots, cbf):
            for patch in bars['boxes']:
                patch.set_facecolor(color)
    

 
    ax[0,0].set_ylabel("Number of haplotypes not found", fontsize=20)
    ax[0,0].set_xlabel("True number of haplotypes",fontsize=20)   
    ax[0,0].set_ylim(bottom=0)
    ax[0,0].tick_params(axis='both', which='major', labelsize=17)  # Adjust the labelsize as needed



    ax[0,1].set_ylabel("Mean relative abundance error", fontsize=20)
    ax[0,1].set_xlabel("True number of haplotypes", fontsize=20)
    ax[0,1].set_ylim(bottom=0)
    ax[0,1].tick_params(axis='both', which='major', labelsize=17)  # Adjust the labelsize as needed


    ax[1,0].set_ylabel("Mean edit distance", fontsize=20)
    ax[1,0].set_xlabel("True number of haplotypes", fontsize=20)
    ax[1,0].set_ylim(bottom=0)
    ax[1,0].tick_params(axis='both', which='major', labelsize=17)  # Adjust the labelsize as needed


    ax[1,1].set_ylabel("Haplotype assemblies\nnot assigned to a true haplotype", fontsize=20)
    ax[1,1].set_xlabel("True number of haplotypes", fontsize=20)
    ax[1,1].set_ylim(bottom=0)
    ax[1,1].tick_params(axis='both', which='major', labelsize=17)  # Adjust the labelsize as needed


    # show only x ticks that are in the haplotypes list
    ax[0,0].set_xticks(haplotypes)
    # xlimit
    ax[0,0].set_xlim(0, 22)
    ax[0,0].set_ylim(0, 21)
    ax[0,1].set_xlim(0, 22)
    # ax[1].set_ylim(0, 20)
    ax[1,0].set_xlim(0, 22)
    # ax[2].set_ylim(0, 2)
    ax[1,1].set_xlim(0, 22)
    # ax[3].set_ylim(0, 3)

    ax[0,1].set_xticks(haplotypes)
    ax[1,0].set_xticks(haplotypes)
    ax[1,1].set_xticks(haplotypes)

    plt.tight_layout()
    plt.savefig(os.path.join(args.directory, "evaluation_results_3.pdf"))

if __name__ == "__main__":
    sys.exit(main())