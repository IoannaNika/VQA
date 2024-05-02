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
        # load meta results
        meta_results = pd.read_csv(meta_results, sep="\t", header=0)

        results[h] = dict()

        for gr in meta_results["genomic_region"].unique():
            results[h][gr] = dict()
            gr_results = meta_results[meta_results["genomic_region"] == gr]
            adjusted_abundances = gr_results["abs_relative_ab_error"].values
            editdistances = gr_results["edit_distance"].values
            n_consensus = len(gr_results)
            extra_haplotypes = n_consensus- h
            print(h, gr, extra_haplotypes, h, n_consensus)

            results[h][gr]["extra_haplotypes"] = extra_haplotypes
            results[h][gr]["mean_abs_relative_ab_error"] = sum(adjusted_abundances)/n_consensus
            results[h][gr]["mean_edit_distance"] = sum(editdistances)/n_consensus

    # plot results
    
    # each haplotype number will be a group of bars in the plot
    extra_haplotypes = dict()
    mean_abs_relative_ab_error = dict()
    mean_edit_distance = dict()

   # colordblind palette
    cbf = ["#DDCC77", '#88CCEE', '#D2A9B0', '#7DB1A8', '#D7D1F7']


    for h in haplotypes:
        extra_haplotypes[h] = [results[h][gr]["extra_haplotypes"] for gr in results[h]]
        mean_abs_relative_ab_error[h] = [results[h][gr]["mean_abs_relative_ab_error"] for gr in results[h]]
        mean_edit_distance[h] = [results[h][gr]["mean_edit_distance"] for gr in results[h]]
    
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    bplots1 =[]
    bplots2 =[]
    bplots3 =[]
    for h in haplotypes:
        # each bar for each number of haplotypes should be next to each other        # plot mean and std per number of haplotypes, one bar per metric
        # ax[0].bar(h, sum(extra_haplotypes[h])/len(extra_haplotypes[h]), label=str(h))
        # overlay a box plot
        bplot1 = ax[0].boxplot(extra_haplotypes[h], positions=[h], widths=0.6, patch_artist=True, vert=True)
        bplots1.append(bplot1)
        # ax[1].bar(h, sum(mean_abs_relative_ab_error[h])/len(mean_abs_relative_ab_error[h]), label=str(h))
        bplot2 = ax[1].boxplot(mean_abs_relative_ab_error[h],  positions=[h], widths=0.6, patch_artist=True, vert=True)
        bplots2.append(bplot2)
        # ax[2].bar(h, sum(mean_edit_distance[h])/len(mean_edit_distance[h]), label=str(h))
        bplot3 = ax[2].boxplot(mean_edit_distance[h],  positions=[h], widths=0.6, patch_artist=True, vert=True)
        bplots3.append(bplot3)

    for bplots in [bplots1, bplots2, bplots3]:
        for bars,color in zip(bplots, cbf):
            for patch in bars['boxes']:
                patch.set_facecolor(color)

    # ax[0].set_title("Mean haplotype number\nunderestimation (-) or overestimation (+)\n(over genomic regions)")
    ax[0].set_ylabel("Number of haplotypes not covered\n(over the genomic regions)")
    ax[0].set_xlabel("True number of haplotypes")           

    # ax[1].set_title("Absolute relative abundance error\n(over genomic regions)")
    ax[1].set_ylabel("Absolute relative abundance error\n(averages per genomic region)")
    ax[1].set_xlabel("True number of haplotypes")

    # ax[2].set_title("Mean edit distance\n(over genomic regions)")
    ax[2].set_ylabel("Edit distance\n(averages per genomic region)")
    ax[2].set_xlabel("True number of haplotypes")

    ax[0].set_xticks(haplotypes)
    ax[0].set_xlim(0, 22)
    ax[1].set_xticks(haplotypes)
    ax[1].set_xlim(0, 22)
    ax[2].set_xticks(haplotypes)
    ax[2].set_xlim(0, 22)


    plt.tight_layout()
    plt.savefig(os.path.join(args.directory, "evaluation_results.pdf"))

if __name__ == "__main__":
    sys.exit(main())