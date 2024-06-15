import argparse
import sys
import pandas as pd
import editdistance
import matplotlib.pyplot as plt
import os 
import pickle 


def main(): 
    parser = argparse.ArgumentParser(description="")
    args = parser.parse_args()

    seeds = [str(i) for i in range(1, 50)]

    results = {}
    results_abs = {}

    for seed in seeds:
        file_path = "Experiments/lumc_subsample/03_50/gr_25712/seed_{}/consensus_lumc_comparison_post_processed.tsv".format(seed)

        df = pd.read_csv(file_path, sep="\t", header=0)

        # n haplotypes is the number of entries (rows)
        n_haplotypes = len(df)

        results[seed] = n_haplotypes
        rels_abs = list(df["Relative_abundance"])
        # make sure they sum to 1, if not, normalize
        if sum(rels_abs) != 1:
            for i in range(len(rels_abs)):
                rels_abs[i] = rels_abs[i] / sum(rels_abs)

        rels_abs = [round(i, 2) for i in rels_abs]
    
        results_abs[seed] = rels_abs
    
    # all abs
    all_abs = []
    for seed in seeds:
        all_abs.extend(results_abs[seed])

    # Plot
    plt.figure(figsize=(10, 5))
    # plot a histogram of the number of haplotypes
    # define bins form 1 to 10
    bins = [i - 0.5 for i in range(1, 11)]
    plt.hist(results.values(), bins=bins, alpha=0.65, color='grey')
    # xticks
    plt.xticks(range(1, 11, 1))
    plt.xlabel("Number of haplotypes across seeds", fontsize=17)
    # increase size of ticks 
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.ylabel("Frequency", fontsize=17)
    plt.tight_layout()
    plt.savefig("Experiments/lumc_subsample/03_50/gr_25712/n_haplotypes_across_seeds.pdf", format="pdf")

    # plot the relative abundances across seeds in a histogram
    plt.figure(figsize=(10, 5))
    # define bins from 0 to 1
    bins = [i * 0.01 for i in range(100)]
    plt.hist(all_abs, bins=bins, alpha=0.65, color='grey')
    # xticks
    plt.xticks([i * 0.1 for i in range(11)])
    plt.xlabel("Relative abundance of reconstructed haplotypes across seeds", fontsize=17)
    # increase size of ticks 
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.ylabel("Frequency", fontsize=17)
    plt.tight_layout()
    plt.savefig("Experiments/lumc_subsample/03_50/gr_25712/rel_abundance_across_seeds.pdf", format="pdf")

    
if __name__ == "__main__":
    sys.exit(main())
