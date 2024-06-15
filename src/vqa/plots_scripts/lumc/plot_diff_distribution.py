import argparse
import sys
import pandas as pd
import editdistance
import matplotlib.pyplot as plt
import os 
import pickle 


def main(): 
    parser = argparse.ArgumentParser(description="Plot the difference between edit distance of reads from consensus and contigs")
    args = parser.parse_args()
    
    # load the data
    with open("Experiments/lumc_subsample/03_50/gr_25712/diff_per_seed.pkl", "rb") as f:
        diff_per_seed = pickle.load(f)
    
    diff_per_seed_avges = {k: sum(v)/len(v) for k,v in diff_per_seed.items()}

    # plot distribution of values per seed in one plot and save it
    fig, ax = plt.subplots()
    ax.hist(diff_per_seed_avges.values(), bins=20, alpha=0.65, color='grey')
    ax.set_xlabel("Average difference between\nthe edit distance of reads from the true haplotypes\nand the reconstructed haplotypes\nper random seed")
    ax.set_ylabel("Frequency")
    # plt.xticks(range(-2, 2))
    plt.tight_layout()
    fig.savefig("src/vqa/plots_scripts/lumc/03_50/diff_per_seed_avges.pdf", format="pdf", bbox_inches="tight")

    all_values = []
    for k,v in diff_per_seed.items():
        all_values.extend(v)

    # plot distribution of values per seed in one plot and save it
    fig, ax = plt.subplots()
    #  define the bins around integers and make them line up such that the center of the bin is at an integer position
    bins  = [i - 0.5 for i in range(-5, 6)]
    # histogram
    ax.hist(all_values, bins=bins, alpha=0.65, color='grey')

   
    ax.set_xlabel("Difference between\nthe edit distance of reads from the true haplotypes\nand the reconstructed haplotypes\nacross random seeds")
    ax.set_ylabel("Frequency")
    # limit x-axis to -5 to 5
    ax.set_xlim(-4, 4)
    plt.tight_layout()
    plt.savefig("src/vqa/plots_scripts/lumc/03_50/diff_per_seed_avges_all.pdf", format="pdf", bbox_inches="tight")


if __name__ == "__main__":
    sys.exit(main())
