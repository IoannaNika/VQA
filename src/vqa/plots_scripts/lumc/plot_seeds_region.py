import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.patches as mpatches

def main():

    directory = "Experiments/lumc_subsample/03_50/gr_25712/seed_"
    seeds = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

    new_tsv = pd.DataFrame(columns=["Genomic_region",	"Consensus",	"Relative_abundance",	"Edit_distance_Wuhan",	"Edit_distance_BA1",	"Min_edit_distance",	"Min_consensus",	"Edit_distance_between_consensus"])


    for seed in seeds:
        consensus_file = directory + str(seed) + "/consensus.tsv"
        consensus_tsv = pd.read_csv(consensus_file, sep='\t', header=0)

        # write the sequences to the new tsv
        for index, row in consensus_tsv.iterrows():
            new_tsv = new_tsv._append(row, ignore_index=True)
    
    # write the new tsv to a file
    new_tsv.to_csv("Experiments/lumc_subsample/03_50/gr_25712/consensus.tsv", sep="\t", index=False)
    
    os.system("python evaluate_lumc_consensus.py --input_file Experiments/lumc_subsample/03_50/gr_25712/consensus.tsv --outdir Experiments/lumc_subsample/03_50/gr_25712")
    os.system("python plots_scripts/lumc/post_process.py --directory Experiments/lumc_subsample/03_50/gr_25712")

if __name__ == '__main__':
    main()
