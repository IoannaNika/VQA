import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def init_editdistances(samples):   
    edit_distances = {}

    for sample in samples:
        edit_distances[sample] = []
    
    return edit_distances


def main():
    parser = argparse.ArgumentParser()
    args = parser.parse_args()
    directory = "Experiments/lumc_subsample/"

    samples = ["01_100", "02_100", "03_50", "04_75", "05_90","07_98", "08_0", "09_0"]
    edit_distances = init_editdistances(samples)

    for sample in samples:
        sample_dir = os.path.join(directory, sample)
        input_file = os.path.join(sample_dir, "consensus_lumc_comparison_post_processed.tsv")
        df = pd.read_csv(input_file, sep="\t", header=0)
        for index, row in df.iterrows():
            edit_distances[sample].append(row["Min_edit_distance"])
        
    fig = plt.figure(figsize=(10, 5))
    data = []
    
    for sample in samples:
        data.append(edit_distances[sample])

    box = plt.boxplot(data, patch_artist=True)
    colors = sns.color_palette("colorblind", 10)

    # Set colors for each boxplot
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
    
    plt.xticks(range(1, len(samples) + 1), samples, fontsize =13)
    plt.ylabel("Minimum edit distance between\ntrue and reconstructed haplotypes", fontsize=15)
    plt.yticks(fontsize=13)
    plt.xlabel("Sample name", fontsize=15)
    plt.tight_layout()
    plt.savefig("Experiments/lumc_subsample/edit_distances_per_sample.pdf")
    

if __name__ == '__main__':
    main()