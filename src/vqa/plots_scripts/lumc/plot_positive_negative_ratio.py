
import argparse
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--predictions', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    args = parser.parse_args()

    # load predictions
    predictions = pd.read_csv(args.predictions, sep="\t", header=0)

    # calculate positive and negative ratio per genomic region
    genomic_regions = predictions["Genomic_region"].unique()

    results = dict()

    for gr in genomic_regions:
        results[gr] = dict()

        gr_predictions = predictions[predictions["Genomic_region"] == gr]
        n_positives = len(gr_predictions[gr_predictions["True_label"] == 1])
        n_negatives = len(gr_predictions[gr_predictions["True_label"] == 0])
        gr.replace("_", "-")
        results[gr]["positive"] = n_positives
        results[gr]["negative"] = n_negatives


    xt = []
    for gr in results.keys():
        xt.append(gr.replace("_", "-"))

    # plot results
    x = np.arange(len(results.keys()))
    fig, ax = plt.subplots()
    negatives = [results[gr]["negative"] for gr in results.keys()]
    positives = [results[gr]["positive"] for gr in results.keys()]

    ax.bar(x+0.3, negatives, label="Negatives", color = "#E1BE6A", width = 0.3)
    ax.bar(xt, positives, label="Positives", color = "#40B0A6", width = 0.3)

    # plot legend outside of the plot
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.xlabel("Genomic region")
    plt.ylabel("Number of positive & negative samples")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "pos_neg_ratio.pdf"))
        





if __name__ == "__main__":
    sys.exit(main())


    