import argparse
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import editdistance
from Bio import SeqIO


def parse_nb(file): 

    with open(file, "r") as f:
        lines = f.readlines()

        nb_dict =  {}

        # Genomic region: 54 - 1183
        # Accuracy:  0.5315068493150685

        cnt = 0
        while cnt < len(lines):
            line = lines[cnt]
            gr = line.split(":")[1].strip().split("-")
            gr = gr[0].strip() + "_" + gr[1].strip()
            cnt += 1
            line = lines[cnt]
            accuracy = line.split(":")[1].strip()
            if accuracy == "NA": 
                accuracy = 0.0
            else:
                accuracy = float(accuracy)
            nb_dict[gr] = accuracy
            cnt += 1
    
    return nb_dict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    args = parser.parse_args()

    file = os.path.join(args.input_dir, "predictions.tsv")
    nb_predictions = os.path.join(args.input_dir, "nb_predictions.out")

    nb_dict = parse_nb(nb_predictions)



    # read tsv file
    df = pd.read_csv(file, sep='\t', header=0)

    results = {}


    for gr in df['Genomic_region'].unique():
        results[gr] = {}
        results[gr]['accuracy'] = 0
        results[gr]['precision'] = 0
        results[gr]['recall'] = 0
        # calculate accuracy with Predicted_label and True_label per Genomic_region
        for index, row in df.iterrows():
            if row['Genomic_region'] == gr:
                if row['Predicted_label'] == row['True_label']:
                    results[gr]['accuracy'] += 1
                if row['Predicted_label'] == 1 and row['True_label'] == 1:
                    results[gr]['precision'] += 1
                    results[gr]['recall'] += 1
            
        results[gr]['accuracy']  = results[gr]['accuracy'] / len(df[df['Genomic_region'] == gr])
        results[gr]['precision']  = results[gr]['precision'] / len(df[(df['Genomic_region'] == gr) & (df['Predicted_label'] == 1)])
        results[gr]['recall'] = results[gr]['recall'] / len(df[(df['Genomic_region']== gr) & (df['True_label'] == 1)])


    wuhan = "data/data/lumc_consensus/consensus_nCoV-2019-Mlb-1Passage3NIETP2.fasta"
    omicron = "data/data/lumc_consensus/consensus_SARS-CoV-2-O-71084_2021BA1passage2.fasta"

    wuhan_seq = SeqIO.read(wuhan, "fasta")
    omicron_seq = SeqIO.read(omicron, "fasta")

    # put the two sequences in a file named "wuhan_omicron.fasta"
    wuhan_omicron = "data/data/lumc_consensus/wuhan_omicron.fasta"
    SeqIO.write([wuhan_seq, omicron_seq], wuhan_omicron, "fasta")

    # align the two sequences using mafft
    os.system("mafft --auto --quiet --thread 4 data/data/lumc_consensus/wuhan_omicron.fasta > data/data/lumc_consensus/msa_alignment.fasta")

    # parse multiple sequence alignment
    alignment = SeqIO.parse("data/data/lumc_consensus/msa_alignment.fasta", "fasta")

    # find the differences between the two sequences
    wuhan_seq = next(alignment).seq
    omicron_seq = next(alignment).seq

    bins = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                                (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                                (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                                (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                                (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                                (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                                (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]

    # find amount of differences in each bin
    differences = []
    genomic_regions = []
    for bin in bins:
        wuhan_bin = wuhan_seq[bin[0]:bin[1]]
        omicron_bin = omicron_seq[bin[0]:bin[1]]
        differences.append(editdistance.eval(str(wuhan_bin), str(omicron_bin)))
        genomic_regions.append(str(bin[0]) + "-" + str(bin[1]))

    m_diff = max(differences)

    # create dict with genomic regions as keys and differences as values
    differences_dict = {}
    for gr, diff in zip(genomic_regions, differences):
        differences_dict[gr] = diff


    # plot barplot with accuracy, precision and recall per Genomic_region
    
    bar_accuracy = []

    for gr in results.keys():
        bar_accuracy.append(results[gr]['accuracy'])
    
    bar_precision = []

    for gr in results.keys():
        bar_precision.append(results[gr]['precision'])

    bar_recall = []

    for gr in results.keys():
        bar_recall.append(results[gr]['recall'])

    fig, ax = plt.subplots()
    #plot bars next to each other
    width = 0.25

    x = np.arange(len(results.keys()))

    xticks = []
    for gr in results.keys():
        xticks.append(gr.replace('_', '-'))
        
    plt.xticks(x, xticks, rotation=90, fontsize=9)
    plt.legend(fontsize=9)

    # ax.bar(x, bar_accuracy, width, label='GoViral model', color = "deepskyblue")
    # ax.bar(x + 0.3, [nb_dict[key] for key in results.keys()], width, label='Naive Bayes model (baseline)', color= "#FE6100", alpha=0.7)
    # ax.set_ylim([0.0, 1.1])
    # ax.set_ylabel("Accuracy")
    # ax.set_xlabel('Genomic region')
    # ax.legend()

    ax.bar(x, bar_recall, width, label='Recall', color = "#56B4E9")
    ax.bar(x + 0.3, bar_precision, width, label='Precision', color= "#CC79A7")
    ax.set_ylim([0.0, 1.3])
    ax.set_ylabel("Score")
    ax.set_xlabel('Genomic region')
    ax.legend()

    ax2 = ax.twinx() 
    # plt.bar(x+0.2, bar_precision, width, label='Precision', color = "#FFB000")
    # plt.bar(x+0.4, bar_recall, width, label='Recall', color = "#648FFF")
    ax2.plot(x, [differences_dict[key.replace("_", "-")] for key in results.keys()], color="purple", linestyle='', marker='o', alpha=1.0)
    ax2.set_ylim([0.0, 50])
    ax2.set_ylabel('Edit distance between\ntrue local haplotypes', color = "purple")

    # rotate x-axis labels

    plt.ylim(bottom=0)
    # plt.yticks(fontsize=12)

    plt.tight_layout()

    plt.savefig(os.path.join(args.output_dir, "accuracy_predictions_goViral_recall_precision.pdf"))
        

if __name__ == "__main__":
    sys.exit(main())