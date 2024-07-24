import pandas as pd
import matplotlib.pyplot as plt
import editdistance
from Bio import SeqIO
import argparse
import os
import sys
import numpy as np


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

# create dict with genomic regions as keys and differences as values
differences_dict = {}
for gr, diff in zip(genomic_regions, differences):
    differences_dict[gr] = diff

mutation_results = "Experiments/lumc_baselines/mutations/predictions.tsv"
edit_distance_results = "Experiments/lumc_baselines/edit_distance/predictions.tsv"
goviral_results = "Experiments/lumc_baselines/goviral/predictions.tsv"

mutation_df = pd.read_csv(mutation_results, sep="\t", header=0)
edit_distance_df = pd.read_csv(edit_distance_results, sep="\t", header=0)
goviral_df = pd.read_csv(goviral_results, sep="\t", header=0)

# calculate accuracy per genomic region for goviral
goviral_accuracy = {}
edit_distance_accuracy = {}
mutation_accuracy = {}

for region in goviral_df["Genomic_region"].unique():
    region_df = goviral_df[goviral_df["Genomic_region"] == region]
    correct = region_df[region_df["Predicted_label"] == region_df["True_label"]].shape[0]
    total = region_df.shape[0]
    accuracy = correct / total
    region_key = region.split("_")[0] + "-" + region.split("_")[1]
    goviral_accuracy[region_key] = accuracy
    mutation_accuracy[region_key] = mutation_df[mutation_df["Genomic_region"] == region_key]["Accuracy"].values[0]
    edit_distance_accuracy[region_key] = edit_distance_df[edit_distance_df["Genomic_region"] == region_key]["Accuracy"].values[0]


# plot them as groups of 3 bars per genomic region
fig, ax = plt.subplots()
bar_width = 0.2
opacity = 1

regions = list(goviral_accuracy.keys())
# sort the regions based on the first number
regions.sort(key=lambda x: int(x.split("-")[0]))
goviral_accuracies = [goviral_accuracy[region] for region in regions]
mutation_accuracies = [mutation_accuracy[region] for region in regions]
edit_distance_accuracies = [edit_distance_accuracy[region] for region in regions]

index = range(len(regions))
plt.xticks([i + bar_width for i in index], regions)
# rotate the x labels
plt.xticks(rotation=90)
bar1 = ax.bar(index, goviral_accuracies, bar_width, alpha=opacity, color='#88CCEE', label='GoViral')
bar2 = ax.bar([i + bar_width for i in index], mutation_accuracies, bar_width, alpha=opacity, color='#DDCC77', label='NB-Mutation')
bar3 = ax.bar([i + 2*bar_width for i in index], edit_distance_accuracies, bar_width, alpha=opacity, color='#CC6677', label='NB-Edit distance')
ax.legend()

ax2 = ax.twinx() 
x = np.arange(len(goviral_accuracy.keys()))
ax2.plot(x, [differences_dict[key] for key in goviral_accuracy.keys()], color="purple", linestyle='', marker='o', alpha=1.0)
ax2.set_ylim([0.0, 50])
ax2.set_ylabel('Edit distance between\ntrue local haplotypes', color = "purple")

plt.xlabel('Genomic region')
ax.set_ylabel('Accuracy')
# ylim to 1.1
ax.set_ylim(0, 1.2)
plt.tight_layout()
# save the plot
plt.savefig("Experiments/lumc_baselines/accuracy_per_genomic_region.png")