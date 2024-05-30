
import argparse
import sys
import os
from Bio import SeqIO
import editdistance
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



def main():
 
    n_operations = [1,2,3,4,5,6,7,8,9,10,20,30,40,50, 60, 70, 80, 90, 100]

    insertions_dir = "Experiments/insertions"
    mutations_dir = "Experiments/mutations"
    deletions_dir = "Experiments/deletions"
    all_dir = "Experiments/indel_n_mutations"

    results_dict = {}

    for n in n_operations: 
        results_dict[n] = {}
        n_str = str(n)
        insertion_predictions = os.path.join(insertions_dir, n_str,"predictions.tsv")
        deletion_predictions = os.path.join(deletions_dir, n_str,"predictions.tsv")
        mutations_predictions = os.path.join(mutations_dir, n_str,"predictions.tsv")
        all_predictions = os.path.join(all_dir, n_str,"predictions.tsv")

        insertions_pred_tsv = pd.read_csv(insertion_predictions, sep='\t', header=0)
        deletions_pred_tsv = pd.read_csv(deletion_predictions, sep='\t', header=0)
        mutations_pred_tsv   = pd.read_csv(mutations_predictions, sep='\t', header=0)
        all_pred_tsv = pd.read_csv(all_predictions, sep='\t', header=0)


        # headers: Genomic_region	Sequence_1_id	Sequence_1	Sequence_2_id	Sequence_2	Predicted_label	Predicted_probability	True_label

        accuracy_insertions = sum(insertions_pred_tsv['Predicted_label'] == insertions_pred_tsv['True_label'])/len(insertions_pred_tsv)
        accuracy_deletions = sum(deletions_pred_tsv['Predicted_label'] == deletions_pred_tsv['True_label'])/len(deletions_pred_tsv)
        accuracy_mutations = sum(mutations_pred_tsv['Predicted_label'] == mutations_pred_tsv['True_label'])/len(mutations_pred_tsv)
        accuracy_all = sum(all_pred_tsv['Predicted_label'] == all_pred_tsv['True_label'])/len(all_pred_tsv)
        results_dict[n]["insertions"] = accuracy_insertions
        results_dict[n]["deletions"] = accuracy_deletions
        results_dict[n]["mutations"] = accuracy_mutations
        results_dict[n]['all'] = accuracy_all

    

    # plot in bar plot per n_operations groups of insertions, deletions and mutations
        
    fig, ax = plt.subplots()
    barWidth = 0.22
    r1 = np.arange(len(n_operations))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    r4 = [x + barWidth for x in r3]

    insertions = [results_dict[n]["insertions"] for n in n_operations]
    deletions = [results_dict[n]["deletions"] for n in n_operations]
    mutations = [results_dict[n]["mutations"] for n in n_operations]
    alls =  [results_dict[n]["all"] for n in n_operations]

    ax.bar(r1, insertions, color='#88CCEE', width=barWidth, edgecolor='white', label='Insertions')
    ax.bar(r2, deletions, color='#CC6677', width=barWidth, edgecolor='white', label='Deletions')
    ax.bar(r3, mutations, color='#AA4499', width=barWidth, edgecolor='white', label='Mutations')
    ax.bar(r4, alls, color='#44AA99', width=barWidth, edgecolor='white', label='All edits')


    ax.set_xlabel('Number of edits\n(insertions, deletions, substitutions)')
    ax.set_ylabel('Accuracy')
    ax.set_xticks([r + barWidth for r in range(len(n_operations))])
    ax.set_xticklabels(n_operations)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4)
    plt.subplots_adjust(top=0.9)
    plt.savefig("Experiments/mutations/accuracy_plot.pdf", bbox_inches='tight')

if __name__ == "__main__":
    sys.exit(main())