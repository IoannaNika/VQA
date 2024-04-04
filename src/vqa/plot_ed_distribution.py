
import pandas as pd
import editdistance
import os
import matplotlib.pyplot as plt
import sys
import argparse
import json 
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--ref_set', dest = 'ref_set', required=True, type=str, help="tsv file with data")
    parser.add_argument('--directory', dest = 'directory', required=True, type=str, help="read directory")
    parser.add_argument('--triplet', dest = 'triplet', required=False, default=False, type=bool, help="triplet dataset or not")
    args = parser.parse_args()

    reference_set = pd.read_csv(args.ref_set, sep='\t', header=0)
    directory = args.directory

    pos = []
    neg = []

    for index in range(len(reference_set)):
        print(index, args.triplet)

        items = reference_set.iloc[index]

        if args.triplet == False:
            label = items["label"]
            target = 1 if label == 'positive' else 0
            read_id_1 = items["read_1"]
            read_id_2 = items["read_2"]

            
            fasta_1_path = os.path.join(directory,"reads",  '{}.fasta'.format(read_id_1))
            fasta_2_path = os.path.join(directory, "reads", '{}.fasta'.format(read_id_2))
            
            with open(fasta_1_path, 'r') as fasta_1_file:
                fasta_1 = fasta_1_file.readlines()[1].strip()
            
            
            with open(fasta_2_path, 'r') as fasta_2_file:
                fasta_2 = fasta_2_file.readlines()[1].strip()

            ed = editdistance.eval(fasta_1, fasta_2)

            if target == 1: 
                pos.append(ed)
            else:
                neg.append(ed)

        else: 
            
            if args.triplet == True: 

                read_path_anch = os.path.join(args.directory, "reads", items["read_anch"] + ".fasta")
                read_path_pos = os.path.join(args.directory, "reads", items["read_pos"] + ".fasta")
                read_path_neg = os.path.join(args.directory, "reads", items["read_neg"] + ".fasta")


                with open(read_path_anch, 'r') as read_f:
                    read_anch = read_f.readlines()[1].strip()

                with open(read_path_pos, 'r') as read_f:
                    read_pos = read_f.readlines()[1].strip()
                
                with open(read_path_neg, 'r') as read_f:
                    read_neg = read_f.readlines()[1].strip()

                ed_pos = editdistance.eval(read_anch, read_pos)
                ed_neg = editdistance.eval(read_anch, read_neg)

                pos.append(ed_pos)
                neg.append(ed_neg)

    print(len(pos), len(neg))

    with open("pos", "w") as fp:
        json.dump(pos, fp)
    with open("neg", "w") as fp:
        json.dump(neg, fp)
    
    
    f= open("pos")
    pos = json.load(f)
    print(len(pos))
    pos = np.asarray(pos)

    f= open("neg")
    neg = json.load(f)
    print(len(neg))
    neg = np.asarray(neg)


    # plot distribution of edit distances for the positive and negative samples
    # plot them on the same plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n1, bins1, _ = ax.hist(pos, bins=200, alpha=0.5, label='positive')
    n2, bins2, _ = ax.hist(neg, bins=200, alpha=0.5, label='negative')
    ax.set_xlabel('Edit distance')
    ax.set_ylabel('Frequency')
    ax.set_title('Edit distance distribution')
    ax.legend(loc='upper right')
    # save the plot
    fig.savefig('triplet_edit_distance_distribution_316.png')

    rng = min(pos.min(), neg.min()), max(pos.max(), neg.max())
    intersection = np.minimum(n1, n2)
    overlap_area = intersection.sum()
    pos_area = sum(np.diff(bins1)*n1)
    neg_area = sum(np.diff(bins2)*n2)
    total_area = pos_area + neg_area - overlap_area
    overlap_percentage = overlap_area/total_area
    print("Overlap percentage: ", overlap_percentage)


    
if __name__ == "__main__":
    sys.exit(main())