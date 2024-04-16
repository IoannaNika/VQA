import argparse
import sys
import os
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="File with consensus sequences")
    
    args = parser.parse_args()
    consensus_file = os.path.join(args.dir, "consensus.tsv")
    out_file_path = os.path.join(args.dir, "consensus.fasta")

    # if it exists, delete it
    if os.path.exists(out_file_path):
        os.remove(out_file_path)
    
    # create the file
    f = open(out_file_path, "x")

    # load the consensus sequences
    consensus = pd.read_csv(consensus_file, sep="\t", header=0)

    for index, row in consensus.iterrows():
        identifier = str(row["Community"]) + "_" + row["Genomic_region"] + "_" + str(row["Relative_abundance"])
        sequence = row['Consensus']
        # open file and write the sequence
        with open(out_file_path, "a") as f:
            f.write(">" + identifier + "\n")
            f.write(sequence + "\n")
            f.close()

if __name__ == "__main__":
    sys.exit(main())
