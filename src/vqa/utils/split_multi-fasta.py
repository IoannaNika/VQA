
import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="Split multi-fasta files into single fasta files")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="path to data directory")
    args = parser.parse_args()

    # metadata path src/vqa/data/data/hcov_global_2023-11-16_09-28
    metadata = args.dir + "/metadata.tsv"
    # sequence path
    seq_path = args.dir + "/sequences.fasta"

    # read metadata.tsv file
    metadata_df = pd.read_csv(metadata, sep='\t')

    # read sequences.fasta file
    fasta_sequences = SeqIO.parse(open(seq_path),'fasta')
    # create output directory in data directory
    parent_out_dir = args.dir + "/split_fasta"
    if not os.path.exists(parent_out_dir):
        os.mkdir(parent_out_dir)

    # total number of sequences considered
    total_seqs = 0
    # total number of sequences with no gisaid_epi_isl
    no_gisaid_epi_isl = 0
    # no pango_lineage found
    no_pango_lineage = 0

    for fasta in fasta_sequences:
        strain, sequence = fasta.id, str(fasta.seq)

        # find gisaid_epi_isl and pango_lineage from strain name
        try:
            gisaid_id = metadata_df.loc[metadata_df['strain'] == strain, 'gisaid_epi_isl'].iloc[0]
        except:
            print ("gisaid_epi_isl not found for strain: " + strain)
            no_gisaid_epi_isl += 1
            continue

        pango_lineage = metadata_df.loc[metadata_df['strain'] == strain, 'pango_lineage'].iloc[0]
        
        if pango_lineage == "?":
            print ("pango_lineage not found for strain: " + strain)
            no_pango_lineage += 1
            continue
        
        total_seqs += 1

        # place new sequence file in directory per pango lineage
        output_dir = parent_out_dir + "/" + pango_lineage
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        # create output file name
        output_file = output_dir + "/" + gisaid_id + ".fasta"
        # if file exists, delete and rewrite
        if os.path.exists(output_file):
            os.remove(output_file)
        # write to file
        with open(output_file, "w") as out_file:
            out_file.write(">" + strain + "\n" + sequence + "\n")

    print("Total number of sequences considered: " + str(total_seqs))
    print("Total number of sequences with no gisaid_epi_isl: " + str(no_gisaid_epi_isl)) 
    print("Total number of sequences with no pango_lineage: " + str(no_pango_lineage))     

if __name__ == "__main__":
    sys.exit(main())





