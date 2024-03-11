import argparse
import sys
import pandas as pd
import igraph as ig
import random
import os
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def make_output_file(output_file_name):
    if os.path.exists(output_file_name):
        os.remove(output_file_name)
    
    # open file to write the consensus
    file = open(output_file_name, "w")
    file.write("Community" + "\t" + "Genomic_region" + "\t" + "Consensus" + "\t" + "Sequence_id" + "\n")
    return file

def clean_temp_files():
    if os.path.exists("temp_input.fasta"):
        os.remove("temp_input.fasta")
    if os.path.exists("temp_output.fasta"):
        os.remove("temp_output.fasta")
    return

def mafft_alignment(sequences):
    # clear any previous temp files
    if os.path.exists("temp_input.fasta"):
        os.remove("temp_input.fasta")
    if os.path.exists("temp_output.fasta"):
        os.remove("temp_output.fasta")

    with open("temp_input.fasta", "w") as file:
        for i, seq in enumerate(sequences):
            file.write(">" + str(i) + "\n")
            file.write(seq + "\n")
    # run mafft
    os.system("mafft --auto --quiet --thread 4 temp_input.fasta > temp_output.fasta")
    return


def create_consensus(community, genomic_region, results):
    # get consensus sequence
    sequences = []
    sequence_ids = []
    # get rows that correspond to the community and genomic region
    community_results = results[results["Community"] == community]
    community_results = community_results[community_results["Genomic_regions"] == genomic_region]

    # get sequences 
    for i in range(len(community_results)):
        sequences.append(community_results.iloc[i]["Sequence"])
        sequence_ids.append(community_results.iloc[i]["Sequence_id"])

    print("Number of sequences: ", len(sequences))
    print("Mafft alignment to be performed...")
    mafft_alignment(sequences)
    print("Mafft alignment done.")
    alignment = AlignIO.read("temp_output.fasta", "fasta")  
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.gap_consensus(ambiguous='N', threshold=0.5)
    # to upper case
    consensus = consensus.upper()
    # # get string representation of the consensus
    consensus = str(consensus).replace("-", "")
    return consensus

def get_results_per_community_and_genomic_region(community, genomic_region, results):
    results_per_community = results[results["Community"] == int(community)]
    results_per_community_and_gr = results_per_community[results_per_community["Genomic_regions"] == genomic_region]
    return results_per_community_and_gr

def main():
    parser = argparse.ArgumentParser(description="Create consensus")
    parser.add_argument('--communities', dest = 'results', required=True, type=str, help="tsv file with communities")
    parser.add_argument('--output', dest = 'output', required=True, type=str, help="output file")
    args = parser.parse_args()

    results = pd.read_csv(args.results, sep='\t', header=0)
    output = args.output

    file = make_output_file(output)

    # get genomic regions
    genomic_regions = results["Genomic_regions"].unique()

    # get communities per genomic region
    for genomic_region in genomic_regions:
        communities = results[results["Genomic_regions"] == genomic_region]
        communities = communities["Community"].unique()

        for community in communities:
            consensus = create_consensus(community, genomic_region, results)
            community = str(community)
            results_per_community_and_gr = get_results_per_community_and_genomic_region(community, genomic_region, results)
            try: 
                sequence_ids = results_per_community_and_gr["Sequence_id"].tolist()
                print(sequence_ids)
                final_sequence_ids = [ seq_id.split("_")[0] for seq_id in sequence_ids]
                # get the most frequent sequence id
                sequence_id = max(set(final_sequence_ids), key=final_sequence_ids.count)
                # write to file
            except:
                sequence_id = "NA"
            file.write(community + "\t" + genomic_region + "\t" + consensus + "\t" + sequence_id + "\n")
        
    
    file.close()
    clean_temp_files()

if __name__ == "__main__":
    sys.exit(main())