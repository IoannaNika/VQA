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
import sys
import editdistance

def read_ground_truth(gt_dir, sequence_id, genomic_region):
    # read tsv file with ground truth sequences
    gt_file = gt_dir + "/" + sequence_id + ".template"
    # open fasta template file
    with open(gt_file, "r") as file:
        lines = file.readlines()

    id_to_be_found = ">+_" + sequence_id + ":" + genomic_region + ":0"

    # find the line with the sequence
    for i in range(len(lines)):
        if lines[i].startswith(id_to_be_found):
            sequence = lines[i+1].strip()
            return sequence

def make_output_file(output_file_name):
    if os.path.exists(output_file_name):
        os.remove(output_file_name)
    
    # open file to write the consensus
    file = open(output_file_name, "w")
    file.write("Community" + "\t" + "Genomic_region" + "\t" + "Consensus" + "\t" + "Sequence_id" + "\t" + "Nreads" + "\t" + "Relative_abundance" + "\n")
    file.close()
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

def check_and_act_if_consensus_exists(output, genomic_region, consensus, results_per_community_and_gr, results):
    # if the consensus exists already in the output file then simply add the abundance of this consensus to the existing one
    # iterate through the output file and check if the consensus exists for this genomic region
    # if it exists then add the abundance of this consensus to the existing one
    # if it does not exist then write the new consensus to the output file
    output_file = pd.read_csv(output, sep='\t', header=0)
    output_file_gr = output_file[output_file["Genomic_region"] == genomic_region]
    for consensus_in_output in output_file_gr["Consensus"]:
        if editdistance.eval(consensus_in_output, consensus) == 0:
            # get the number of sequences for this consensus
            number_of_sequences = len(results_per_community_and_gr)
            relative_abundance = number_of_sequences / len(results[results["Genomic_regions"] == genomic_region])
            # add the relative abundance to the existing one                    
            output_file.loc[(output_file["Consensus"] == consensus )& (output_file["Genomic_region"] == genomic_region), "Relative_abundance"] += relative_abundance
            output_file.to_csv(output, sep='\t', index=False)
            return True  
    return False

def main():
    parser = argparse.ArgumentParser(description="Create consensus")
    parser.add_argument('--communities', dest = 'results', required=True, type=str, help="tsv file with communities")
    parser.add_argument('--gt_dir', dest = 'gt_dir', required=False, type=str, help="directory with ground truth sequences")
    parser.add_argument('--output', dest = 'output', required=True, type=str, help="output file")
    parser.add_argument('--seq_ids', dest = 'seq_ids', required=False, type=str, help="comma separated list of sequence ids")

    args = parser.parse_args()

    results = pd.read_csv(args.results, sep='\t', header=0)
    output = args.output

    file = make_output_file(output)
    sequence_ids_input = args.seq_ids.split(",")
    print(args.gt_dir)
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

            if check_and_act_if_consensus_exists(output, genomic_region, consensus, results_per_community_and_gr, results):
                continue
            
    
            gt_sequence = read_ground_truth(args.gt_dir,sequence_ids_input[0], genomic_region)
            edit_distance_1 = editdistance.eval(consensus, gt_sequence)
            gt_sequence = read_ground_truth(args.gt_dir,sequence_ids_input[1], genomic_region)
            edit_distance_2 = editdistance.eval(consensus, gt_sequence)

            if edit_distance_1 < edit_distance_2:
                sequence_id_chosen = sequence_ids_input[0]
            elif edit_distance_1 > edit_distance_2:
                sequence_id_chosen = sequence_ids_input[1]
            else:
                sequence_id_chosen = "ambiguous"

            number_of_sequences = len(results_per_community_and_gr)
            relative_abundance = number_of_sequences / len(results[results["Genomic_regions"] == genomic_region])
            if number_of_sequences >= 3: 
                with open(output, "a") as file:
                    file.write(community + "\t" + genomic_region + "\t" + consensus + "\t" + sequence_id_chosen + "\t" + str(number_of_sequences) + "\t" + str(relative_abundance) + "\n")
        
    file.close()
    clean_temp_files()

    

if __name__ == "__main__":
    sys.exit(main())