import argparse
import sys
import pandas as pd
import os
import Levenshtein
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

def create_output_file(output_file_name):
    if os.path.exists(output_file_name):
        os.remove(output_file_name)
    
    # open file to write the consensus
    file = open(output_file_name, "w")
    file.write("Genomic_region" + "\t" + "Sequence_id" + "\t" + "Levenshtein_ratio" + "\t" + "Edit_distance" +"\t" + "Consensus" + "\t" + "Ground_truth" + "\n")
    return file

def main():
    parser = argparse.ArgumentParser(description="Evaluate consensus sequences")
    parser.add_argument('--consensus_seqs', dest = 'consensus_seqs', required=True, type=str, help="tsv file with consensus sequences")
    parser.add_argument('--gt_dir', dest = 'gt_dir', required=False, type=str, help="directory with ground truth sequences")
    parser.add_argument('--output', dest = 'output', required=True, type=str, help="output file")
    args = parser.parse_args()

    file = create_output_file(args.output)

    # read consensus sequences
    consensus_seqs = pd.read_csv(args.consensus_seqs, sep='\t', header=0)
    
    for i in range(len(consensus_seqs)):
        genomic_region = consensus_seqs.iloc[i]["Genomic_region"]
        consensus = consensus_seqs.iloc[i]["Consensus"]
        sequence_id = consensus_seqs.iloc[i]["Sequence_id"]

        if sequence_id == "ambiguous":
            # if the consensus is ambiguous choose randomly one of the sequences
            seqs = consensus_seqs['Sequence_id'].unique().tolist()
            seqs.remove("ambiguous")
            temp_seq_id = seqs[0]
            gt_sequence = read_ground_truth(args.gt_dir, temp_seq_id, genomic_region)
            l_ratio = Levenshtein.ratio(consensus, gt_sequence)
            l_distance = editdistance.eval(consensus, gt_sequence)
            file.write(genomic_region + "\t" + sequence_id + "\t" + str(l_ratio) + "\t" +  str(l_distance)  + "\t" + consensus + "\t" + gt_sequence + "\n")
            continue
        # read ground truth sequence
        gt_sequence = read_ground_truth(args.gt_dir, sequence_id, genomic_region)
        # calculate edit distance
        l_ratio = Levenshtein.ratio(consensus, gt_sequence)
        l_distance = editdistance.eval(consensus, gt_sequence)

        
        file.write(genomic_region + "\t" + sequence_id + "\t" + str(l_ratio) + "\t" +  str(l_distance)  + "\t" + consensus + "\t" + gt_sequence + "\n")
    

    file.close()


if __name__ == "__main__":
    sys.exit(main())