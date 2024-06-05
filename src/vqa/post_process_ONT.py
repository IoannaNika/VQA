import argparse
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import editdistance


def remove_short_sequences(df, min_length):
    # remove sequences that are shorter than min_length
    df = df[df["Consensus"].apply(lambda x: len(x)) >= min_length]
    return df

def adjust_relative_abundance(df):
    # adjust the relative abundance
    total_abundance = df["Relative_abundance"].sum()
    for i in range(len(df)):
        df["Relative_abundance"].iloc[i] = df["Relative_abundance"].iloc[i] / total_abundance
    return df

def remove_identical_sequences(df):
    # remove identical sequences
    sequences = df["Consensus"].values
    unique_sequences = []
    unique_indices = []
    to_delete = []
    for i, seq in enumerate(sequences):
        if seq not in unique_sequences:
            unique_sequences.append(seq)
            unique_indices.append(i)
        else:
            # if the sequence is already in the unique sequences, add the relative abundance
            index = unique_sequences.index(seq)
            print(index, i,  len(df))
            df["Relative_abundance"].iloc[index] += df["Relative_abundance"].iloc[i]
            # drop the row
            to_delete.append(df.index[i])
    
    df = df.drop(to_delete)
    
    return df
    


def get_consensus_dicts(al_consensus_wuhan, al_consensus_ba1):

    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                        (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                        (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                        (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                        (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                        (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                        (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]

    consensus_wuhan_dict = {}
    consensus_ba1_dict = {}

    for i, region in enumerate(genomic_regions):
        gr = str(region[0]) + "_" + str(region[1])
        consensus_wuhan_dict[gr] = al_consensus_wuhan[region[0]:region[1]].replace("-", "").upper()
        consensus_ba1_dict[gr] = al_consensus_ba1[region[0]:region[1]].replace("-", "").upper()

    return consensus_wuhan_dict, consensus_ba1_dict

def mafft_alignment(input_temp_file, output_temp_file, wuhan, ba1):
    # # clear any previous temp files
    # if os.path.exists(input_temp_file):
    #     os.remove(input_temp_file)
    # if os.path.exists(output_temp_file):
    #     os.remove(output_temp_file)

    with open(input_temp_file, "w") as file:
        file.write(">Wuhan\n")
        file.write(str(wuhan.seq) + "\n")
        file.write(">BA.1\n")
        file.write(str(ba1.seq))
    # run mafft
    os.system("mafft --auto --quiet --thread 4 {} > {}".format(input_temp_file, output_temp_file))
    return

def align_sequences(directory, sequences):
    # create a temporary file with the consensus sequences
    input_file = os.path.join(directory, "temp_input.fasta")
    output_file = os.path.join(directory, "temp_output.fasta")
    # make the temporary input file if it does not exist
    with open(input_file, "w") as f:
        for i, consensus_seq in enumerate(sequences):
            f.write(">{}\n".format(i))
            f.write("{}\n".format(consensus_seq))
    
    os.system("mafft --auto --quiet {} > {}".format(input_file, output_file))
    print("Alignment done")

    # read the aligned sequences
    aligned_sequences = []
    with open(output_file, "r") as f:
        lines = f.readlines()
        index = 0
        while index < len(lines):
            if lines[index].startswith(">"):
                aligned_sequences.append(lines[index+1].strip())
                index += 2
            else:
                index += 1

    
    return aligned_sequences

def do_they_only_differ_by_Ns(seq1, seq2):
    # if the sequences differ only by Ns: 
    # case 1: the sequence has an N where the other sequence has a nucleotide
    # case 2: the sequence has a nucleotide where the other sequence has an N
    # case 3: the sequence has a gap where the other sequence has an N
    # case 4: the sequence has an N where the other sequence has a gap

    for n1, n2 in zip(seq1, seq2):
        if n1 != n2 and (n1 == "N" or n2 == "N"):
            continue
        else:
            return False
    return True

def make_output_file(directory, output_file):
    # make the output file
    with open(output_file, "w") as f:
        f.write("Genomic_region\tConsensus\tRelative_abundance\tEdit_distance_Wuhan\tEdit_distance_BA1\tMin_edit_distance\tMin_consensus\tEdit_distance_between_consensus\n")

    return

def correct_Ns(seq1, seq2):
    corrected_seq = ""
    for n1, n2 in zip(seq1, seq2):            
        if n1 != n2 and n1 == "N" and n2 != "-" :
            corrected_seq += n2
        else: 
            if n1 != n2 and n1 != "-" and n2 == "N":
                corrected_seq += n1
            else:
                if n1 != n2 and n1 == "N" and n2  == "-":
                    corrected_seq += n1
                else:
    
                    if n1 != n2 and n1 == "-" and n2 == "N":
                        corrected_seq += n2
                    
                    else:
                        # if the nucleotides are the same
                        corrected_seq += n1
     
    return corrected_seq

def ed_ba1_wuhan(directory, seq, consensus_wuhan, consensus_ba1):

    temp_out = os.path.join(directory, "temp_out.fasta")
    temp_in = os.path.join(directory, "temp_in.fasta")
     # align seq with consensus to trim the ends
    with open(temp_in, "w") as f:
        f.write(">seq\n")
        f.write(seq)
        f.write("\n>consensus\n")
        f.write(str(consensus_wuhan))
    
    os.system("mafft --auto --quiet {} > {}".format(temp_in, temp_out))

    alignments = list(SeqIO.parse(temp_out, "fasta"))
    seq_aligned = alignments[0].seq
    consensus_aligned_wuhan = alignments[1].seq

    with open(temp_in, "w") as f:
        f.write(">seq\n")
        f.write(seq)
        f.write("\n>consensus\n")
        f.write(str(consensus_ba1))
  
    os.system("mafft --auto --quiet {} > {}".format(temp_in, temp_out))
   
    # read the aligned sequences
    alignments = list(SeqIO.parse(temp_out, "fasta"))
    seq_aligned = alignments[0].seq
    consensus_aligned_ba1 = alignments[1].seq
        

    
    # trim both sequences so no gaps are present
    index_start = 0
    index_end = len(seq_aligned)
    for i in range(len(seq_aligned)):
        if  consensus_aligned_ba1[i] != "-":
            index_start = i
            break

    index_end = len(seq_aligned)
    for i in range(len(seq_aligned)-1, -1, -1):
        if  consensus_aligned_ba1[i] != "-":
            index_end = i
            break


    seq_trimmed = seq_aligned[index_start:index_end+1].upper().replace("-", "")
    consensus_trimmed = consensus_aligned_wuhan[index_start:index_end+1].upper().replace("-", "")

    wh_ed = editdistance.eval(seq_trimmed, consensus_trimmed)

    index2_start = 0
    index2_end = len(seq_aligned)
    for i in range(len(seq_aligned)):
        if consensus_aligned_ba1[i] != "-":
            index2_start = i
            break
    
    
    index2_end = len(seq_aligned)
    for i in range(len(seq_aligned)-1, -1, -1):
        if consensus_aligned_ba1[i] != "-":
            index2_end = i
            break
    
    
    seq_trimmed = seq_aligned[index2_start:index2_end+1].upper().replace("-", "")
    consensus_trimmed = consensus_aligned_ba1[index2_start:index2_end+1].upper().replace("-", "")

    ed_ba1 = editdistance.eval(seq_trimmed, consensus_trimmed)

    return seq_trimmed, wh_ed, ed_ba1, seq_trimmed


    
    
    

def get_edit_distance_from_consensus(directory, omicron_dict, wuhan_dict, seq, genomic_region):

    wuhan_seq = wuhan_dict[genomic_region]
    ba1_seq = omicron_dict[genomic_region]

    seq, ed_wuhan, ed_ba1, seq_trimmed = ed_ba1_wuhan(directory, seq, wuhan_seq, ba1_seq)

    if ed_wuhan < ed_ba1:
        min_consensus = "Wuhan"
        min_ed = ed_wuhan
    else:
        if ed_wuhan > ed_ba1:
            min_consensus = "BA1"
            min_ed = ed_ba1
        else:
            min_consensus = "Both"
            min_ed = ed_wuhan


    results = {"Wuhan": ed_wuhan, "BA1": ed_ba1 , "Min": min_ed, "Min_consensus": min_consensus, "Seq": seq_trimmed}
  
    return results

def post_process(directory, df, output_file):

    # get consensus dictionary
    consensus_wuhan_dict, consensus_ba1_dict = get_conensus_regions(directory)
    
    # get consensus
    consensus_seqs = df["Consensus"].values
    # if the consensus differ only by Ns, then the consensus is the same

    # align the consensus sequences using MAFFT
    aligned_sequences = align_sequences(directory, consensus_seqs)

    index_i = 0
    index_j = 0
    while index_i < len(aligned_sequences):
        seq1 = aligned_sequences[index_i]
        index_j = index_i + 1
        while index_j < len(aligned_sequences):
            seq2 = aligned_sequences[index_j]
            if do_they_only_differ_by_Ns(seq1, seq2):
                print("Sequences {} and {} differ only by Ns".format(index_i, index_j))
                # correct the Ns
                corrected_seq = correct_Ns(seq1, seq2)
                aligned_sequences[index_i] = corrected_seq
                # modify the dataframe
                df["Consensus"].iloc[index_i] = corrected_seq
                # change the relative abundance
                df["Relative_abundance"].iloc[index_i] += df["Relative_abundance"].iloc[index_j]
                # replace the edit distances
                results = get_edit_distance_from_consensus(directory, consensus_ba1_dict, consensus_wuhan_dict, corrected_seq, df["Genomic_region"].iloc[index_i])
                df["Edit_distance_Wuhan"].iloc[index_i] = results["Wuhan"]
                df["Edit_distance_BA1"].iloc[index_i] = results["BA1"]
                df["Min_edit_distance"].iloc[index_i] = results["Min"]
                df["Min_consensus"].iloc[index_i] = results["Min_consensus"]
                # delete the second sequence
                aligned_sequences.pop(index_j)     
                # modify the dataframe
                df = df.drop(df.index[index_j])     
            else:
                index_j += 1
        
        index_i += 1

    # go through the sequences in the dataframe and get the edit distance from the consensus
    for i in range(len(df)):
        seq = df["Consensus"].iloc[i]
        genomic_region = df["Genomic_region"].iloc[i]
        results = get_edit_distance_from_consensus(directory, consensus_ba1_dict, consensus_wuhan_dict, seq, genomic_region)
        df["Edit_distance_Wuhan"].iloc[i] = results["Wuhan"]
        df["Edit_distance_BA1"].iloc[i] = results["BA1"]
        df["Min_edit_distance"].iloc[i] = results["Min"]
        df["Min_consensus"].iloc[i] = results["Min_consensus"]
        df["Consensus"].iloc[i] = results["Seq"]
        df["Edit_distance_between_consensus"].iloc[i] = editdistance.eval(consensus_wuhan_dict[genomic_region], consensus_ba1_dict[genomic_region])
    
    return  df
    

def get_conensus_regions(outdir):
    consensus_wuhan_file = "data/data/lumc_consensus/consensus_nCoV-2019-Mlb-1Passage3NIETP2.fasta"
    consensus_ba1_file = "data/data/lumc_consensus/consensus_SARS-CoV-2-O-71084_2021BA1passage2.fasta"

    # read the consensus sequences
    consensus_wuhan = SeqIO.read(consensus_wuhan_file, "fasta")
    consensus_ba1 = SeqIO.read(consensus_ba1_file, "fasta")

    in_temp_file = os.path.join(outdir, "in_temp.fasta")
    out_temp_file = os.path.join(outdir, "out_temp.fasta")
    # align the consensus sequences
    mafft_alignment(in_temp_file, out_temp_file, consensus_wuhan, consensus_ba1)
    alignments = list(SeqIO.parse(out_temp_file, "fasta"))

    al_consensus_wuhan = alignments[0].seq
    al_consensus_ba1 = alignments[1].seq

    # get consensus dictionaries
    consensus_wuhan_dict, consensus_ba1_dict = get_consensus_dicts(al_consensus_wuhan, al_consensus_ba1)

    return consensus_wuhan_dict, consensus_ba1_dict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", type=str, required=True)
    args = parser.parse_args()

    directory = args.directory

    input_file = os.path.join(directory, "consensus_lumc_comparison.tsv")
    output_file = os.path.join(directory, "consensus_lumc_comparison_post_processed.tsv")

    # make the output file
    make_output_file(directory, output_file)

    df_input = pd.read_csv(input_file, sep="\t", header=0)

    for genomic_region in df_input["Genomic_region"].unique():
        df_genomic_region = df_input[df_input["Genomic_region"] == genomic_region]
        df_genomic_region = remove_short_sequences(df_genomic_region, 800)
        # post process the data for each genomic region
        if len(df_genomic_region) > 1:
            df_genomic_region= post_process(directory, df_genomic_region, output_file)
            # leave only unique sequences
            df_genomic_region = remove_identical_sequences(df_genomic_region)
            # adjust the relative abundance
            df_genomic_region = adjust_relative_abundance(df_genomic_region)

            # write dataframe for the genomic region
            with open(output_file, "a") as f:
                for index, row in df_genomic_region.iterrows():
                        f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(row["Genomic_region"],  row["Consensus"], row["Relative_abundance"], row["Edit_distance_Wuhan"], row["Edit_distance_BA1"], row["Min_edit_distance"], row["Min_consensus"], row["Edit_distance_between_consensus"]))
        else:
            if len(df_genomic_region) != 0:
                df_genomic_region = adjust_relative_abundance(df_genomic_region)

                # write the dataframe for the genomic region
                with open(output_file, "a") as f:
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(df_genomic_region["Genomic_region"].iloc[0], df_genomic_region["Consensus"].iloc[0], df_genomic_region["Relative_abundance"].iloc[0], df_genomic_region["Edit_distance_Wuhan"].iloc[0], df_genomic_region["Edit_distance_BA1"].iloc[0], df_genomic_region["Min_edit_distance"].iloc[0], df_genomic_region["Min_consensus"].iloc[0], df_genomic_region["Edit_distance_between_consensus"].iloc[0]))





if __name__ == "__main__":
    main()