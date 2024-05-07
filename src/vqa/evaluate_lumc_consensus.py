import os
import argparse
import sys
from Bio import SeqIO
import pandas as pd
from Bio import pairwise2
from Bio.Seq import Seq
import editdistance


def get_consensus_dicts(al_consensus_wuhan, al_consensus_ba1, genomic_regions):

    consensus_wuhan_dict = {}
    consensus_ba1_dict = {}

    for i, region in enumerate(genomic_regions):
        gr = str(region[0]) + "_" + str(region[1])
        consensus_wuhan_dict[gr] = al_consensus_wuhan[region[0]:region[1]].replace("-", "").upper()
        consensus_ba1_dict[gr] = al_consensus_ba1[region[0]:region[1]].replace("-", "").upper()

    return consensus_wuhan_dict, consensus_ba1_dict

def mafft_alignment(input_temp_file, output_temp_file, wuhan, ba1):
    # clear any previous temp files
    if os.path.exists(input_temp_file):
        os.remove(input_temp_file)
    if os.path.exists(output_temp_file):
        os.remove(output_temp_file)

    with open(input_temp_file, "w") as file:
        file.write(">Wuhan\n")
        file.write(str(wuhan.seq) + "\n")
        file.write(">BA.1\n")
        file.write(str(ba1.seq))
    # run mafft
    os.system("mafft --auto --quiet --thread 4 {} > {}".format(input_temp_file, output_temp_file))
    return


def make_output_file(out_file_name):

    # delete the file if it exists
    if os.path.exists(out_file_name):
        os.remove(out_file_name)
    
    # create the file
    f = open(out_file_name, "x")
    f.write("Genomic_region\tConsensus\tRelative_abundance\tEdit_distance_Wuhan\tEdit_distance_BA1\tMin_edit_distance\tMin_consensus\tEdit_distance_between_consensus\n")
    f.close()
    return 



def main(): 
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--input_file', dest = 'input_file', required=True, type=str, help="Input file")
    parser.add_argument('--outdir', dest = 'outdir', required=True, type=str, help="output directory")
    args = parser.parse_args()

    in_temp_file = os.path.join(args.outdir, "in_temp.fasta")
    out_temp_file = os.path.join(args.outdir, "out_temp.fasta")

    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                            (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                            (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                            (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                            (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                            (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                            (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]

    consensus_wuhan_file = "data/data/lumc_consensus/consensus_nCoV-2019-Mlb-1Passage3NIETP2.fasta"
    consensus_ba1_file = "data/data/lumc_consensus/consensus_SARS-CoV-2-O-71084_2021BA1passage2.fasta"


    # read the consensus sequences

    consensus_wuhan = SeqIO.read(consensus_wuhan_file, "fasta")
    consensus_ba1 = SeqIO.read(consensus_ba1_file, "fasta")

    # align the consensus sequences
    mafft_alignment(in_temp_file, out_temp_file, consensus_wuhan, consensus_ba1)

    # read the alignment from the output file
    alignments = list(SeqIO.parse(out_temp_file, "fasta"))

    al_consensus_wuhan = alignments[0].seq
    al_consensus_ba1 = alignments[1].seq
    
    # get consensus dictionaries
    consensus_wuhan_dict, consensus_ba1_dict = get_consensus_dicts(al_consensus_wuhan, al_consensus_ba1, genomic_regions)

    # load consensus (input) file tsv 
    consensus_tsv = pd.read_csv(args.input_file, sep="\t", header=0)
    # Community	Genomic_region	Consensus	Sequence_id	Nreads	Relative_abundance

    # create output file
    out_file_name = os.path.join(args.outdir, "consensus_lumc_comparison.tsv")
    make_output_file(out_file_name)


    for index, row in consensus_tsv.iterrows():
        gr = row["Genomic_region"]
        consensus = row["Consensus"]
        sequence_id = row["Sequence_id"]

        consensus_wuhan = consensus_wuhan_dict[gr]
        consensus_ba1 = consensus_ba1_dict[gr]

        edit_distance_wuhan = editdistance.eval(consensus, consensus_wuhan)
        edit_distance_ba1 = editdistance.eval(consensus, consensus_ba1)
        edit_dist_consensus = editdistance.eval(consensus_wuhan, consensus_ba1)

        min_edit_distance = min(edit_distance_wuhan, edit_distance_ba1)
        min_consensus = "Wuhan" if edit_distance_wuhan < edit_distance_ba1 else "BA.1"
        if edit_distance_ba1 == edit_distance_wuhan:
            min_consensus = "Both"

        f = open(out_file_name, "a")
        f.write(f"{gr}\t{consensus}\t{row['Relative_abundance']}\t{edit_distance_wuhan}\t{edit_distance_ba1}\t{min_edit_distance}\t{min_consensus}\t{edit_dist_consensus}\n")
        f.close()

if __name__ == "__main__":
    sys.exit(main())