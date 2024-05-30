import argparse
import sys
import pandas as pd
import editdistance
import matplotlib.pyplot as plt
import os
from Bio import SeqIO   
import pickle 


def mafft_alignment_wuhan_ba1(input_temp_file, output_temp_file, wuhan, ba1):
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

def get_consensus_dicts(al_consensus_wuhan, al_consensus_ba1, genomic_regions):

    consensus_wuhan_dict = {}
    consensus_ba1_dict = {}

    for i, region in enumerate(genomic_regions):
        gr = str(region[0]) + "_" + str(region[1])
        consensus_wuhan_dict[gr] = al_consensus_wuhan[region[0]:region[1]].replace("-", "").upper()
        consensus_ba1_dict[gr] = al_consensus_ba1[region[0]:region[1]].replace("-", "").upper()

    return consensus_wuhan_dict, consensus_ba1_dict

def align_read_to_contig(contig, read, gr):
   # align read using mafft
    mafft_alignment(contig, read, gr)
    # read the aligned sequences from the output file, there are two sequences in the file
    seqs = list(SeqIO.parse("temp_output_{}.fasta".format(gr), "fasta"))
    contig_aligned = str(seqs[0].seq)
    read_aligned = str(seqs[1].seq)
    
    start = 0
    end = 0

    # find the start and end of the read in the consensus
    for i in range(len(contig_aligned)):
        if contig_aligned[i] != '-' and read_aligned[i] != '-':
            start = i
            break
    for i in range(len(contig_aligned)-1, -1, -1):
        if contig_aligned[i] != '-' and read_aligned[i] != '-':
            end = i
            break
    
    contig_trimmed = contig_aligned[start:end+1].replace("-", "").upper()
    read_trimmed = read_aligned[start:end+1].replace("-", "").upper()

    ed = editdistance.eval(contig_trimmed, read_trimmed)
    
    return ed

def get_wuhan_omicron_contigs():
      
    in_temp_file = "temp_input_wuhan_ba1.fasta"
    out_temp_file = "temp_output_wuhan_ba1.fasta"
    
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
    mafft_alignment_wuhan_ba1(in_temp_file, out_temp_file, consensus_wuhan, consensus_ba1)

    # read the alignment from the output file
    alignments = list(SeqIO.parse(out_temp_file, "fasta"))

    al_consensus_wuhan = str(alignments[0].seq)
    al_consensus_ba1 = str(alignments[1].seq)
    
    # get consensus dictionaries
    consensus_wuhan_dict, consensus_ba1_dict = get_consensus_dicts(al_consensus_wuhan, al_consensus_ba1, genomic_regions)

    return consensus_wuhan_dict, consensus_ba1_dict

def mafft_alignment(consensus, read, gr):
    # # clear any previous temp files
    if os.path.exists("temp_input_{}.fasta".format(gr)):
        os.remove("temp_input_{}.fasta".format(gr))
    if os.path.exists("temp_output_{}.fasta".format(gr)):
        os.remove("temp_output_{}.fasta".format(gr))

    with open("temp_input_{}.fasta".format(gr), "w") as file:
        for i, seq in enumerate([consensus, read]):
            file.write(">" + str(i) + "\n")
            file.write(seq + "\n")
    # run mafft
    os.system("mafft --auto --quiet --thread 4 temp_input_{}.fasta > temp_output_{}.fasta".format(gr, gr))
    return


def main(): 
    parser = argparse.ArgumentParser(description="This script will calculate the consensus distance from reads")
    args = parser.parse_args()

    edit_distances_per_sample_consensus = {}
    edit_distances_per_sample_contigs = {}

    samples = ["01_100", "02_100", "03_50", "04_75", "05_90", "06_95", "07_98", "08_0", "09_0"]
    consensus_wuhan_dict, consensus_ba1_dict = get_wuhan_omicron_contigs()

    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                            (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                            (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                            (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                            (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                            (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                            (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]
    
    genomic_regions = [str(region[0]) + "_" + str(region[1]) for region in genomic_regions]

    for sample in samples:
        consensus_file = "Experiments/lumc_subsample/" + sample + '/consensus.tsv'
        reads_file = "Experiments/lumc_subsample/" + sample + '/communities.tsv'

        consensus = pd.read_csv(consensus_file, sep='\t', header=0)
        reads = pd.read_csv(reads_file, sep='\t', header=0)

        edit_distances_per_sample_consensus[sample] = []
        edit_distances_per_sample_contigs[sample] = []
        
        for gr in genomic_regions:
            wuhan_gr = consensus_wuhan_dict[gr]
            ba1_gr = consensus_ba1_dict[gr]
            # get reads for this region
            reads_region = reads[reads["Genomic_regions"] == gr]
            contigs_region = consensus[consensus["Genomic_region"] == gr]

            for i, read in reads_region.iterrows():
                read_seq = read["Sequence"]
                ed_wuhan = align_read_to_contig(wuhan_gr, read_seq, gr)
                ed_ba1 = align_read_to_contig(ba1_gr, read_seq, gr)
                edit_distances_per_sample_consensus[sample].append(min(ed_wuhan, ed_ba1))

                min_contig_ed = float("inf")
                for j, contig in contigs_region.iterrows():
                    contig_seq = contig["Consensus"]
                    ed = align_read_to_contig(contig_seq, read_seq, gr)
                    if ed < min_contig_ed:
                        min_contig_ed = ed
                
                edit_distances_per_sample_contigs[sample].append(min_contig_ed)
    
    f_consensus = open("Experiments/lumc_subsample/edit_distances_per_sample_consensus.pkl","wb")
    f_contigs = open("Experiments/lumc_subsample/edit_distances_per_sample_contigs.pkl", "wb")

    pickle.dump(edit_distances_per_sample_consensus, f_consensus)
    f_consensus.close()
    pickle.dump(edit_distances_per_sample_contigs, f_contigs)
    f_contigs.close()

    for gr in genomic_regions:
        # cleanup temp files
        if  os.path.exists("temp_input_{}.fasta".format(gr)):
            os.remove("temp_input_{}.fasta".format(gr))
        
        if  os.path.exists("temp_output_{}.fasta".format(gr)):
            os.remove("temp_output_{}.fasta".format(gr))

            
if __name__ == "__main__":
    sys.exit(main())
