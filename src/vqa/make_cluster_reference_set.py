from vqa.utils import sample_examples
import os
import pandas as pd
import json
import random


n = 10000
# metadata file
output_data_dir = "data/data/hcov_global_2022-08-26_05-12"
data_dir = 'data/data/hcov_global_2022-08-26_05-12/split_fasta'

lineages = os.listdir(data_dir)
lineages = [lineage for lineage in lineages if lineage != "cluster_reference_set"]
print(lineages)

new_path = output_data_dir + "/" + "cluster_reference_set"
# create a directory to store the reference set
if not os.path.exists(new_path):
    os.makedirs(new_path)
else: 
    # remove all files in the directory
    for file in os.listdir(new_path):
        os.remove(new_path + "/" + file)

reference_set = set()

while len(reference_set) < n:
    # sample a clade uniformly
    # sample a lineage uniformly from the clade
    lineage = random.choice(lineages)   

    # sample a positive example from the lineage
    read1_seq, read2_seq, read1_id, read2_id, id, _ = sample_examples.sample_positive(lineage, data_dir=data_dir)
    reference_set.add(read1_id)
    reference_set.add(read2_id)
    # write read1 to a file
    with open(new_path + "/" + read1_id + ".fasta", "w") as outfile:
        outfile.write(">{}\n".format(read1_id))
        outfile.write("{}\n".format(read1_seq))
    # write read2 to a file
    with open(new_path + "/" + read2_id + ".fasta", "w") as outfile:
        outfile.write(">{}\n".format(read2_id))
        outfile.write("{}\n".format(read2_seq))


    # find a lineage with at least 2 samples
    lineage1 = random.choice(lineages)
    # check amount of files that end in .fasta in the lineage1 directory
    while len([name for name in os.listdir(data_dir + "/" + lineage1) if name.endswith(".fasta")]) < 2:
        lineage1 = random.choice(lineages)

    # sample a negative example from the lineage
    read1_seq, read2_seq, read1_id, read2_id, id, _ = sample_examples.sample_negative(True, lineage1, data_dir= data_dir)
    reference_set.add(read1_id)
    reference_set.add(read2_id)
    # write read1 to a file
    with open(new_path + "/" + read1_id + ".fasta", "w") as outfile:
        outfile.write(">{}\n".format(read1_id))
        outfile.write("{}\n".format(read1_seq))
    # write read2 to a file
    with open(new_path + "/" + read2_id + ".fasta", "w") as outfile:
        outfile.write(">{}\n".format(read2_id))
        outfile.write("{}\n".format(read2_seq))


    lineage2 = random.choice(lineages)
    while lineage2 == lineage1: 
        lineage2 = random.choice(lineages)

    # sample a negative example from the lineage
    read1_seq, read2_seq, read1_id, read2_id, id, _ = sample_examples.sample_negative(False,lineage1, lineage2, data_dir= data_dir)
    reference_set.add(read1_id)
    reference_set.add(read2_id)
    # write read1 to a file
    with open(new_path + "/" + read1_id + ".fasta", "w") as outfile:
        outfile.write(">{}\n".format(read1_id))
        outfile.write("{}\n".format(read1_seq))
    
    # write read2 to a file
    with open(new_path + "/" + read2_id + ".fasta", "w") as outfile:
        outfile.write(">{}\n".format(read2_id))
        outfile.write("{}\n".format(read2_seq))





