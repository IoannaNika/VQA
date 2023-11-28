

from vqa.utils import sample_examples
import os
import pandas as pd
import json
import argparse
import random

def main():
    parser = argparse.ArgumentParser(description="make triplet reference set")
    parser.add_argument('--dir', dest = 'directory', default="data/data/hcov_global_2023-11-16_09-28", required=False, type=str, help="data directory")
    parser.add_argument('--n', dest = 'n', default= 10, required=False, type=int, help="number of samples to be generated")
    args = parser.parse_args()

    n = args.n
    directory = args.directory
    # metadata file
    metadata_file = os.path.join(directory,'metadata.tsv')

    metadata_dataframe = pd.read_csv(metadata_file, sep="\t")

    data_dir = os.path.join(directory, 'split_fasta')
    lineages = os.listdir(data_dir)

    # remove samples from the metadata file that their pango_lineage is not in lineages
    metadata_dataframe = metadata_dataframe[metadata_dataframe['pango_lineage'].isin(lineages)]

    # keep only samples with lineage information
    metadata_dataframe = metadata_dataframe[metadata_dataframe['pango_lineage'] != '?']
    # keep only samples with clade information
    metadata_dataframe = metadata_dataframe[metadata_dataframe['GISAID_clade'] != '?']

    clade_lineage_info = metadata_dataframe[['GISAID_clade', 'pango_lineage']] 

    new_path = os.path.join(directory, 'triplet_reference_set')
    # create a directory to store the reference set
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    new_path = os.path.join(new_path, 'reads')
    if not os.path.exists(new_path):
        os.makedirs(new_path)
    
    # cleanup new_path
    for file in os.listdir(new_path):
        os.remove(os.path.join(new_path, file))

    reference_set_dict = {}


    # sample uniformly positive examples from each clade
    while len(reference_set_dict.keys()) < n:
        # sample a clade uniformly
        clade = clade_lineage_info['GISAID_clade'].sample(1).values[0]
        # sample a lineage uniformly from the clade
        lineage = clade_lineage_info[clade_lineage_info['GISAID_clade'] == clade]['pango_lineage'].sample(1).values[0]
        
        # sample anchor and positive
        seq_anchor, seq_positive, id_read_anchor, id_read_positive, id_anchor, id_positive = sample_examples.sample_positive(lineage)

        # sample negative
        lineage_options = [True, False]
        same_lineage = lineage_options[random.randint(0,1)]
        
        # check if the lineage directory has more than two .fasta files
        lineage_dir = os.path.join(data_dir, lineage)
        lineage_files = os.listdir(lineage_dir)
        lineage_fasta_files = [file for file in lineage_files if file.endswith('.fasta')]
        if len(lineage_fasta_files) < 2:
            same_lineage = False


        if same_lineage:
           anchor_seq, negative_seq, id_read_anchor, id_read_negative, anchor_sample_id, negative_sample_id = sample_examples.sample_single_negative(lineage,lineage, id_read_anchor)
        
        else:
            # sample clade uniformly
            clade_negative = clade_lineage_info['GISAID_clade'].sample(1).values[0]
            # sample lineage uniformly from the clade
            lineage_negative = clade_lineage_info[clade_lineage_info['GISAID_clade'] == clade_negative]['pango_lineage'].sample(1).values[0]

            while lineage_negative == lineage:
                clade_negative = clade_lineage_info['GISAID_clade'].sample(1).values[0]
                # sample lineage uniformly from the clade
                lineage_negative = clade_lineage_info[clade_lineage_info['GISAID_clade'] == clade_negative]['pango_lineage'].sample(1).values[0]
            
            anchor_seq, negative_seq, id_read_anchor, id_read_negative, anchor_sample_id, negative_sample_id = sample_examples.sample_single_negative(lineage,lineage_negative, id_read_anchor)

        # if any of the keys contains all the ids, then skip this iteration
        for key in reference_set_dict.keys():
            if id_read_anchor in key and id_read_positive in key and id_read_negative in key:
                print("sample already exists")
                continue
        # write anchor read to a file
        with open(new_path + "/{}.fasta".format(id_read_anchor), "w") as outfile:
            outfile.write(">{}\n".format(id_read_anchor))
            outfile.write("{}\n".format(anchor_seq))
        
        # write positive read to a file
        with open(new_path + "/{}.fasta".format(id_read_positive), "w") as outfile:
            outfile.write(">{}\n".format(id_read_positive))
            outfile.write("{}\n".format(seq_positive))

        # write negative to a file
        with open(new_path + "/{}.fasta".format(id_read_negative), "w") as outfile:
            outfile.write(">{}\n".format(id_read_negative))
            outfile.write("{}\n".format(negative_seq))

        reference_set_dict[id_read_anchor + "_" + id_read_positive + "_" + id_read_negative] = {
            "anchor": id_read_anchor,
            "positive": id_read_positive,
            "negative": id_read_negative
        }

    # write reference set dictionary to a json file
    with open(directory + "/triplet_reference_set/reference_set.json", "w") as outfile: 
        json.dump(reference_set_dict, outfile)

if __name__ == "__main__":
    main()

