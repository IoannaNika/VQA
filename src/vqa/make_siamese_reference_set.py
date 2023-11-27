from vqa.utils import sample_examples
import os
import pandas as pd
import json


n = 100
# metadata file
metadata_file = 'data/data/hcov_global_2023-11-16_09-28/metadata.tsv'

metadata_dataframe = pd.read_csv(metadata_file, sep="\t")

data_dir = 'data/data/hcov_global_2023-11-16_09-28/split_fasta'
lineages = os.listdir(data_dir)

# remove samples from the metadata file that their pango_lineage is not in lineages
metadata_dataframe = metadata_dataframe[metadata_dataframe['pango_lineage'].isin(lineages)]

# keep only samples with lineage information
metadata_dataframe = metadata_dataframe[metadata_dataframe['pango_lineage'] != '?']
# keep only samples with clade information
metadata_dataframe = metadata_dataframe[metadata_dataframe['GISAID_clade'] != '?']

clade_lineage_info = metadata_dataframe[['GISAID_clade', 'pango_lineage']] 

counts = clade_lineage_info.value_counts()

# set counts to zero 
ref_set_stats = counts.reset_index(name='count')
ref_set_stats['count'] = 0

# create a directory to store the reference set
if not os.path.exists('data/data/hcov_global_2023-11-16_09-28/siamese_reference_set'):
    os.makedirs('data/data/hcov_global_2023-11-16_09-28/siamese_reference_set')

if not os.path.exists('data/data/hcov_global_2023-11-16_09-28/siamese_reference_set/reads'):
    os.makedirs('data/data/hcov_global_2023-11-16_09-28/siamese_reference_set/reads')

reference_set_dict = {}

# sample uniformly positive examples from each clade
for i in range(int(n/2)):
    # sample a clade uniformly
    clade = clade_lineage_info['GISAID_clade'].sample(1).values[0]
    # sample a lineage uniformly from the clade
    lineage = clade_lineage_info[clade_lineage_info['GISAID_clade'] == clade]['pango_lineage'].sample(1).values[0]

    # sample a positive example from the lineage
    read1_seq, read2_seq, read1_id, read_2_id, id, _ = sample_examples.sample_positive(lineage)

    while reference_set_dict.keys().__contains__("{}_{}".format(read1_id, read_2_id)) or reference_set_dict.keys().__contains__("{}_{}".format(read_2_id, read1_id)):
        read1_seq, read2_seq, read1_id, read_2_id, id, _ = sample_examples.sample_positive(lineage)

    # write read1 to a file
    with open("data/data/hcov_global_2023-11-16_09-28/siamese_reference_set/reads/{}.fasta".format(read1_id), "w") as outfile:
        outfile.write(">{}\n".format(read1_id))
        outfile.write("{}\n".format(read1_seq))
    
    # write read2 to a file
    with open("data/data/hcov_global_2023-11-16_09-28/siamese_reference_set/reads/{}.fasta".format(read_2_id), "w") as outfile:
        outfile.write(">{}\n".format(read_2_id))
        outfile.write("{}\n".format(read2_seq))

    reference_set_dict["{}_{}".format(read1_id, read_2_id)] ={
    'clade':clade,
    'lineage': lineage,
    'label': "positive",
    'read1': read1_id,
    'read2':read_2_id}

    ref_set_stats.loc[(ref_set_stats['GISAID_clade'] == clade) & (ref_set_stats['pango_lineage'] == lineage), 'count'] += 2
# sample uniformly negative examples from the same  clade and lineage
for i in range(int(n/4)): 

    # sample a clade uniformly
    clade = clade_lineage_info['GISAID_clade'].sample(1).values[0]
    # sample a lineage uniformly from the clade
    lineage = clade_lineage_info[clade_lineage_info['GISAID_clade'] == clade]['pango_lineage'].sample(1).values[0]

    # check how many samples are in the lineage
    lineage_dir = os.path.join(data_dir, lineage)
    lineage_count = len([name for name in os.listdir(lineage_dir) if name.endswith(".fasta")])
    while lineage_count < 2:
        clade = clade_lineage_info['GISAID_clade'].sample(1).values[0]
        lineage = clade_lineage_info[clade_lineage_info['GISAID_clade'] == clade]['pango_lineage'].sample(1).values[0]
        # count number of files that end in .fasta in the lineage directory
        lineage_dir = os.path.join(data_dir, lineage)
        lineage_count = len([name for name in os.listdir(lineage_dir) if name.endswith(".fasta")])        
    
    read1_seq, read2_seq, read1_id, read_2_id, id1, id2 =  sample_examples.sample_negative(True, lineage)


    while reference_set_dict.keys().__contains__("{}_{}".format(read1_id, read_2_id)) or reference_set_dict.keys().__contains__("{}_{}".format(read_2_id, read1_id)):
        read1_seq, read2_seq, read1_id, read_2_id, id1, id2 = sample_examples.sample_negative(True, lineage)
    
     # write read1 to a file
    with open("data/data/hcov_global_2023-11-16_09-28/siamese_reference_set/reads/{}.fasta".format(read1_id), "w") as outfile:
        outfile.write(">{}\n".format(read1_id))
        outfile.write("{}\n".format(read1_seq))
    
    # write read2 to a file
    with open("data/data/hcov_global_2023-11-16_09-28/siamese_reference_set/reads/{}.fasta".format(read_2_id), "w") as outfile:
        outfile.write(">{}\n".format(read_2_id))
        outfile.write("{}\n".format(read2_seq))

    reference_set_dict["{}_{}".format(read1_id, read_2_id)] ={
    'clade':clade,
    'lineage': lineage,
    'label': "negative",
    'read1': read1_id,
    'read2':read_2_id}

    ref_set_stats.loc[(ref_set_stats['GISAID_clade'] == clade) & (ref_set_stats['pango_lineage'] == lineage), 'count'] += 2

# sample uniformly negative examples from different clade/lineage
for i in range(int(n/4)):
    # sample a clade uniformly
    clade_1 = clade_lineage_info['GISAID_clade'].sample(1).values[0]
    # sample a lineage uniformly from the clade
    lineage_1 = clade_lineage_info[clade_lineage_info['GISAID_clade'] == clade_1]['pango_lineage'].sample(1).values[0]

    clade_2 = clade_lineage_info['GISAID_clade'].sample(1).values[0]
    # sample a lineage uniformly from the clade
    lineage_2 = clade_lineage_info[clade_lineage_info['GISAID_clade'] == clade_2]['pango_lineage'].sample(1).values[0]

    if lineage_1 == lineage_2:

        clade_2 = clade_lineage_info['GISAID_clade'].sample(1).values[0]
         # sample a lineage uniformly from the clade
        lineage_2 = clade_lineage_info[clade_lineage_info['GISAID_clade'] == clade_2]['pango_lineage'].sample(1).values[0]


    read1_seq, read2_seq, read1_id, read_2_id, id1, id2 =  sample_examples.sample_negative(False, lineage_1, lineage_2)

    while reference_set_dict.keys().__contains__("{}_{}".format(read1_id, read_2_id)) or reference_set_dict.keys().__contains__("{}_{}".format(read_2_id, read1_id)):
        read1_seq, read2_seq, read1_id, read_2_id, id1, id2 =  sample_examples.sample_negative(False, lineage_1, lineage_2)
    
     # write read1 to a file
    with open("data/data/hcov_global_2023-11-16_09-28/siamese_reference_set/reads/{}.fasta".format(read1_id), "w") as outfile:
        outfile.write(">{}\n".format(read1_id))
        outfile.write("{}\n".format(read1_seq))
    
    # write read2 to a file
    with open("data/data/hcov_global_2023-11-16_09-28/siamese_reference_set/reads/{}.fasta".format(read_2_id), "w") as outfile:
        outfile.write(">{}\n".format(read_2_id))
        outfile.write("{}\n".format(read2_seq))

    reference_set_dict["{}_{}".format(read1_id, read_2_id)] ={
    'clade1':clade_1,
    'lineage1': lineage_1,
    'clade2':clade_2,
    'lineage2': lineage_2,
    'label': "negative",
    'read1': read1_id,
    'read2':read_2_id}

    ref_set_stats.loc[(ref_set_stats['GISAID_clade'] == clade_1) & (ref_set_stats['pango_lineage'] == lineage_1), 'count'] += 1
    ref_set_stats.loc[(ref_set_stats['GISAID_clade'] == clade_2) & (ref_set_stats['pango_lineage'] == lineage_2), 'count'] += 1

# write reference set dictionary to a json file
with open("data/data/hcov_global_2023-11-16_09-28/siamese_reference_set/reference_set.json", "w") as outfile: 
    json.dump(reference_set_dict, outfile)

# write reference set stats to a csv file
ref_set_stats.to_csv('data/data/hcov_global_2023-11-16_09-28/siamese_reference_set/reference_set_stats.csv', index=False)




