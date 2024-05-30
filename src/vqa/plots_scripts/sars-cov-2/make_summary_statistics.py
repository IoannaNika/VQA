import argparse
import json
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
import editdistance


def get_templates(data_dir, haplotypes, outdir, genomic_regions):
    templates = dict()

    for haplotype in haplotypes:
       
        templates[haplotype] = dict()
       
        template_dir = os.path.join(data_dir, str(haplotype), "NCBI_processed")

        dirs = os.listdir(template_dir)
        # filter out the files that do not end with .template
        dirs = [file for file in dirs if file.endswith(".template")]
        for file in dirs:
            identifier = file[:-9]
            for gr in genomic_regions:
                gr = str(gr[0]) + "_" + str(gr[1])

                if gr not in templates[haplotype].keys(): 
                    templates[haplotype][gr] = dict()

                os.system("grep -A 1 +_{}:{}:0 {} > {}".format(identifier, gr, os.path.join(template_dir, file), os.path.join(outdir, "temp.fasta")))
                
                try: 
                    template = next(SeqIO.parse(os.path.join(outdir, "temp.fasta"), "fasta"))
                except: 
                    continue
                templates[haplotype][gr][identifier] = template.seq  

    return templates 


# def n_true_haplotypes_func(templates): 
#     n_true_haplotypes = dict()
#     n_true_haplotypes_n = dict()
    
#     # get the number of true haplotypes, per genomic region if sequences are the same they are considered the same haplotype

#     for h in templates.keys():
#         n_true_haplotypes[h]  = dict()
#         n_true_haplotypes_n[h] = dict()

#         for gr in templates[h].keys():
#             n_true_haplotypes[h][gr] = set(templates[h][gr].values())
#             n_true_haplotypes_n[h][gr] = len(n_true_haplotypes[h][gr])
    
#     avg_over_genomic_regions  = dict()
#     for h in n_true_haplotypes_n.keys():
#         print(n_true_haplotypes_n[h].values())
#         avg_over_genomic_regions[h] = np.mean(list(n_true_haplotypes_n[h].values()))

#     return avg_over_genomic_regions


def get_accuracy(templates, haplotypes, outdir):

    overall_accuracy = dict()

    for h in haplotypes:
        overall_accuracy[h] = dict()
        # initialize the dictionary
        
        overall_accuracy[h] = 0
        predictions_file = os.path.join(outdir, str(h), "predictions.tsv")

        predictions = pd.read_csv(predictions_file, sep="\t", header=0)

        # iterate over the rows of the predictions file
        for index, row in predictions.iterrows():
            gr = row["Genomic_region"]

            seq_id_1 = row["Sequence_1_id"].split("_")[0]
            seq_id_2 = row["Sequence_2_id"].split("_")[0]

            template_1 = templates[h][gr][seq_id_1]
            template_2 = templates[h][gr][seq_id_2]

            if row["Predicted_label"] == row["True_label"] and row["True_label"] == 1:
                overall_accuracy[h] += 1
            
            if row["Predicted_label"] == row["True_label"] and row["True_label"] == 0 and editdistance.eval(template_1, template_2) != 0:
                overall_accuracy[h] += 1

            if row["Predicted_label"] != row["True_label"] and row["True_label"] == 0 and editdistance.eval(template_1, template_2) == 0:
                overall_accuracy[h] += 1
        

        overall_accuracy[h] = overall_accuracy[h] / len(predictions)
        print(h,overall_accuracy[h])

    return overall_accuracy



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', dest = 'directory', required=True, type=str, help="")
    parser.add_argument('--outdir', dest = 'outdir', required=True, type=str, help="Output folder")
    parser.add_argument('--data_dir', dest = 'data_dir', required=True, type=str, help="Folder with templates")
    args = parser.parse_args()

    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                        (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                        (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                        (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                        (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                        (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                        (26766, 27872), (27808, 28985), (28699, 29768)]

    haplotypes = [1, 5, 10, 15, 20]
    
    templates = get_templates(args.data_dir, haplotypes, args.outdir, genomic_regions)


    
    edit_distances = dict()
    rel_ab_error = dict()
    n_reconstructed_haplotypes = dict()
    accuracy = get_accuracy(templates, haplotypes, args.outdir)

    for h in haplotypes:
        meta_results = os.path.join(args.directory, str(h), "meta_results_2.tsv")
        # load meta results
        meta_results = pd.read_csv(meta_results, sep="\t", header=0)

        # edit_distances[h] = np.mean(meta_results["edit_distance"])
        # rel_ab_error[h] = np.mean(meta_results["abs_relative_ab_error"])

        haplotypes_across_regions = dict()
        rel_ab_error_across_regions = dict()
        edit_distances_across_regions = dict()
        for gr in genomic_regions:
            gr = str(gr[0]) + "_" + str(gr[1])
            haplotypes_across_regions[gr] = len(meta_results[meta_results["genomic_region"] == gr])
            rel_ab_error_across_regions[gr] = np.mean(list(meta_results[meta_results["genomic_region"] == gr]["abs_relative_ab_error"]) + [1] * min(h - haplotypes_across_regions[gr], 0))
            edit_distances_across_regions[gr] = np.mean(list(meta_results[meta_results["genomic_region"] == gr]["edit_distance"]))
        
        n_reconstructed_haplotypes[h] = np.mean(list(haplotypes_across_regions.values()))
        rel_ab_error[h] = np.mean(list(rel_ab_error_across_regions.values()))
        edit_distances[h] = np.mean(list(edit_distances_across_regions.values()))


    # write to csv 
    with open(os.path.join(args.outdir, "summary_statistics.csv"), 'w') as f:
        f.write("True_haplotypes,Edit_distance,Relative_abundance_error,Number_of_reconstructed_haplotypes\n")
        for h in haplotypes:
            f.write("{},{},{},{},{}\n".format(h, edit_distances[h], rel_ab_error[h], accuracy[h], n_reconstructed_haplotypes[h]))
    

        
    
  

if __name__ == "__main__":
    sys.exit(main())