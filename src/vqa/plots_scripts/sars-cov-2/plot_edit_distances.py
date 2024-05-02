import argparse
import sys
import os
from Bio import SeqIO
import editdistance
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', dest = 'outdir', required=True, type=str, help="Output folder")
    parser.add_argument('--data_dir', dest = 'data_dir', required=True, type=str, help="Folder with templates")
    args = parser.parse_args()

    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                        (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                        (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                        (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                        (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                        (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                        (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]
    
    haplotypes = [5, 10, 15, 20]
    templates = dict()

    for haplotype in haplotypes:
       
        templates[haplotype] = dict()
       
        template_dir = os.path.join(args.data_dir, str(haplotype), "NCBI_processed")

        dirs = os.listdir(template_dir)
        # filter out the files that do not end with .template
        dirs = [file for file in dirs if file.endswith(".template")]
       
        for file in dirs:
            identifier = file[:-9]
            for gr in genomic_regions:
                gr = str(gr[0]) + "_" + str(gr[1])

                if gr not in templates[haplotype].keys(): 
                    templates[haplotype][gr] = dict()

                os.system("grep -A 1 +_{}:{}:0 {} > {}".format(identifier, gr, os.path.join(template_dir, file), os.path.join(args.outdir, "temp.fasta")))
                
                try:
                    template = next(SeqIO.parse(os.path.join(args.outdir, "temp.fasta"), "fasta"))
                except: 
                    continue
                templates[haplotype][gr][identifier] = template.seq  

    
    # calculate edit distances
    edit_distances = dict()

    for haplotype in haplotypes:
        edit_distances[haplotype] = dict()
        for gr in genomic_regions:
            gr = str(gr[0]) + "_" + str(gr[1])
            edit_distances[haplotype][gr] = []

            for i, identifier in enumerate(list(set(templates[haplotype][gr].keys()))):

                for j, identifier2 in enumerate(list(set(templates[haplotype][gr].keys()))):
                    if j > i and identifier != identifier2:
                        edit_distances[haplotype][gr].append(editdistance.eval(templates[haplotype][gr][identifier], templates[haplotype][gr][identifier2]))

    # plot results
    fig, ax = plt.subplots(1, len(haplotypes), figsize=(50, 10))
    grs = []
    for gr in genomic_regions:
        gr = str(gr[0]) + "_" + str(gr[1])
        grs.append(gr)

    # for each genomic region have a boxplot
    for i, haplotype in enumerate(haplotypes):

        for gr in genomic_regions:
            gr = str(gr[0]) + "_" + str(gr[1])
        
            ax[i].boxplot(edit_distances[haplotype][gr], positions=[grs.index(gr)], showfliers=False)


        ax[i].set_title("Number of true haplotypes {}".format(haplotype))
        ax[i].set_xlabel("Genomic region")
        ax[i].set_ylabel("Edit distance")
        ax[i].tick_params(labelrotation=90)
        ax[i].set_xticklabels(grs)


    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "edit_distances.png"))

if __name__ == "__main__":
    sys.exit(main())
        
        

