import argparse
import sys
import os
from Bio import SeqIO
import editdistance
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--template_dir', dest = 'template_dir', required=True, type=str, help="Folder with templates")
    parser.add_argument('--data_dir', dest = 'data_dir', required=True, type=str, help="Folder with data, will be used as output folder too")
    parser.add_argument('--n_true_haplotypes', dest = 'n_true_haplotypes', required=True, type=int, help="Number of true haplotypes")
    args = parser.parse_args()

    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                        (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                        (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                        (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                        (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                        (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                        (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]

    consensus = os.path.join(args.data_dir, "consensus.fasta")

    # find templates
    templates = dict()

    for gr in genomic_regions:
        gr = str(gr[0]) + "_" + str(gr[1])
        templates[gr] = {}

    os.system("touch {}".format(os.path.join(args.data_dir, "temp.fasta")))

    for file in os.listdir(args.template_dir):
        if file.endswith(".template"):
            identifier = file[:-9]
            for gr in genomic_regions:
                gr = str(gr[0]) + "_" + str(gr[1])

                os.system("grep -A 1 +_{}:{}:0 {} > {}".format(identifier, gr, os.path.join(args.template_dir, file), os.path.join(args.data_dir, "temp.fasta")))
                
                try:
                    template = next(SeqIO.parse(os.path.join(args.data_dir, "temp.fasta"), "fasta"))
                except: 
                    continue
               
                templates[gr][identifier] = template.seq
    

    # load consensus
    consensus_sequences = dict()

    for gr in genomic_regions:
        gr =  str(gr[0]) + "_" + str(gr[1])
        consensus_sequences[gr] = []

    for record in SeqIO.parse(consensus, "fasta"):
        community = record.id.split("_")[0]
        # >0_54_1183_0.3541666666666667
        genomic_region =str(record.id.split("_")[1]) + "_" + str(record.id.split("_")[2])
        predicted_abundance = float(record.id.split("_")[3])
        consensus_sequences[genomic_region].append((predicted_abundance, record))

    # output tsv file with results
    # header: genomic_region, template_id, consensus_id ,predicted_abundance, similarity
    output_file = os.path.join(args.data_dir, "meta_results_2.tsv")
    with open(output_file, "w") as f:
        f.write("genomic_region\ttemplate_id\tconsensus_id\tpredicted_abundance\tabs_relative_ab_error\tedit_distance\tnr_consensus\n")

    # for each reconstructed haplotype find the  templates it represents best
    # temporarily predicted abundance per template is the predicted abundance per consensus
    # temporarily adjusted abundance is the predicted abundance of the consensus / number of templates it represents best
    # temporarily the rel_ab_error is calculated with the adjusted abundance (mentioned above)
    for gr in genomic_regions: 
        gr = str(gr[0]) + "_" + str(gr[1])
        
        for consensus_sequence in consensus_sequences[gr]:
            predicted_abundance, consensus_sequence = consensus_sequence

            min_distance = float("inf")
            min_template = None
            
            result = {}

            result[min_distance] = set()
            result[min_distance].add(min_template)  

            for template in templates[gr]:
                template_sequence = templates[gr][template]
                # compare sequences
                similarity = editdistance.eval(consensus_sequence.seq, template_sequence)
                if similarity < min_distance:
                    min_distance = similarity
                    min_template = template

                    if min_distance not in result:
                        result[min_distance] = set()
                
                if similarity == min_distance:
                    result[min_distance].add(template)

            # find the key with the minimum edit distance
            min_distance = min(result.keys())

            for template in result[min_distance]:  
                with open(output_file, "a") as f:
                    true_abundance = 1/args.n_true_haplotypes
                    adjusted_ab = "{:.3f}".format((float(predicted_abundance)/len(result[min_distance])))
                    abs_relative_ab_error = "{:.3f}".format(abs(float(adjusted_ab) - true_abundance)/true_abundance)
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gr, template, consensus_sequence.id, predicted_abundance, abs_relative_ab_error, min_distance, 1))
                    f.close()

    # open the output file and per genomic region if a template is given more than once, keep the one with the smallest edit distance
    # load meta results
    meta_results = pd.read_csv(output_file, sep="\t", header=0)
    for gr in genomic_regions:
        gr = str(gr[0]) + "_" + str(gr[1])
        gr_results = meta_results[meta_results["genomic_region"] == gr]
        for template in gr_results["template_id"].unique():
            template_results = gr_results[gr_results["template_id"] == template]
            if len(template_results) > 1:
                min_distance = min(template_results["edit_distance"])
                for index, row in meta_results.iterrows():
                    if row["edit_distance"] != min_distance and row["template_id"] == template and row["genomic_region"] == gr:
                        meta_results.drop(index, inplace=True)
        
    
    meta_results.to_csv(output_file, sep="\t", index=False)

    # open it again and calculate and replace the predicted abundance with the adjusted abundance and the absolute relative abundance error
    # load meta results
    meta_results = pd.read_csv(output_file, sep="\t", header=0)

    for gr in genomic_regions:
        gr = str(gr[0]) + "_" + str(gr[1])
        gr_results = meta_results[meta_results["genomic_region"] == gr]
        for consensus_id in gr_results["consensus_id"].unique():
            consensus_results = gr_results[gr_results["consensus_id"] == consensus_id]
            n_templates = len(consensus_results)
            true_abundance = 1/args.n_true_haplotypes
            for index, row in meta_results.iterrows():
                if row["consensus_id"] == consensus_id and row["genomic_region"] == gr:
                    adjusted_ab = "{:.3f}".format((float(row["predicted_abundance"])/n_templates))
                    meta_results.at[index, "predicted_abundance"] = adjusted_ab
    
    meta_results.to_csv(output_file, sep="\t", index=False)

    # if a template is represented more than once with the same edit distance add the predicted abundances from all consensus it is represented
    # and keep only one row, the number of consensus and the consensus ids
    for gr in genomic_regions:
        gr = str(gr[0]) + "_" + str(gr[1])
        gr_results = meta_results[meta_results["genomic_region"] == gr]
        
        for template in gr_results["template_id"].unique():
            template_results = gr_results[gr_results["template_id"] == template]
            if len(template_results) > 1:
                predicted_abundance = sum(template_results["predicted_abundance"].values)
                for index, row in meta_results.iterrows():
                    if row["template_id"] == template and row["genomic_region"] == gr:
                        meta_results.at[index, "predicted_abundance"] = predicted_abundance
                        meta_results.at[index, "nr_consensus"] = len(template_results)
                        # concatenate consensus ids
                        conc = ""
                        for consensus_id in template_results["consensus_id"].values:
                            conc += consensus_id + "-"
                        meta_results.at[index, "consensus_id" ] = conc[:-1]
                meta_results.drop(template_results.index[1:], inplace=True)
    

    meta_results.to_csv(output_file, sep="\t", index=False)

    # normalize relative abundance to sum to 1 per genomic region and then calculate relative abundance error
    # load meta results
    meta_results = pd.read_csv(output_file, sep="\t", header=0)
    for gr in genomic_regions:
        gr = str(gr[0]) + "_" + str(gr[1])
        gr_results = meta_results[meta_results["genomic_region"] == gr]
       
        sum_abundance = sum(gr_results["predicted_abundance"].values)
        for index, row in meta_results.iterrows():
            if row["genomic_region"] == gr:
                true_abundance = 1/args.n_true_haplotypes
                adjusted_ab = "{:.3f}".format((float(row["predicted_abundance"])/sum_abundance))
                meta_results.at[index, "predicted_abundance"] = adjusted_ab
                meta_results.at[index, "abs_relative_ab_error"] = "{:.3f}".format(abs(float(adjusted_ab) - true_abundance)/true_abundance)

    meta_results.to_csv(output_file, sep="\t", index=False)

    # remove temp.fasta file
    os.system("rm {0}".format(os.path.join(args.data_dir, "temp.fasta")))

 
if __name__ == "__main__":
    sys.exit(main())