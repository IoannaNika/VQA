import argparse
import sys
import pandas as pd
import editdistance
import matplotlib.pyplot as plt
import os
from Bio import SeqIO   
import pickle 
import numpy as np

def is_the_coverage_sufficient(reads, gr, low_limit=100):
    reads_region = reads[reads["Genomic_regions"] == gr]
    reads_coverage = len(reads_region)
    if reads_coverage >= low_limit:
        return True
    else:
        return False

def true_abundances():

    true_abundances = { 
        "01_100": {"BA.1":0, "Wuhan": 1},
        "02_100": {"BA.1":0, "Wuhan": 1},
        "03_50": {"BA.1":0.5, "Wuhan": 0.5},
        "04_75": {"BA.1":0.25, "Wuhan": 0.75},
        "05_90": {"BA.1":0.1, "Wuhan": 0.9},
        "06_95": {"BA.1":0.05, "Wuhan": 0.95},
        "07_98": {"BA.1":0.02, "Wuhan": 0.98},
        "08_0": {"BA.1":1.0, "Wuhan": 0.0},
        "09_0": {"BA.1":1.0, "Wuhan": 0.0}
    }

    return true_abundances


    
def regions_with_sufficie_coverage_in_sample(reads, low_limit=100):
    genomic_regions = reads["Genomic_regions"].unique()
    regions_with_sufficient_coverage = []
    for gr in genomic_regions:
        if is_the_coverage_sufficient(reads, gr, low_limit):
            regions_with_sufficient_coverage.append(gr)
    return regions_with_sufficient_coverage

def get_expected_n_of_strains(genomic_regions, samples):
    expected_n_strains = {}
    for sample in samples:
        expected_n_strains[sample] = {}
        for region in genomic_regions:
            if sample in ["03_50", "04_75", "05_90", "06_95", "07_98"]:
                expected_n_strains[sample][region] = 2
                if region in ["1128_2244", "3166_4240","15634_16698", "16647_17732", "19604_20676", "29768_29790"]:
                    expected_n_strains[region] = 1
            elif sample in ["01_100", "02_100", "08_0", "09_0"]:
                expected_n_strains[sample][region] = 1
    return expected_n_strains

def get_grs_with_sufficient_coverage(samples, low_limit=100):
    sufficient_coverage_per_sample = {}

    for sample in samples:
        reads_file = "Experiments/lumc_subsample/" + sample + '/communities.tsv'
        reads = pd.read_csv(reads_file, sep="\t", header=0)
        sufficient_coverage_per_sample[sample] = regions_with_sufficie_coverage_in_sample(reads, low_limit=100)
    
    return sufficient_coverage_per_sample


def calculate_rel_ab_error(true_abundance_per_sample, rel_abundances, sufficient_coverage_per_sample,  directory="Experiments/lumc_subsample/"):
    rel_ab_errors = {}
    for sample in rel_abundances.keys():
        sample_dir = os.path.join(directory, sample)
        input_file = os.path.join(sample_dir, "consensus_lumc_comparison_post_processed.tsv")
        df = pd.read_csv(input_file, sep="\t", header=0)
        df_same = list(df[df["Edit_distance_between_consensus"] == 0]['Genomic_region'])
        rel_ab_errors[sample] = {}
        for strain in true_abundance_per_sample[sample].keys():
            rel_ab_errors[sample][strain] = {}
            for gr in rel_abundances[sample][strain].keys():
                if gr not in sufficient_coverage_per_sample[sample]:
                    continue

                rel_ab_errors[sample][strain][gr] = abs(true_abundance_per_sample[sample][strain] - rel_abundances[sample][strain][gr]) / ((true_abundance_per_sample[sample][strain]) + 1e-6)
                
                if gr in df_same:
                    # if the abundance of one of the strains is 0 and the other is 1, the error is 0
                    min_distance = np.argmin([abs(1 - rel_abundances[sample][strain][gr]), abs(0 - rel_abundances[sample][strain][gr])]) 
                    if min_distance == 0:
                        rel_ab_errors[sample][strain][gr] = abs(1 - rel_abundances[sample][strain][gr]) /((true_abundance_per_sample[sample][strain]) + 1e-6)
                    else:
                        rel_ab_errors[sample][strain][gr] = abs(0 - rel_abundances[sample][strain][gr]) / ((true_abundance_per_sample[sample][strain]) + 1e-6)
                   
    return rel_ab_errors


def main(): 
    parser = argparse.ArgumentParser(description="This script will calculate the consensus distance from reads")
    args = parser.parse_args()

    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                            (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                            (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                            (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                            (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                            (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                            (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]
    
    gr_merged = [str(gr[0]) + "_" + str(gr[1]) for gr in genomic_regions]


    # samples = ["01_100", "02_100", "08_0", "09_0", "03_50", "04_75", "05_90", "06_95", "07_98"]
    # samples = ["03_50", "04_75", "05_90", "06_95", "07_98"]
    samples = ["01_100", "02_100",  "08_0", "09_0"]
    expected_number_of_strains = get_expected_n_of_strains(gr_merged, samples)
    true_abundance_per_sample = true_abundances()
    sufficient_coverage_per_sample = get_grs_with_sufficient_coverage(samples)

    with open('Experiments/lumc_subsample/sufficient_coverage_per_sample.pickle', 'wb') as file:
        pickle.dump(sufficient_coverage_per_sample, file)
    
    min_edit_distances_go_viral = {}
    number_of_haplotypes_go_viral = {}
    relative_abundance_error_go_viral = {}

    for sample in samples:
        consensus_file = "Experiments/lumc_subsample/" + sample + '/consensus_lumc_comparison_post_processed.tsv'
        consensus = pd.read_csv(consensus_file, sep="\t", header=0)

        min_edit_distances_go_viral[sample] = {}
        number_of_haplotypes_go_viral[sample] = {}
        relative_abundance_error_go_viral[sample] = {}

        for gr in consensus["Genomic_region"].unique():
            if gr not in sufficient_coverage_per_sample[sample]:
                continue

            if gr not in min_edit_distances_go_viral[sample]:
                min_edit_distances_go_viral[sample][gr] = []
                number_of_haplotypes_go_viral[sample][gr] = len(consensus[consensus["Genomic_region"] == gr])
                relative_abundance_error_go_viral[sample][gr] = []

            consensus_region = consensus[consensus["Genomic_region"] == gr]
            for i, haplotype_row in consensus_region.iterrows():
                min_edit_distances_go_viral[sample][gr].append(haplotype_row["Min_edit_distance"])
                

    # get average of min edit distances across genomic regions and samples 
    avg_min_ed = 0

    for sample in min_edit_distances_go_viral.keys():
        avg_min_ed += np.mean([np.mean(min_edit_distances_go_viral[sample][gr]) for gr in min_edit_distances_go_viral[sample].keys()])

    avg_min_ed = avg_min_ed / len(samples)

    print("Average min edit distance: ", avg_min_ed)

    # get average number of expected and predicted number of samples

    avg_n_haplotypes_predicted = 0
    avg_n_haplotypes_expected = 0

    for sample in number_of_haplotypes_go_viral.keys():
        avg_n_haplotypes_predicted += np.mean([number_of_haplotypes_go_viral[sample][gr] for gr in number_of_haplotypes_go_viral[sample].keys()])
        avg_n_haplotypes_expected += np.mean([expected_number_of_strains[sample][gr] for gr in expected_number_of_strains[sample].keys()])

    avg_n_haplotypes_predicted = avg_n_haplotypes_predicted / len(samples)
    avg_n_haplotypes_expected = avg_n_haplotypes_expected / len(samples)

    print("Average number of haplotypes predicted: ", avg_n_haplotypes_predicted)
    print("Average number of haplotypes expected: ", avg_n_haplotypes_expected)

    # # load relative abundances from Experiments/lumc_subsample/rel_abundances_mixture.pkl
    # rel_abundances = pickle.load(open("Experiments/lumc_subsample/rel_abundances_mixture.pkl", "rb"))

    # load relative abundances from Experiments/lumc_subsample/rel_abundances_mixture.pkl
    rel_abundances = pickle.load(open("Experiments/lumc_subsample/rel_abundances_pures.pkl", "rb"))

    # include the relative abundances of the pure strains
    # rel_abundances_pures = pickle.load(open("Experiments/lumc_subsample/rel_abundances_pures.pkl", "rb"))
    # rel_abundances.update(rel_abundances_pures)
    # calculate relative abundance error
    rel_ab_errors = calculate_rel_ab_error(true_abundance_per_sample, rel_abundances, sufficient_coverage_per_sample)

    avg_rel_ab_error = 0
    for sample in rel_ab_errors.keys():
        avg_rel_ab_error_strain_1 = np.mean([rel_ab_errors[sample]["BA.1"][gr] for gr in rel_ab_errors[sample]["BA.1"].keys()])
        avg_rel_ab_error_strain_2 = np.mean([rel_ab_errors[sample]["Wuhan"][gr] for gr in rel_ab_errors[sample]["Wuhan"].keys()])
        avg_rel_ab_error += (avg_rel_ab_error_strain_1 + avg_rel_ab_error_strain_2) / 2
           

    avg_rel_ab_error = avg_rel_ab_error / len(samples)   
    print("Average relative abundance error: ", avg_rel_ab_error)

    # calculate the average difference between expected and predicted number of strains
    avg_diff_n_strains = 0
    for sample in samples:
        avg_diff_n_strains += np.mean([abs(expected_number_of_strains[sample][gr] -number_of_haplotypes_go_viral[sample][gr]) for gr in number_of_haplotypes_go_viral[sample].keys()])

        # if sample in ["03_50", "02_100"]:
        #     avg_diff_n_strains += np.mean([abs(expected_number_of_strains[sample][gr] -0) for gr in number_of_haplotypes_go_viral[sample].keys()])
        # else: 
        #     avg_diff_n_strains += np.mean([abs(expected_number_of_strains[sample][gr] -1) for gr in number_of_haplotypes_go_viral[sample].keys()])


    avg_diff_n_strains = avg_diff_n_strains / len(samples)

    print("Average difference between expected and predicted number of strains: ", avg_diff_n_strains)
            
    
  
            
if __name__ == "__main__":
    sys.exit(main())
