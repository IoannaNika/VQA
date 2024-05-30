import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def init_n_strains(genomic_regions):
    n_strains = {}
    for region in genomic_regions:
        n_strains[region] = 0

    return n_strains

def process_df(df):
    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                       (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                        (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                         (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                        (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                         (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                          (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]


    genomic_regions = [ str(gr[0])  + "_" + str(gr[1]) for gr in genomic_regions]
    
    n_strains = init_n_strains(genomic_regions)


    for index, row in df.iterrows():
        region = row["Genomic_region"]
        n_strains[region] += 1
    
    return n_strains

def get_expected_n_of_strains(genomic_regions):
    expected_n_strains = {}
    for region in genomic_regions:
        expected_n_strains[region] = 2
        if region in ["1128_2244", "3166_4240","15634_16698", "16647_17732", "19604_20676", "29768_29790"]:
            expected_n_strains[region] = 1
    return expected_n_strains

def plot_swarmplot(n_strains, sample_name, clr, offset):
    # n_strains is a dictionary with the genomic region as key and the number of strains as value
    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                       (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                        (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                         (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                        (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                         (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                          (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]


    genomic_regions = [ str(gr[0])  + "_" + str(gr[1]) for gr in genomic_regions]
    
    true_n_strains = n_strains.values()

    x = [i + offset for i in range(len(genomic_regions))]
    swarm = plt.scatter(x, list(true_n_strains), color=clr, label=sample_name, s=20)

    return swarm

def main(): 
    colors = sns.color_palette("Set2", 5)
    input_03 = "Experiments/lumc_subsample/03_50/consensus_lumc_comparison_post_processed.tsv"

    df_03 = pd.read_csv(input_03, sep="\t", header=0)
    sample_name = "03_50"
    n_strains_03 = process_df(df_03)
    swarm = plot_swarmplot(n_strains_03, sample_name, colors[0], -0.2)

    input_04 = "Experiments/lumc_subsample/04_75/consensus_lumc_comparison_post_processed.tsv"

    df_04 = pd.read_csv(input_04, sep="\t", header=0)
    sample_name = "04_75"
    n_strains_04 = process_df(df_04)
    swarm = plot_swarmplot(n_strains_04, sample_name, colors[1], -0.05)

    input_05 = "Experiments/lumc_subsample/05_90/consensus_lumc_comparison_post_processed.tsv"

    df_05 = pd.read_csv(input_05, sep="\t", header=0)
    sample_name = "05_90"
    n_strains_05 = process_df(df_05)
    swarm = plot_swarmplot(n_strains_05, sample_name, colors[2], 0.05)

    input_06 = "Experiments/lumc_subsample/06_95/consensus_lumc_comparison_post_processed.tsv"

    df_06 = pd.read_csv(input_06, sep="\t", header=0)
    sample_name = "06_95"
    n_strains_06 = process_df(df_06)
    swarm = plot_swarmplot(n_strains_06, sample_name, colors[4], 0.05)


    input_07 = "Experiments/lumc_subsample/07_98/consensus_lumc_comparison_post_processed.tsv"
    
    df_07 = pd.read_csv(input_07, sep="\t", header=0)
    sample_name = "07_98"
    n_strains_07 = process_df(df_07)
    swarm = plot_swarmplot(n_strains_07, sample_name, colors[3], 0.3)


    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                       (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                        (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                         (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                        (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                         (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                          (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]


    genomic_regions = [ str(gr[0])  + "_" + str(gr[1]) for gr in genomic_regions]
    


    expected_n_strains = get_expected_n_of_strains(genomic_regions)

    
    # plot expected number of strains
    sns.lineplot(list(expected_n_strains.values()), linestyle ='--', label="Expected number of strains", color="black", linewidth=0.5,  marker=".", markersize=7)
    # plt.setp(swarm.lines, zorder=1000)

    # xticks
    xticks = [region.replace("_", "-")for region in genomic_regions]
    # rotate xticks
    plt.xticks(ticks=np.arange(len(xticks)), labels=xticks, rotation=90)
    plt.ylim(-0.5, 7)
    # y ticks   
    plt.yticks(np.arange(0, 7, 1))
    plt.xlabel("Genomic region")
    plt.ylabel("Number of haplotypes")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.10), ncol=6, fontsize=7)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.grid(axis='x', linestyle='--', alpha=0.5)
    # y-axis should start from 0 but a bit more up to have a better view
    plt.subplots_adjust(bottom=0.7)
    plt.tight_layout()
    plt.savefig(f"Experiments/lumc_subsample/n_strains.pdf", bbox_inches='tight')







if __name__ == "__main__":
    main()