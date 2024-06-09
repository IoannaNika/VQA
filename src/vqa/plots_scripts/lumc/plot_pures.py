import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def is_the_coverage_sufficient(reads, gr, low_limit=100):
    reads_region = reads[reads["Genomic_regions"] == gr]
    reads_coverage = len(reads_region)
    if reads_coverage >= low_limit:
        return True
    else:
        return False
    
def regions_with_sufficie_coverage_in_sample(reads, low_limit=100):
    genomic_regions = reads["Genomic_regions"].unique()
    regions_with_sufficient_coverage = []
    for gr in genomic_regions:
        if is_the_coverage_sufficient(reads, gr, low_limit):
            regions_with_sufficient_coverage.append(gr)
    return regions_with_sufficient_coverage

def get_grs_with_sufficient_coverage(samples, low_limit=100):
    sufficient_coverage_per_sample = {}

    for sample in samples:
        reads_file = "Experiments/lumc_subsample/" + sample + '/communities.tsv'
        reads = pd.read_csv(reads_file, sep="\t", header=0)
        sufficient_coverage_per_sample[sample] = regions_with_sufficie_coverage_in_sample(reads, low_limit=100)
    
    return sufficient_coverage_per_sample


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
        expected_n_strains[region] = 1
    return expected_n_strains

def plot_swarmplot(n_strains, sample_name, clr, offset, high_coverage_regions):
    # n_strains is a dictionary with the genomic region as key and the number of strains as value
    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                       (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                        (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                         (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                        (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                         (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                          (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]


    genomic_regions = [ str(gr[0])  + "_" + str(gr[1]) for gr in genomic_regions]

    # remove the regions with low coverage
    n_strains = [n_strains[gr] if gr in high_coverage_regions else np.nan for gr in genomic_regions]
    # plot swarmplot the x-axis is the genomic region and the y-axis is the number of strains
    x = [i + offset for i in range(len(genomic_regions))]
    swarm = plt.scatter(x, list(n_strains), color=clr, label=sample_name, s=20)

    # swarm = plt.scatter(x=list(genomic_regions), y=list(true_n_strains), color=clr, s=50, label=sample_name)

    # # move the swarm dots a bit to the left on x-axis
    # for i in range(len(swarm.collections)):
    #     swarm.collections[i].set_offsets(swarm.collections[i].get_offsets() + np.array([i+offset, 0]))
    
    return swarm

def main(): 
    grs_high_coverage = get_grs_with_sufficient_coverage(["01_100", "02_100", "08_0", "09_0"], low_limit=100)
    colors = sns.color_palette("Set2", 4)
    input_08 = "Experiments/lumc_subsample/08_0/consensus_lumc_comparison_post_processed.tsv"

    df_08 = pd.read_csv(input_08, sep="\t", header=0)
    sample_name = "BA.1: 100%"
    n_strains_08 = process_df(df_08)
    offset_080 = -0.3
    swarm = plot_swarmplot(n_strains_08, sample_name, colors[0], offset_080, grs_high_coverage["08_0"])

    input_01100 = "Experiments/lumc_subsample/01_100/consensus_lumc_comparison_post_processed.tsv"
    df_01100 = pd.read_csv(input_01100, sep="\t", header=0)
    sample_name = "Wuhan: 100%"
    n_strains_01100 = process_df(df_01100)
    offset_01100 = -0.10
    swarm = plot_swarmplot(n_strains_01100, sample_name, colors[1], offset_01100, grs_high_coverage["01_100"])

    input_02100 = "Experiments/lumc_subsample/02_100/consensus_lumc_comparison_post_processed.tsv"
    df_02100 = pd.read_csv(input_02100, sep="\t", header=0)
    sample_name = "Wuhan: 100%"
    n_strains_02100 = process_df(df_02100)
    offset_02100 = 0.3
    swarm = plot_swarmplot(n_strains_02100, sample_name, colors[3], offset_02100, grs_high_coverage["02_100"])

    input_090 = "Experiments/lumc_subsample/09_0/consensus_lumc_comparison_post_processed.tsv"
    df_090 = pd.read_csv(input_090, sep="\t", header=0)
    sample_name = "BA.1: 100%"
    n_strains_090 = process_df(df_090)
    offset_090 = 0.15
    swarm = plot_swarmplot(n_strains_090, sample_name, colors[2], offset_090, grs_high_coverage["09_0"])
    
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
    plt.gca().spines['bottom'].set_position(('data', 0))
    # xticks
    xticks = [region.replace("_", "-")for region in genomic_regions]
    # rotate xticks
    plt.xticks(ticks=np.arange(len(xticks)), labels=xticks, rotation=90)
    # y ticks   
    plt.yticks(np.arange(0, 5, 1))
    plt.xlabel("Genomic regions")
    plt.ylabel("Number of haplotypes")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.135), ncol=3, fontsize=7)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.grid(axis='x', linestyle='--', alpha=0.5)
    plt.ylim(0, 3)
    # y-axis should start from 0 but a bit more up to have a better view
    
    plt.tight_layout()
    plt.savefig(f"Experiments/lumc_subsample/n_strains_pure_v2.pdf", bbox_inches='tight')


if __name__ == "__main__":
    main()