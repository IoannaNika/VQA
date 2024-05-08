import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def gt(): 
    return {"08_0": "Omicron", "09_0": "Omicron", "01_100": "Wuhan", "02_100": "Wuhan"}

def init_results(genomic_regions):
   
    results = {}
    for region in genomic_regions:

        results[region] = {}
        results[region]["Wuhan"] = {}
        results[region]["Omicron"] = {}
        results[region]["Wuhan"]["Relative_abundance"]  = 0
        results[region]["Omicron"]["Relative_abundance"]  = 0
        results[region]["Wuhan"]["Edit_distance"]  = []
        results[region]["Omicron"]["Edit_distance"]  = []
    return results

def init_n_strains(genomic_regions):
    n_strains = {}
    for region in genomic_regions:
        n_strains[region] = 0

    return n_strains

def process_df(df, gt):
    genomic_regions = df["Genomic_region"].unique()

    results = init_results(genomic_regions)
    n_strains = init_n_strains(genomic_regions)


    for index, row in df.iterrows():
        region = row["Genomic_region"]
        strain = row["Min_consensus"]
        if strain == "BA.1":
            strain = "Omicron"
        if strain == "Both":
            strain = gt
        relative_abundance = row["Relative_abundance"]
        edit_distance = row["Min_edit_distance"]

        results[region][strain]["Relative_abundance"] += relative_abundance
        results[region][strain]["Edit_distance"].append(edit_distance)
        n_strains[region] += 1
    
    return results, n_strains

def get_expected_n_of_strains(genomic_regions):
    expected_n_strains = {}
    for region in genomic_regions:
        expected_n_strains[region] = 1
    return expected_n_strains

def plot_swarmplot(n_strains, sample_name, clr, offset):
    # n_strains is a dictionary with the genomic region as key and the number of strains as value
    genomic_regions = n_strains.keys()
    true_n_strains = n_strains.values()
    # plot swarmplot the x-axis is the genomic region and the y-axis is the number of strains
    x = [i + offset for i in range(len(genomic_regions))]
    swarm = plt.scatter(x, list(true_n_strains), color=clr, label=sample_name, s=20)

    # swarm = plt.scatter(x=list(genomic_regions), y=list(true_n_strains), color=clr, s=50, label=sample_name)

    # # move the swarm dots a bit to the left on x-axis
    # for i in range(len(swarm.collections)):
    #     swarm.collections[i].set_offsets(swarm.collections[i].get_offsets() + np.array([i+offset, 0]))
    
    return swarm
def main(): 
    colors = sns.color_palette("Set2", 4)
    input_08 = "Experiments/lumc/mixtures/08_0/consensus_lumc_comparison.tsv"

    df_08 = pd.read_csv(input_08, sep="\t", header=0)
    sample_name = "08_0"
    gt_sample = gt()[sample_name]
    results_08, n_strains_08 = process_df(df_08, gt_sample)
    offset_080 = -0.3
    swarm = plot_swarmplot(n_strains_08, sample_name, colors[0], offset_080)

    input_01100 = "Experiments/lumc/mixtures/01_100/consensus_lumc_comparison.tsv"
    df_01100 = pd.read_csv(input_01100, sep="\t", header=0)
    sample_name = "01_100"
    results_01100, n_strains_01100 = process_df(df_01100, gt_sample)
    offset_01100 = -0.15
    swarm = plot_swarmplot(n_strains_01100, sample_name, colors[1], offset_01100)

    input_090 = "Experiments/lumc/mixtures/09_0/consensus_lumc_comparison.tsv"
    df_090 = pd.read_csv(input_090, sep="\t", header=0)
    sample_name = "09_0"
    gt_sample = gt()[sample_name]
    results_090, n_strains_090 = process_df(df_090, gt_sample)
    offset_090 = 0.15
    swarm = plot_swarmplot(n_strains_090, sample_name, colors[2], offset_090)
    genomic_regions = df_08["Genomic_region"].unique()

    expected_n_strains = get_expected_n_of_strains(genomic_regions)

    
    # plot expected number of strains
    sns.lineplot(list(expected_n_strains.values()), linestyle ='--', label="Expected number of strains", color="black", linewidth=0.5,  marker=".", markersize=7)
    # plt.setp(swarm.lines, zorder=1000)

    # xticks
    xticks = [region.replace("_", "-")for region in genomic_regions]
    # rotate xticks
    plt.xticks(ticks=np.arange(len(xticks)), labels=xticks, rotation=90)
    plt.ylim(-0.5, 5)
    # y ticks   
    plt.yticks(np.arange(0, 5, 1))
    plt.xlabel("Genomic region")
    plt.ylabel("Number of haplotypes")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.10), ncol=5, fontsize=6)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.grid(axis='x', linestyle='--', alpha=0.5)
    # y-axis should start from 0 but a bit more up to have a better view
    
    plt.tight_layout()
    plt.savefig(f"Experiments/lumc/mixtures/n_strains_pure.pdf", bbox_inches='tight')







if __name__ == "__main__":
    main()