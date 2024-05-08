import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def gt(): 
    return {"03_50": "Omicron", "04_75":"Wuhan", "05_90":"Omicron", "07_98": "Wuhan"}

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
        expected_n_strains[region] = 2
        if region in ["1128_2244", "3166_4240","15634_16698", "16647_17732", "19604_20676", "29768_29790"]:
            expected_n_strains[region] = 1
    return expected_n_strains

def plot_swarmplot(n_strains, sample_name, clr, offset):
    # n_strains is a dictionary with the genomic region as key and the number of strains as value
    genomic_regions = n_strains.keys()
    true_n_strains = n_strains.values()

    x = [i + offset for i in range(len(genomic_regions))]
    swarm = plt.scatter(x, list(true_n_strains), color=clr, label=sample_name, s=20)

    # plot swarmplot the x-axis is the genomic region and the y-axis is the number of strains
    # swarm = sns.scatterplot(x=list(genomic_regions), y=list(true_n_strains), color=clr, s=50, label=sample_name)
    # set legend label for the specific swarmplot
    # swarm.set_label(sample_name)
    # # show label for the swarmplot
    # swarm.legend()

    # move the swarm dots a bit to the left on x-axis
    # for i in range(len(swarm.collections)):
    #     swarm.collections[i].set_offsets(swarm.collections[i].get_offsets() + np.array([i+offset, 0]))
    return swarm
def main(): 
    colors = sns.color_palette("Set2", 4)
    input_03 = "Experiments/lumc/mixtures/03_50/consensus_lumc_comparison.tsv"

    df_03 = pd.read_csv(input_03, sep="\t", header=0)
    sample_name = "03_50"
    gt_sample = gt()[sample_name]
    results_03, n_strains_03 = process_df(df_03, gt_sample)
    swarm = plot_swarmplot(n_strains_03, sample_name, colors[0], -0.2)

    input_04 = "Experiments/lumc/mixtures/04_75/consensus_lumc_comparison.tsv"

    df_04 = pd.read_csv(input_04, sep="\t", header=0)
    sample_name = "04_75"
    gt_sample = gt()[sample_name]
    results_04, n_strains_04 = process_df(df_04, gt_sample)
    swarm = plot_swarmplot(n_strains_04, sample_name, colors[1], -0.05)

    input_05 = "Experiments/lumc/mixtures/05_90/consensus_lumc_comparison.tsv"

    df_05 = pd.read_csv(input_05, sep="\t", header=0)
    sample_name = "05_90"
    gt_sample = gt()[sample_name]
    results_05, n_strains_05 = process_df(df_05, gt_sample)
    swarm = plot_swarmplot(n_strains_05, sample_name, colors[2], 0.05)

    input_07 = "Experiments/lumc/mixtures/07_98/consensus_lumc_comparison.tsv"
    
    df_07 = pd.read_csv(input_07, sep="\t", header=0)
    sample_name = "07_98"
    gt_sample = gt()[sample_name]
    results_07, n_strains_07 = process_df(df_07, gt_sample)
    swarm = plot_swarmplot(n_strains_07, sample_name, colors[3], 0.2)


    genomic_regions = ["54_1183", "1128_2244", "2179_3235", "3166_4240", "4189_5337", "5286_6358", "6307_7379", "7328_8363", "8282_9378", "9327_10429", "10370_11447", "11394_12538", "12473_13599", "13532_14619", "14568_15713", "15634_16698", "16647_17732", "17649_18684", "18618_19655", "19604_20676", "21562_22590", "24658_25768", "25712_26835", "26766_27872", "27808_28985", "28699_29768"]



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
    plt.subplots_adjust(bottom=0.7)
    plt.tight_layout()
    plt.savefig(f"Experiments/lumc/mixtures/n_strains.pdf", bbox_inches='tight')







if __name__ == "__main__":
    main()