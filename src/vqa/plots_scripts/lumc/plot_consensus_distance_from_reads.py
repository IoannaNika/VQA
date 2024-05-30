import os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

def main(): 


    genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                            (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                            (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                            (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                            (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                            (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                            (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]
    
    genomic_regions = [str(region[0]) + "_" + str(region[1]) for region in genomic_regions]

    samples = ["01_100", "02_100", "03_50", "04_75", "05_90", "06_95", "07_98", "08_0", "09_0"]

    # load the dictionaries
    with open("Experiments/lumc_subsample/edit_distances_per_sample_consensus.pkl", "rb") as file:
        edit_distances_per_sample_consensus = pickle.load(file)
    with open("Experiments/lumc_subsample/edit_distances_per_sample_contigs.pkl", "rb") as file:
        edit_distances_per_sample_contigs = pickle.load(file)

    # plot the edit distances as boxplots next to each other
    fig, ax = plt.subplots()
    positions_bplot_1 = [x-0.2 for x in range(1, len(samples)+1)]
    positions_bplot_2 = [x+0.2 for x in range(1, len(samples)+1)]
    colors_set_2 = sns.color_palette("Set2")
    bplot_1 = ax.boxplot([edit_distances_per_sample_consensus[sample] for sample in samples], positions=positions_bplot_1, patch_artist=True, showfliers=False, widths=0.2)
    bplot_2 = ax.boxplot([edit_distances_per_sample_contigs[sample] for sample in samples], positions=positions_bplot_2, patch_artist=True, showfliers=False, widths=0.2)

    for patch in bplot_1['boxes']:
        patch.set_facecolor(colors_set_2[2])
    for patch in bplot_2['boxes']:
        patch.set_facecolor(colors_set_2[3])

    ax.set_xticks(range(1, len(samples)+1))
    ax.set_xticklabels(samples)
    ax.set_xlabel("Samples")
    ax.set_ylabel("Edit distance")
    legend = [plt.Rectangle((0,0),1,1,fc=colors_set_2[i]) for i in [2,3]]
    ax.legend(legend, ['True haplotypes', 'Contigs'], loc='upper right')
    plt.ylim(0, 7)
    plt.savefig("Experiments/lumc_subsample/edit_distances_consensus_contigs.pdf", bbox_inches='tight')


if __name__ == '__main__':
    main()