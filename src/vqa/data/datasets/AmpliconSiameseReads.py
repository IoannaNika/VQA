from torch.utils.data import Dataset
import json
import os
import pandas as pd
import torch

class SiameseReads(Dataset):
    def __init__(self, directory: str, transform=None, genomic_region=None):
        self.directory = directory
        self.transform = transform
        self.genomic_region = genomic_region
        # read tsv
        self.reference_set = pd.read_csv(os.path.join(self.directory, 'lumc_dataset.tsv'), sep='\t', header=0)
        if self.genomic_region != None:
            try:
                start = self.genomic_region.split("_")[0]
                end = self.genomic_region.split("_")[1]
                self.reference_set = self.reference_set[self.reference_set["start"] == int(start)]
            except:
                self.reference_set = self.reference_set[self.reference_set["genomic_regions"] == self.genomic_region]

        print("Positives: ", len(self.reference_set[self.reference_set["label"] == "positive"]))
        print("Negatives: ", len(self.reference_set[self.reference_set["label"] != "positive"]))
        self.length = len(self.reference_set)
        if self.length > 0:
            print("Accuracy if it predicts all positive: ", len(self.reference_set[self.reference_set["label"] == "positive"])/len(self.reference_set))

        # edit distance to list 
        # self.edit_distances_pos = self.reference_set[self.reference_set["label"] == "positive"]["edit_distance"].to_list()
        # self.edit_distances_neg = self.reference_set[self.reference_set["label"] != "positive"]["edit_distance"].to_list()

        # print(self.edit_distances_pos)
        # save as list  to file that if it exists append to it otherwise create it
        # if os.path.exists("edit_distances_positives.txt"):
        #     with open("edit_distances_positives.txt", "a") as f:
        #         f.write("\n{}\n".format(self.genomic_region))
        #         f.write(str(self.edit_distances_pos))
        # else:
        #     with open("edit_distances_positives.txt", "w") as f:
        #         f.write("\n{}\n".format(self.genomic_region))
        #         f.write(str(self.edit_distances_pos))

        # if os.path.exists("edit_distances_negatives.txt"):
        #     with open("edit_distances_negatives.txt", "a") as f:
        #         f.write("\n{}\n".format(self.genomic_region))
        #         f.write(str(self.edit_distances_neg))
        # else:
        #     with open("edit_distances_negatives.txt", "w") as f:
        #         f.write("\n{}\n".format(self.genomic_region))
        #         f.write(str(self.edit_distances_neg))


    def __getitem__(self, index: int):

        items = self.reference_set.iloc[index]
        # print(items)
        label = items["label"]
        read1 = items["read1"]
        read2 = items["read2"]
        
        target = 1 if label == 'positive' else 0

        data = (read1, read2)

        if self.transform:
            data = self.transform(data)

        return data, target
    
    def __len__(self):
        return self.length
