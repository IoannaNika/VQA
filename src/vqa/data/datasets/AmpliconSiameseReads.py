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
