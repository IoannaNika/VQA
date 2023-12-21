from torch.utils.data import Dataset
import json
import os
import pandas as pd
import torch

class SiameseReads(Dataset):
    def __init__(self, directory: str, transform=None):
        self.directory = directory
        self.transform = transform
       
        # read tsv
        self.reference_set = pd.read_csv(os.path.join(self.directory, 'lumc_dataset.tsv'), sep='\t', header=0)
        # print(self.reference_set)
        self.length = len(self.reference_set)

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
