from torch.utils.data import Dataset
import json
import os
import pandas as pd

class ClusterReads(Dataset):
    def __init__(self, directory: str, transform=None, genomic_region = None):
        self.directory = directory
        self.transform = transform
        self.genomic_region = genomic_region
        self.reference_set = pd.read_csv(os.path.join(self.directory, 'single_samples.tsv'), sep='\t', header=0)
        if self.genomic_region != None:
            self.reference_set = self.reference_set[self.reference_set["genomic_region"] == self.genomic_region]
        self.length = len(self.reference_set)

    def __getitem__(self, index: int):

        row = self.reference_set.iloc[index]
        read_path = os.path.join(self.directory, "reads", row["read_id"] + ".fasta")

        with open(read_path, 'r') as read_f:
            read = read_f.readlines()[1].strip()
        
        data = read
        target = row["label"] + "_" + row["genomic_region"]

        if self.transform:
            data = self.transform(data)
        return data, target
    
    def __len__(self):
        return self.length
     