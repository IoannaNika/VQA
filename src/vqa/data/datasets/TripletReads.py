from torch.utils.data import Dataset
import json
import os
import pandas as pd

class TripletReads(Dataset):
    def __init__(self, directory: str, transform=None):
        self.directory = directory
        self.transform = transform
        self.reference_set = pd.read_csv(os.path.join(self.directory, 'triplets.tsv'), sep='\t', header=0)
        self.length = len(self.reference_set)

    def __getitem__(self, index: int):

        row = self.reference_set.iloc[index]

        read_path_anch = os.path.join(self.directory, "reads", row["read_anch"] + ".fasta")
        read_path_pos = os.path.join(self.directory, "reads", row["read_pos"] + ".fasta")
        read_path_neg = os.path.join(self.directory, "reads", row["read_neg"] + ".fasta")


        with open(read_path_anch, 'r') as read_f:
            read_anch = read_f.readlines()[1].strip()

        with open(read_path_pos, 'r') as read_f:
            read_pos = read_f.readlines()[1].strip()
        
        with open(read_path_neg, 'r') as read_f:
            read_neg = read_f.readlines()[1].strip()

        if self.transform:
           data = (self.transform(read_anch), self.transform(read_pos), self.transform(read_neg))
        else: 
            data = (read_anch, read_pos, read_neg)
        
        return data
    
    def __len__(self):
        return self.length
     