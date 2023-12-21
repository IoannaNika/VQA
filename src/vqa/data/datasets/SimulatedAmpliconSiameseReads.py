from torch.utils.data import Dataset
import json
import os
import pandas as pd
import torch

class SiameseReads(Dataset):
    def __init__(self, directory: str, transform=None):
        self.directory = directory
        self.transform = transform
       
        # read tsv file
        self.reference_set = pd.read_csv(self.directory + "/samples.tsv", sep='\t', header=0)
        # print(self.reference_set)
        self.length = len(self.reference_set)

    def __getitem__(self, index: int):

        items = self.reference_set.iloc[index]
        read_id_1 = items["read_1"]
        read_id_2 = items["read_2"]

        label = items["label"]
        
        fasta_1_path = os.path.join(self.directory,"reads",  '{}.fasta'.format(read_id_1))
        fasta_2_path = os.path.join(self.directory, "reads", '{}.fasta'.format(read_id_2))
        
        with open(fasta_1_path, 'r') as fasta_1_file:
            fasta_1 = fasta_1_file.readlines()[1].strip()
        
        
        with open(fasta_2_path, 'r') as fasta_2_file:
            fasta_2 = fasta_2_file.readlines()[1].strip()
            
        target = 1 if label == 'positive' else 0

        data = (fasta_1, fasta_2)

        if self.transform:
            data = self.transform(data)

        return data, target
    
    def __len__(self):
        return self.length
