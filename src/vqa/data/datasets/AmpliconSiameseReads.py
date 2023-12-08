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
        self.reference_set = pd.read_csv(os.path.join(self.directory, 'amplicon_lumc_pairs.tsv'), sep='\t', header=None)
        # print(self.reference_set)
        self.length = len(self.reference_set)

    def __getitem__(self, index: int):

        items = self.reference_set.iloc[index]
        if items[0].startswith('m'):
            read1 = items[0]
            if items[1].startswith('m'):
                read2 = items[1]
                label = items[2]
            else:
                read2 = items[2]
                label = items[1]

        else:
            label = items[0]
            read1 = items[1]
            read2 = items[2]
        
        fasta_1_path = os.path.join(self.directory, '{}.fastq'.format(read1))
        fasta_2_path = os.path.join(self.directory, '{}.fastq'.format(read2))


        with open(fasta_1_path, 'r') as fasta_1_file:
            fasta_1 = fasta_1_file.readlines()[1].strip()
        
        
        with open(fasta_2_path, 'r') as fasta_2_file:
            fasta_2 = fasta_2_file.readlines()[1].strip()
            
        target = 0 if label == 'positive' else 1

        data = (fasta_1, fasta_2)

        if self.transform:
            data = self.transform(data)




        return data, target
    
    def __len__(self):
        return self.length
