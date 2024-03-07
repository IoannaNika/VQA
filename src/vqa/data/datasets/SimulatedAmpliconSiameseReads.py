from torch.utils.data import Dataset
import json
import os
import pandas as pd
import torch
import editdistance

class SiameseReads(Dataset):
    def __init__(self, directory: str, transform=None, genomic_region=None, test_mode=False):
        self.directory = directory
        self.transform = transform
        self.genomic_region = genomic_region
        self.test_mode = test_mode
       
        # read tsv file
        self.reference_set = pd.read_csv(self.directory + "/test_pairs.tsv", sep='\t', header=0)
        # self.reference_set = self.reference_set[self.reference_set["label"] != "positive"]
        # print(self.reference_set)
        if self.genomic_region != None:
            self.reference_set = self.reference_set[self.reference_set["genomic_region"] == self.genomic_region]
        
        # print("Positives: ", len(self.reference_set[self.reference_set["label"] == "positive"]))
        # print("Negatives: ", len(self.reference_set[self.reference_set["label"] != "positive"]))
        self.length = len(self.reference_set)
        # if self.length > 0:
        #     print("Accuracy if it predicts all positive: ", len(self.reference_set[self.reference_set["label"] == "positive"])/len(self.reference_set))


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
            
        target = 1 if label == 'positive' or label == 1 else 0

        data = (fasta_1, fasta_2)

        # ed = editdistance.eval(fasta_1, fasta_2)
        # print("data info: ", ed, target)

        if self.transform:
            data = self.transform(data)
        
        if self.test_mode: 
            ids = (read_id_1, read_id_2)
            data = (fasta_1, fasta_2)
            gr = items["genomic_region"]
            return  ids, data, target, gr

        return data, target
    
    def __len__(self):
        return self.length
