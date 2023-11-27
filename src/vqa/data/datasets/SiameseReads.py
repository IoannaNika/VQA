from torch.utils.data import Dataset
import json
import os

class SiameseReads(Dataset):
    def __init__(self, directory: str, transform=None):
        self.directory = directory
        self.transform = transform
       
       # read reference_set.json into a dictionary
        with open(os.path.join(self.directory, 'reference_set.json')) as json_file:
            self.reference_set = list(json.load(json_file).items())
        self.length = len(self.reference_set)

    def __getitem__(self, index: int):

        key, value = self.reference_set[index]

        fasta_1_id = value['read1']
        fasta_2_id = value['read2']

        fasta_1_path = os.path.join(self.directory, 'reads', '{}.fasta'.format(fasta_1_id))
        fasta_2_path = os.path.join(self.directory, 'reads', '{}.fasta'.format(fasta_2_id))


        with open(fasta_1_path, 'r') as fasta_1_file:
            fasta_1 = fasta_1_file.readlines()[1].strip()
        
        
        with open(fasta_2_path, 'r') as fasta_2_file:
            fasta_2 = fasta_2_file.readlines()[1].strip()

        target = 1 if value['label'] == 'positive' else 0

        data = (fasta_1, fasta_2)

        if self.transform:
            data = self.transform(data)

        return data, target
    
    def __len__(self):
        return self.length
