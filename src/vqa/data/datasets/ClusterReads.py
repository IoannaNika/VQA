from torch.utils.data import Dataset
import json
import os

class ClusterReads(Dataset):
    def __init__(self, directory: str, transform=None):
        self.directory = directory
        self.transform = transform
       
       # get list of .fasta files in the directory
        self.reference_set = []
        for file in os.listdir(directory):
            if file.endswith(".fasta"):
                self.reference_set.append(file)
        self.length = len(self.reference_set)

    def __getitem__(self, index: int):

        read_file = self.reference_set[index]


        read_path = os.path.join(self.directory, read_file)


        with open(read_path, 'r') as read_f:
            read = read_f.readlines()[1].strip()
        
        data = read
        target = read_file.split("_")[:-1]
        target = "_".join(target)

        if self.transform:
            data = self.transform(data)

        return data, target
    
    def __len__(self):
        return self.length
     