import torch
import torch.nn as nn


############################################################################################################################################################################
# Input: max_len (int), read (str)
# Returns: padded read. If read is longer than max_len, then truncate the read.
############################################################################################################################################################################
def pad_read(max_len: int, read: str, pre_or_post: str):
    if len(read) < max_len:
        if pre_or_post == "post":
            read += "N" * (max_len - len(read))
        elif pre_or_post == "pre":
            read = "N" * (max_len - len(read)) + read

    if len(read) > max_len:
        read = read[:max_len]    
    return read

############################################################################################################################################################################
# Input: read (str)
# Returns: one hot encoded read
############################################################################################################################################################################
def one_hot_encode(read: str, iupac):
    # A = [1, 0, 0, 0]
    # C = [0, 1, 0, 0]
    # G = [0, 0, 1, 0]
    # T = [0, 0, 0, 1]
    # N = [0, 0, 0, 0]

    base_to_encoding_dict = {"A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T": [0, 0, 0, 1], "N": [0, 0, 0, 0]}
    
    if iupac:
        base_to_encoding_dict = {"A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T": [0, 0, 0, 1], "N": [0, 0, 0, 0], "R": [0.5, 0, 0.5, 0], "Y": [0, 0.5, 0, 0.5], "S": [0, 0.5, 0.5, 0], "W": [0.5, 0, 0, 0.5], "K": [0, 0, 0.5, 0.5], "M": [0.5, 0.5, 0, 0], "B": [0, 0.33, 0.33, 0.33], "D": [0.33, 0, 0.33, 0.33], "H": [0.33, 0.33, 0, 0.33], "V": [0.33, 0.33, 0.33, 0]}
    
    encoded_read = []
    for base in read:
        try:
            encoded_read.append(base_to_encoding_dict[base])
        except KeyError:
            encoded_read.append(base_to_encoding_dict["N"])
    
    # to tensor
    encoded_read = torch.tensor(encoded_read).transpose(0, 1).float()
    return encoded_read


class PadNOneHot(nn.Module):
    def __init__(self, max_len: int, pre_or_post: str, single_read: bool = False, iupac: bool = False, pad: bool = True):
        super(PadNOneHot, self).__init__()
        self.max_len = max_len
        self.iupac = iupac
        self.pad = pad
        self.pre_or_post = pre_or_post
        self.single_read = single_read

    def forward(self, data):
        if self.single_read == False:
            read_1, read_2 = data
            if self.pad:
                read_1 = pad_read(self.max_len, read_1, self.pre_or_post)
                read_2 = pad_read(self.max_len, read_2, self.pre_or_post)
            read_1 = one_hot_encode(read_1, iupac=self.iupac)
            read_2 = one_hot_encode(read_2, iupac=self.iupac)
            return read_1, read_2
     
        read_1 = data
        if self.pad:
            read_1 = pad_read(self.max_len, read_1, self.pre_or_post)
        read_1 = one_hot_encode(read_1, iupac=self.iupac)

        return read_1

        