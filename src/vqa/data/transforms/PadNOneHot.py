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
def one_hot_encode(read: str):
    # A = [1, 0, 0, 0]
    # C = [0, 1, 0, 0]
    # G = [0, 0, 1, 0]
    # T = [0, 0, 0, 1]
    # N = [0, 0, 0, 0]

    base_to_encoding_dict = {"A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T": [0, 0, 0, 1], "N": [0, 0, 0, 0]}
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
    def __init__(self, max_len: int, pre_or_post: str, single_read: bool = False):
        super(PadNOneHot, self).__init__()
        self.max_len = max_len
        self.pre_or_post = pre_or_post
        self.single_read = single_read

    def forward(self, data):
        if self.single_read == False:
            read_1, read_2 = data
            read_1 = pad_read(self.max_len, read_1, self.pre_or_post)
            read_2 = pad_read(self.max_len, read_2, self.pre_or_post)
            read_1 = one_hot_encode(read_1)
            read_2 = one_hot_encode(read_2)
            return read_1, read_2
     
        read_1 = data
        read_1 = pad_read(self.max_len, read_1, self.pre_or_post)
        read_1 = one_hot_encode(read_1)

        return read_1

        