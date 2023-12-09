import torch
import torch.nn.functional as F
import torch.nn as nn

class ContrastiveLoss(torch.nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, output1, output2, label):

        # lstm is bidirectional so the output is of shape (d, bach_size ,hidden_size) where d is 2 for bidirectional lstm
        # output1 =torch.mean(output1, dim=0)
        # output2 =torch.mean(output2, dim=0)
        bce_loss = nn.BCELoss()
        l1_distance = F.pairwise_distance(output1, output2, keepdim = True, p=1.0)
        m = nn.Sigmoid()
        l1_distance = m(l1_distance)
        loss = bce_loss(l1_distance.unsqueeze(1), label.unsqueeze(1))
        
        # if the negative class probability is greater than 0.5, then the model predicts the negative class so predicted label is 1 else 0
        predicted_labels = torch.where(l1_distance > 0.5, torch.tensor([1]).to('cpu'), torch.tensor([0]).to('cpu'))
        # calculate accuracy
        return loss, predicted_labels
