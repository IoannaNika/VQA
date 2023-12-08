import torch
import torch.nn.functional as F
import torch.nn as nn

class ContrastiveLoss(torch.nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, output1, output2, label):
        bce_loss = nn.BCELoss()
        l1_distance = F.pairwise_distance(output1, output2, keepdim = True, p=1.0)
        m = nn.Sigmoid()
        l1_distance = m(l1_distance)
        loss = bce_loss(l1_distance, label.unsqueeze(1))
        
        # if the negative class probability is greater than 0.5, then the model predicts the negative class so predicted label is 1 else 0
        predicted_labels = torch.where(l1_distance > 0.5, torch.tensor([1]).to('gpu'), torch.tensor([0]).to('gpu'))
        # calculate accuracy
        return loss, predicted_labels
