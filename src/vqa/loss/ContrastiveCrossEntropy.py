import torch
import torch.nn.functional as F
import torch.nn as nn

class ContrastiveLoss(torch.nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, output1, output2, label):
        print(output1.shape, output2.shape, label.shape)
        bce_loss = nn.BCELoss()
        m = nn.Sigmoid()
        l1_distance = F.pairwise_distance(output1, output2, keepdim = True, p=1.0)
        # apply sigmoid
        negative_class_prob = m(l1_distance)
        loss = bce_loss(negative_class_prob, label.unsqueeze(1))
        # if the negative class probability is greater than 0.5, then the model predicts the negative class so predicted label is 1 else 0
        predicted_label = torch.where(negative_class_prob > 0.5, torch.tensor([1.0]).to('cuda:0'), torch.tensor([0.0]).to('cuda:0'))
        # calculate accuracy
        return loss, predicted_label
