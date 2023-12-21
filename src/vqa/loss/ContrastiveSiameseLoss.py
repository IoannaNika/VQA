import torch
import torch.nn.functional as F

class ContrastiveLoss(torch.nn.Module):
  
    def __init__(self, margin):
        super().__init__()
        self.margin = margin

    def forward(self, output1, output2, label):
        euclidean_distance = F.pairwise_distance(output1, output2, keepdim = True)
        # print("output1: ", output1.shape ,"\noutput2: ", output2.shape)
  
        try:
            euclidean_distance = euclidean_distance.squeeze(1)
        except:
            euclidean_distance = euclidean_distance
            
        loss_contrastive = torch.mean(0.5 * (label) * torch.pow(euclidean_distance, 2) + 0.5 * (1-label) * torch.pow(torch.clamp(self.margin - euclidean_distance, min=0.0), 2))
        # print("loss: ", loss_contrastive)
        predicted_labels = torch.where(euclidean_distance > self.margin, 0, 1)
        # print("predicted labels: ", predicted_labels)
        return loss_contrastive, predicted_labels

    