import torch
import torch.nn.functional as F

class TripletLoss(torch.nn.Module):
      
     def __init__(self, margin):
          super().__init__()
          self.margin = margin
    
     def forward(self, anchor, positive, negative):
          distance_positive = self.cosine_distance(anchor, positive)
          distance_negative = selfcosine_distance(anchor, negative)
          loss = torch.mean(torch.clamp(distance_positive - distance_negative + self.margin, min=0.0))
          return loss
     
     def cosine_distance(self, input1, input2):
          distance = 1-  F.cosine_similarity(input1, input2, dim=1)
          return distance