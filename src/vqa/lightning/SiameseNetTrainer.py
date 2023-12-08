import pytorch_lightning as pl
import torch.optim as optim
import torch
import wandb
import vqa.models.cluster_embeddings as cluster_embeddings
import torch.nn.functional as F
import torch.nn as nn
import numpy as npz

class SiameseNetTrainer(pl.LightningModule):
    def __init__(self, model, train_datal, val_datal, test_datal, criterion, optimizer, scheduler):
        super().__init__()
        self.model = model
        self.criterion = criterion
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.train_datal = train_datal
        self.test_datal = test_datal
        self.val_datal = val_datal
    
    def forward(self, input1, input2):
        output1 = self.model(input1)
        output2 = self.model(input2)
        # Compute the distance between the anchor and the unknown, both of shape (batch size, embedding size)
        return output1, output2

  
    def training_step(self, batch, batch_idx):
        (x0, x1) , y = batch
        output1, output2 = self.forward(x0, x1)
        loss = self.criterion(output1, output2, y.to(torch.float))
        euclidean_distance = F.pairwise_distance(output1, output2, keepdim = True)
        predicted_labels = torch.where(euclidean_distance > self.criterion.margin, 0, 1)
        accuracy = torch.mean((predicted_labels == y.to(torch.float)).to(torch.float))
        wandb.log({"train/loss": loss})
        wandb.log({"train/epoch": self.current_epoch})
        wandb.log({"train/accuracy": accuracy})
        return loss

    def validation_step(self, batch, batch_idx):
        (x0, x1) , y = batch
        output1, output2 = self.forward(x0, x1)
        loss = self.criterion(output1, output2, y.to(torch.float))
        euclidean_distance = F.pairwise_distance(output1, output2, keepdim = True)
        predicted_labels = torch.where(euclidean_distance > self.criterion.margin, 0, 1)
        accuracy = torch.mean((predicted_labels == y.to(torch.float)).to(torch.float))
        wandb.log({"val/loss": loss})
        wandb.log({"val/epoch": self.current_epoch})
        wandb.log({"val/accuracy": accuracy})
        # predicted_labels , n_clusters_, homogeneity, completeness  = cluster_embeddings.cluster_embeddings_kmeans(2,output, y)
        # wandb.log({"val/n_clusters": n_clusters_})
        # wandb.log({"val/homogeneity": homogeneity})
        # wandb.log({"val/completeness": completeness})
        return loss
     

    def test_step(self, batch, batch_idx):
        (x0, x1) , y = batch
        output1, output2 = self.forward(x0, x1)
        loss = self.criterion(output1, output2, y.to(torch.float ))
        euclidean_distance = F.pairwise_distance(output1, output2, keepdim = True)
        predicted_labels = torch.where(euclidean_distance > self.criterion.margin, 0, 1)
        accuracy = torch.mean((predicted_labels == y.to(torch.float)).to(torch.float))
        wandb.log({"test/loss": loss})
        wandb.log({"test/epoch": self.current_epoch})
        wandb.log({"test/accuracy": accuracy})
        return loss
        
    
    def configure_optimizers(self):
        # optimizer = optim.Adam(self.model.parameters(), lr=2e-5)
        return self.optimizer
    
    def train_dataloader(self):
        return self.train_datal
    
    def test_dataloader(self):
        return self.test_datal
    
    def val_dataloader(self):
        return self.val_datal
    
