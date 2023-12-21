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
        # change second dimension to the last dimension
        # input1 = input1.transpose(1,2)
        # input2 = input2.transpose(1,2)
        output1 = self.model(input1) # hidden dimension only
        output2 = self.model(input2) # hidden dimension only
        output1 = output1.squeeze(0)
        output2 = output2.squeeze(0)

        # Compute the distance between the anchor and the unknown, both of shape (batch size, embedding size)
        return output1, output2

  
    def training_step(self, batch, batch_idx):
        (x0, x1) , y = batch
        output1, output2 = self.forward(x0, x1)
        # print("output1: ", output1, "output2: ", output2, "y: ", y)
        loss, predicted_labels = self.criterion(output1, output2, y.to(torch.float))
        # predicted_labels = predicted_labels.squeeze(1)
        accuracy = torch.mean((predicted_labels == y.to(torch.float)).to(torch.float))
        wandb.log({"epoch": self.current_epoch, "train/loss": loss})
        wandb.log({"epoch": self.current_epoch, "train/accuracy": accuracy})
        return loss

    # def validation_step(self, batch, batch_idx):
    #     (x0, x1) , y = batch
    #     output1, output2 = self.forward(x0, x1)
    #     print("output1: ", output1, "output2: ", output2)
    #     loss, predicted_labels = self.criterion(output1, output2, y.to(torch.float))
    #     accuracy = torch.mean((predicted_labels == y.to(torch.float)).to(torch.float))
    #     wandb.log({"epoch": self.current_epoch, "val/loss": loss})
    #     wandb.log({"epoch": self.current_epoch, "val/accuracy": accuracy})
    #     return loss
     
    def validation_step(self, batch, batch_idx):
        x, y = batch
        # x = x.transpose(1,2)
        output = self.model(x)
        print("output: ", output)
        output = output.squeeze(0)
        predicted_labels , n_clusters_, homogeneity, completeness  = cluster_embeddings.cluster_embeddings_dbscan(output, y)
        print("predicted_labels" , predicted_labels)
        wandb.log({"epoch": self.current_epoch, "val/n_clusters": n_clusters_})
        wandb.log({"epoch": self.current_epoch, "val/homogeneity": homogeneity})
        wandb.log({"epoch": self.current_epoch, "val/completeness": completeness})
        return homogeneity

    def test_step(self, batch, batch_idx):
        (x0, x1) , y = batch
        output1, output2 = self.forward(x0, x1)
        loss, predicted_labels = self.criterion(output1, output2, y.to(torch.float ))
        # predicted_labels = predicted_labels.squeeze(1)
        accuracy = torch.mean((predicted_labels == y.to(torch.float)).to(torch.float))
        wandb.log({"test/loss": loss})
        wandb.log({"test/accuracy": accuracy})
        return loss
        
    
    def configure_optimizers(self):
        return self.optimizer
    
    def train_dataloader(self):
        return self.train_datal
    
    def test_dataloader(self):
        return self.test_datal
    
    def val_dataloader(self):
        return self.val_datal
    
