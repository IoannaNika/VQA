import pytorch_lightning as pl
import torch.optim as optim
import torch
import wandb
import vqa.models.cluster_embeddings as cluster_embeddings
import torch.nn.functional as F
import torch.nn as nn
import numpy as npz

class TripletNetTrainer(pl.LightningModule):
    def __init__(self, model, train_datal, val_datal, test_datal,optimizer):
        super().__init__()
        self.model = model
        self.margin = 0.3
        self.criterion = nn.TripletMarginLoss(margin=self.margin, p=2, eps=1e-7)
        self.optimizer = optimizer
        self.train_datal = train_datal
        self.test_datal = test_datal
        self.val_datal = val_datal

    def forward(self, input1, input2, input3):
        # change second dimension to the last dimension
        output1 = self.model(input1) # hidden dimension only
        output2 = self.model(input2) # hidden dimension only
        output3 = self.model(input3) # hidden dimension only
        # Compute the distance between the anchor and the unknown, both of shape (batch size, embedding size)
        return output1, output2, output3

  
    def training_step(self, batch, batch_idx):
        x0, x1, x2 = batch
        output1, output2, output3 = self.forward(x0, x1, x2)
        # anchor, positive, negative
        loss = self.criterion(output1, output2, output3)
        accuracy = self.get_accuracy(self.margin, output1, output2, output3)
        wandb.log({"epoch": self.current_epoch, "train/loss": loss})
        wandb.log({"epoch": self.current_epoch, "train/accuracy": accuracy})
        return loss

    def validation_step(self, batch, batch_idx):
        x0, x1, x2 = batch
        # anchor, positive, negative
        output1, output2, output3 = self.forward(x0, x1, x2)
        loss = self.criterion(output1, output2, output3)
        accuracy = self.get_accuracy(self.margin, output1, output2, output3)
        self.log("val_loss", loss,  batch_size = 20) 
        self.log("val_acc", accuracy,  batch_size = 20)
        print("validation accuracy: ", accuracy)
        wandb.log({"epoch": self.current_epoch, "val/loss": loss})
        wandb.log({"epoch": self.current_epoch, "val/accuracy": accuracy})
        return loss
    

    def test_step(self, batch, batch_idx):
        x0, x1, x2 = batch
        output1, output2, output3 = self.forward(x0, x1, x2)
        # anchor, positive, negative
        loss = self.criterion(output1, output2, output3)
        accuracy = self.get_accuracy(self.margin, output1, output2, output3)
        wandb.log({"epoch": self.current_epoch, "val/loss": loss})
        wandb.log({"epoch": self.current_epoch, "val/accuracy": accuracy})
        return loss
      
    def get_accuracy(self, margin, output1, output2, output3):
        euclidean_distance_pos = F.pairwise_distance(output1, output2, keepdim = True) 
        euclidean_distance_neg = F.pairwise_distance(output1, output3, keepdim = True) 

        predicted_labels_pos = torch.where(euclidean_distance_pos > self.margin, 0, 1)
        predicted_labels_neg = torch.where(euclidean_distance_neg > self.margin, 0, 1)

        accuracy_pos = torch.sum(predicted_labels_pos == 1).item() / len(predicted_labels_pos)
        accuracy_neg = torch.sum(predicted_labels_neg == 0).item() / len(predicted_labels_neg)
        accuracy = (accuracy_pos + accuracy_neg) / 2
        return accuracy
    
    def configure_optimizers(self):
        return self.optimizer
    
    def train_dataloader(self):
        return self.train_datal
    
    def test_dataloader(self):
        return self.test_datal
    
    def val_dataloader(self):
        return self.val_datal
    
