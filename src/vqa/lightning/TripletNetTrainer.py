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
        self.criterion = nn..TripletMarginLoss(margin=2.0, p=2, eps=1e-7)
        self.optimizer = optimizer
        self.train_datal = train_datal
        self.test_datal = test_datal
        self.val_datal = val_datal

    def forward(self, input1, input2, input3):
        # change second dimension to the last dimension
        input1 = input1.transpose(1,2)
        input2 = input2.transpose(1,2)
        input3 = input3.transpose(1,2)
        output1 = self.model(input1) # hidden dimension only
        output2 = self.model(input2) # hidden dimension only
        output3 = self.model(input3) # hidden dimension only
        output1 = output1.squeeze(0)
        output2 = output2.squeeze(0)
        output3 = output3.squeeze(0)

        # Compute the distance between the anchor and the unknown, both of shape (batch size, embedding size)
        return output1, output2, output3

  
    def training_step(self, batch, batch_idx):
        (x0, x1, x2) , y = batch
        output1, output2, output3 = self.forward(x0, x1, x2)
        # anchor, positive, negative
        loss, predicted_labels, accuracy = self.criterion(output1, output2, output3)
        wandb.log({"epoch": self.current_epoch, "train/loss": loss})
        wandb.log({"epoch": self.current_epoch, "train/accuracy": accuracy})
        return loss

    def validation_step(self, batch, batch_idx):
        (x0, x1, x2) , y = batch
        # anchor, positive, negative
        output1, output2, output3 = self.forward(x0, x1, x2)
        loss = self.criterion(output1, output2, output3)
        self.log("val_loss", loss,  batch_size = 20) 
        wandb.log({"epoch": self.current_epoch, "val/loss": loss})
        wandb.log({"epoch": self.current_epoch, "val/accuracy": accuracy})
        return {"val_loss": loss}
    

    def test_step(self, batch, batch_idx):
        (x0, x1, x2) , y = batch
        output1, output2, output3 = self.forward(x0, x1, x2)
        # anchor, positive, negative
        loss = self.criterion(output1, output2, output3)
        self.log("val_loss", loss,  batch_size = 20) 
        wandb.log({"epoch": self.current_epoch, "val/loss": loss})
        wandb.log({"epoch": self.current_epoch, "val/accuracy": accuracy})
        return {"val_loss": loss}
      
    def get_accuracy(self, margin):
        euclidean_distanc_pos = F.pairwise_distance(output1, output2, keepdim = True) 
        euclidean_distanc_neg = F.pairwise_distance(output1, output3, keepdim = True) 

        predicted_labels_pos = torch.where(euclidean_distance_pos > self.margin, 0, 1)
        predicted_labels_neg = torch.where(euclidean_distance_pos > self.margin, 0, 1)

        #TODO

        return
    
    def configure_optimizers(self):
        return self.optimizer
    
    def train_dataloader(self):
        return self.train_datal
    
    def test_dataloader(self):
        return self.test_datal
    
    def val_dataloader(self):
        return self.val_datal
    
