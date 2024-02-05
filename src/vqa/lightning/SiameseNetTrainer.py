import pytorch_lightning as pl
import torch.optim as optim
import torch
import wandb
import vqa.models.cluster_embeddings as cluster_embeddings
import torch.nn.functional as F
import torch.nn as nn
import numpy as npz
from torcheval.metrics.functional import multiclass_f1_score


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
    
    def forward(self, input1):
        # change second dimension to the last dimension
        output1 = self.model(input1) # hidden dimension only
        return output1

    def forward_double(self, input1, input2):
        # change second dimension to the last dimension
        output1 = self.model(input1) # hidden dimension only
        output2 = self.model(input2) # hidden dimension only
        # Compute the distance between the anchor and the unknown, both of shape (batch size, embedding size)
        return output1, output2

  
    def training_step(self, batch, batch_idx):
        (x0, x1) , y = batch
        output1, output2 = self.forward_double(x0, x1)
        # print("output1: ", output1, "output2: ", output2, "y: ", y)
        loss, predicted_labels = self.criterion(output1, output2, y.to(torch.float))
        # predicted_labels = predicted_labels.squeeze(1)
        accuracy = torch.mean((predicted_labels == y.to(torch.float)).to(torch.float))
        wandb.log({"epoch": self.current_epoch, "train/loss": loss})
        wandb.log({"epoch": self.current_epoch, "train/accuracy": accuracy})
        return loss

    def validation_step(self, batch, batch_idx):
        (x0, x1) , y = batch
        output1, output2 = self.forward_double(x0, x1)
        print("output1: ", output1, "output2: ", output2)
        loss, predicted_labels = self.criterion(output1, output2, y.to(torch.float))
        accuracy = torch.mean((predicted_labels == y.to(torch.float)).to(torch.float))
        self.log("val_loss", loss,  batch_size = 20) 
        wandb.log({"epoch": self.current_epoch, "val/loss": loss})
        wandb.log({"epoch": self.current_epoch, "val/accuracy": accuracy})
        return {"val_loss": loss}
     
    # def validation_step(self, batch, batch_idx):
    #     x, y = batch
    #     output = self.model(x)
    #     print("output: ", output)
    #     predicted_labels , n_clusters_, homogeneity, completeness  = cluster_embeddings.cluster_embeddings_dbscan(output, y)
    #     print("predicted_labels" , predicted_labels)
    #     wandb.log({"epoch": self.current_epoch, "val/n_clusters": n_clusters_})
    #     wandb.log({"epoch": self.current_epoch, "val/homogeneity": homogeneity})
    #     wandb.log({"epoch": self.current_epoch, "val/completeness": completeness})
    #     return homogeneity

    def test_step(self, batch, batch_idx):
        (x0, x1) , y = batch
        output1, output2 = self.forward_double(x0, x1)
        loss, predicted_labels = self.criterion(output1, output2, y.to(torch.float))
        accuracy = torch.mean((predicted_labels == y.to(torch.float)).to(torch.float))
        self.log('accuracy', accuracy, on_epoch=True)
        f1 = multiclass_f1_score(predicted_labels, y, num_classes=2, average="weighted")
        self.log('f1', f1, on_epoch=True)
        # wandb.log({"epoch": self.current_epoch, "test/loss": loss})
        # wandb.log({"epoch": self.current_epoch, "test/accuracy": accuracy})
        return accuracy

    def predict_step(self, batch, batch_idx):
        x, y = batch
        output = self.model(x)
        # predicted_labels , n_clusters_, homogeneity, completeness  = cluster_embeddings.cluster_embeddings_dbscan(output, y)
        # wandb.log({"epoch": self.current_epoch, "test/n_clusters": n_clusters_})
        # wandb.log({"epoch": self.current_epoch, "test/homogeneity": homogeneity})
        # wandb.log({"epoch": self.current_epoch, "test/completeness": completeness})
        # print("Clusters: ", n_clusters_, "Homogeneity: ", homogeneity, "Completeness: ", completeness)
        return output, y

    def configure_optimizers(self):
        return self.optimizer
    
    def train_dataloader(self):
        return self.train_datal
    
    def test_dataloader(self):
        return self.test_datal
    
    def val_dataloader(self):
        return self.val_datal
    
