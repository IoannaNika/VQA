import pytorch_lightning as pl
import torch.optim as optim
import torch
import wandb
import vqa.models.cluster_embeddings as cluster_embeddings
import torch.nn.functional as F
import torch.nn as nn
import numpy as npz
from transformers import AutoTokenizer

class TransformerBinaryNetTrainer(pl.LightningModule):
    def __init__(self, model, train_datal, val_datal, test_datal,optimizer, batch_size):
        super().__init__()
        self.model = model
        # replace head with a binary classification head
        self.model.classifier = nn.Linear(self.model.config.hidden_size, 2)
        self.tokenizer = AutoTokenizer.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-100m-multi-species", trust_remote_code=True)
        # add | as a special token to the tokenizer
        self.tokenizer.add_tokens("|")
        self.batch_size = batch_size
        self.max_length = 500 # self.tokenizer.model_max_length
        self.optimizer = optimizer
        self.train_datal = train_datal
        self.test_datal = test_datal
        self.val_datal = val_datal

    def forward(self, input1, input2):

        # merge input1 and input2 into one tensor and add a special token in between them 
        inpt = input1 + "|" + input2

        inpt = self.tokenizer.batch_encode_plus(inpt, return_tensors="pt", padding="max_length", max_length = self.max_length)["input_ids"].to("cuda:0")

        output = self.get_binary_output(inpt)
        # print("Output 1: ", output1.shape)
        return output

  
    def training_step(self, batch, batch_idx):
        (x0, x1), y = batch
        output = self.forward(x0, x1)
        # anchor, positive, negative
        loss = F.cross_entropy(output, y)
        accuracy = torch.sum(torch.argmax(output, dim=1) == y).item() / len(y)
        wandb.log({"epoch": self.current_epoch, "train/loss": loss})
        wandb.log({"epoch": self.current_epoch, "train/accuracy": accuracy})
        return loss

    def validation_step(self, batch, batch_idx):
        (x0, x1), y = batch
        output = self.forward(x0, x1)
        loss = F.cross_entropy(output, y)
        accuracy = torch.sum(torch.argmax(output, dim=1) == y).item() / len(y)
        self.log("val_loss", loss,  batch_size = 20) 
        self.log("val_acc", accuracy,  batch_size = 20)
        print("validation accuracy: ", accuracy)
        wandb.log({"epoch": self.current_epoch, "val/loss": loss})
        wandb.log({"epoch": self.current_epoch, "val/accuracy": accuracy})
        return loss
    

    def test_step(self, batch, batch_idx):
        (x0, x1), y = batch
        output = self.forward(x0, x1)
        loss = F.cross_entropy(output, y)
        accuracy = torch.sum(torch.argmax(output, dim=1) == y).item() / len(y)
        wandb.log({"epoch": self.current_epoch, "val/loss": loss})
        wandb.log({"epoch": self.current_epoch, "val/accuracy": accuracy})
        return loss

    def get_binary_output(self, inpt):
        # get probabilities of the two classes
        output = self.model(inpt).logits
        return  output
 
    def configure_optimizers(self):
        return self.optimizer
    
    def train_dataloader(self):
        return self.train_datal
    
    def test_dataloader(self):
        return self.test_datal
    
    def val_dataloader(self):
        return self.val_datal
    
