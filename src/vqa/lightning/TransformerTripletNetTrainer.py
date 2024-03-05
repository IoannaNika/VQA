import pytorch_lightning as pl
import torch.optim as optim
import torch
import wandb
import vqa.models.cluster_embeddings as cluster_embeddings
import torch.nn.functional as F
import torch.nn as nn
import numpy as np
from transformers import AutoTokenizer
from torcheval.metrics.functional import multiclass_f1_score

class TripletNetTrainer(pl.LightningModule):
    def __init__(self, model, train_datal, val_datal, test_datal,optimizer, margin, criterion=None):
        super().__init__()
        self.model = model
        self.margin = margin
        if criterion == None: 
            criterion=nn.TripletMarginLoss(margin=self.margin, p=2, eps=1e-7)
        self.criterion =criterion
        self.tokenizer = AutoTokenizer.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", trust_remote_code=True, cache_dir = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/cache")
        self.max_length = 400 
        self.optimizer = optimizer
        self.train_datal = train_datal
        self.test_datal = test_datal
        self.val_datal = val_datal

    def forward(self, input1, input2, input3):
        input1 = self.tokenizer.batch_encode_plus(input1, return_tensors="pt", padding="max_length", max_length = self.max_length, truncation=True)["input_ids"].to("cuda")
        input2 = self.tokenizer.batch_encode_plus(input2, return_tensors="pt", padding="max_length", max_length = self.max_length, truncation=True)["input_ids"].to("cuda")
        input3 = self.tokenizer.batch_encode_plus(input3, return_tensors="pt", padding="max_length", max_length = self.max_length, truncation=True)["input_ids"].to("cuda")

        # print("Input 1: ", input1.shape)
        # print("Input 1: ", input1)
        output1 = self.get_embeddings(input1)
        output2 = self.get_embeddings(input2)
        output3 = self.get_embeddings(input3)
        # print("Output 1: ", output1.shape)
        return output1, output2, output3
    
    def forward_double(self, input1, input2): 
        input1 = self.tokenizer.batch_encode_plus(input1, return_tensors="pt", padding="max_length", max_length = self.max_length, truncation=True)["input_ids"].to("cuda")
        input2 = self.tokenizer.batch_encode_plus(input2, return_tensors="pt", padding="max_length", max_length = self.max_length, truncation=True)["input_ids"].to("cuda")
        output1 = self.get_embeddings(input1)
        output2 = self.get_embeddings(input2)
        return output1, output2
  
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
        accuracy = self.get_accuracy(self.margin, output1, output2, output3)
        loss = self.criterion(output1, output2, output3)
        self.log("val_loss", loss,  batch_size = 20) 
        self.log("val_acc", accuracy,  batch_size = 20)
        print("validation accuracy: ", accuracy)
        wandb.log({"epoch": self.current_epoch, "val/loss": loss})
        wandb.log({"epoch": self.current_epoch, "val/accuracy": accuracy})
        return loss
    

    def test_step(self, batch, batch_idx):
        (x0, x1), y = batch
        output1, output2 = self.forward_double(x0, x1)
        loss, predicted_labels = self.criterion(output1, output2, y.to(torch.float))
        accuracy = torch.mean((predicted_labels == y.to(torch.float)).to(torch.float))
        self.log('accuracy', accuracy, on_epoch=True)
        f1 = multiclass_f1_score(predicted_labels, y, num_classes=2, average="weighted")
        self.log('f1', f1, on_epoch=True)
        return loss
      
    def get_accuracy(self, margin, output1, output2, output3):
        euclidean_distance_pos = F.pairwise_distance(output1, output2, keepdim = True) 
        # print("Euclidean pos: ", euclidean_distance_pos)
        euclidean_distance_neg = F.pairwise_distance(output1, output3, keepdim = True) 
        # print("Euclidean neg: ", euclidean_distance_neg)
        predicted_labels_pos = torch.where(euclidean_distance_pos > self.margin, 0, 1)
        # print("Predicted labels pos:", predicted_labels_pos)
        predicted_labels_neg = torch.where(euclidean_distance_neg > self.margin, 0, 1)
        # print("Predicted labels neg:", predicted_labels_neg)
        accuracy_pos = torch.sum(predicted_labels_pos == 1).item() / len(predicted_labels_pos)
        # print("Accuracy pos: ", accuracy_pos)
        accuracy_neg = torch.sum(predicted_labels_neg == 0).item() / len(predicted_labels_neg)
        # print("Accuracy neg:", accuracy_neg)
        accuracy = (accuracy_pos + accuracy_neg) / 2

        return accuracy

    def get_embeddings(self, tokens_ids):
        attention_mask = tokens_ids != self.tokenizer.pad_token_id
        
        output = self.model(tokens_ids,
            attention_mask=attention_mask,
            encoder_attention_mask=attention_mask,
            output_hidden_states=True)
        
        embeddings = output['hidden_states'][-1]
        # print("Embeddings: ", embeddings.shape)
        attention_mask = torch.unsqueeze(attention_mask, dim=-1)
        # print("Attention mask: ", attention_mask, attention_mask.shape)
        mean_sequence_embeddings = torch.sum(attention_mask*embeddings, axis=-2)/torch.sum(attention_mask, axis=1)
        # print("Mean sequence embeddings: ", mean_sequence_embeddings.shape)
        return mean_sequence_embeddings



    
    def configure_optimizers(self):
        return self.optimizer
    
    def train_dataloader(self):
        return self.train_datal
    
    def test_dataloader(self):
        return self.test_datal
    
    def val_dataloader(self):
        return self.val_datal
    
