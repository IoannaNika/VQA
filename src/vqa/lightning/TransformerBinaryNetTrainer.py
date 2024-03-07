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
    def __init__(self, model, train_datal, val_datal, test_datal, optimizer, batch_size, checkpoint_dir, max_length=800, device="cpu"):
        super().__init__()
        self.model = model
        self.tokenizer = AutoTokenizer.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", trust_remote_code=True)
        self.batch_size = batch_size
        self.max_length = max_length # self.tokenizer.model_max_length
        self.optimizer = optimizer
        self.train_datal = train_datal
        self.test_datal = test_datal
        self.val_datal = val_datal
        self.criterion = nn.CrossEntropyLoss()
        self.device = device
        self.checkpoint_dir = checkpoint_dir
        self.predictions_file_name = "predictions.tsv"

        # open file to write the predictions
        self.file = open(self.predictions_file_name, 'w')
        self.file.write("Genomic_region\tSequence_1_id\tSequence_1\tSequence_2_id\tSequence_2\tPredicted_label\tPredicted_probability\tTrue_label\n")


    def forward(self, input1, input2):

        inpt = []

        for i in range(len(input1)):
            inpt.append(input1[i].strip() + "<mask><mask><mask><mask><mask><mask><mask><mask><mask><mask><mask><mask><mask><mask>" + input2[i].strip())

        inpt = self.tokenizer.batch_encode_plus(inpt, return_tensors="pt", padding="max_length", max_length = self.max_length, truncation=True)["input_ids"].to(self.device)
        
        attention_mask = inpt != self.tokenizer.pad_token_id
        
        output = self.model(inpt,
            attention_mask=attention_mask,
            # encoder_attention_mask=attention_mask,
            output_hidden_states=True)

        return  output.logits
   

  
    def training_step(self, batch, batch_idx):
        (x0, x1), y = batch
        output = self.forward(x0, x1)
        loss = self.criterion(output, y)
        accuracy = torch.sum(torch.argmax(output, dim=1) == y).item() / len(y)
        wandb.log({"epoch": self.current_epoch, "train/loss": loss})
        wandb.log({"epoch": self.current_epoch, "train/accuracy": accuracy})
        return loss

    def validation_step(self, batch, batch_idx):
        self.model.save_pretrained("{}/{}".format(self.checkpoint_dir, self.current_epoch))
        (x0, x1), y = batch
        # y = torch.Tensor(y)
        # print("Y: ", y)
        output = self.forward(x0, x1)
        loss = self.criterion(output, y)
        accuracy = torch.sum(torch.argmax(output, dim=1) == y).item() / len(y)
        self.log("val_loss", loss,  batch_size = self.batch_size) 
        self.log("val_acc", accuracy,  batch_size = self.batch_size)
        print("validation accuracy: ", accuracy)
        wandb.log({"epoch": self.current_epoch, "val/loss": loss})
        wandb.log({"epoch": self.current_epoch, "val/accuracy": accuracy})
        return loss
    

    def test_step(self, batch, batch_idx):
        (read_id_1,read_id_2), (x0, x1), y, gr = batch
        output = self.forward(x0, x1)
        loss = self.criterion(output, y)
        accuracy = torch.sum(torch.argmax(output, dim=1) == y).item() / len(y)
        self.log("val_acc", accuracy,  batch_size = self.batch_size)
        # write in file the pair the prediction and the true label
        for i in range(len(y)):
            print("Y: ", y)
            self.file.write(gr[i] + "\t" + read_id_1[i] + "\t" + x0[i] + "\t" + read_id_2[i] + "\t" + x1[i] + "\t" + str(torch.argmax(output[i]).item()) + "\t" + str(torch.max(F.softmax(output[i]), dim=0).values.item()) + "\t" + str(y[i].item()) + "\n")  
        return loss
 
    def configure_optimizers(self):
        return self.optimizer
    
    def train_dataloader(self):
        return self.train_datal
    
    def test_dataloader(self):
        return self.test_datal
    
    def val_dataloader(self):
        return self.val_datal

    def device(self): 
        return self.device
    
