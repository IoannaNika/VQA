import pytorch_lightning as pl
import torch.optim as optim
import torch
import wandb
import vqa.models.cluster_embeddings as cluster_embeddings


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
        # output1 = self.model(input1)
        # output2 = self.model(input2)
        concat_two_reads = torch.cat((input1, input2), dim=0)
        output = self.model(concat_two_reads)
        output1, output2 = torch.split(output,input1.shape[0])
        return output1, output2

    def training_step(self, batch, batch_idx):
        (x0, x1) , y = batch
        output1,output2 = self.forward(x0, x1)
        loss = self.criterion(output1,output2, y)
        self.log('train_loss', loss, prog_bar=True)
        wandb.log({"train/loss": loss})
        return loss

    def validation_step(self, batch, batch_idx):
        x , y = batch
        output = self.model(x)
        predicted_labels , n_clusters_, homogeneity, completeness  = cluster_embeddings.cluster_embeddings(output, y)
        wandb.log({"val/n_clusters": n_clusters_})
        wandb.log({"val/homogeneity": homogeneity})
        wandb.log({"val/completeness": completeness})
        return homogeneity

    def test_step(self, batch, batch_idx):
        (x0, x1) , y = batch
        output1,output2 = self.forward(x0, x1)
        loss = self.criterion(output1,output2, y)
        self.log('test_loss', loss, prog_bar=True)
        wandb.log({"test/loss": loss})
        return loss
    
    def configure_optimizers(self):
        optimizer = optim.Adam(self.model.parameters(), lr=2e-2)
        return optimizer
    
    def train_dataloader(self):
        return self.train_datal
    
    def test_dataloader(self):
        return self.test_datal
    
    def val_dataloader(self):
        return self.val_datal
    
