import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
import torch.optim as optim
from pytorch_lightning.loggers import WandbLogger
import pytorch_lightning as pl
import wandb
from vqa.lightning.SiameseNetTrainer import SiameseNetTrainer
from vqa.data.datasets.SiameseReads import SiameseReads
from vqa.data.datasets.ClusterReads import ClusterReads
from vqa.data.transforms.PadNOneHot import PadNOneHot
from vqa.loss.ContrastiveSiameseLoss import ContrastiveLoss
from vqa.loss.ContrastiveCrossEntropy import ContrastiveLoss as CBCE_loss
import os
from tcn_lib import TCN
from pytorch_lightning.callbacks import ModelCheckpoint

max_length = 29900

model = TCN(4, -1, [96]*13, 3, batch_norm = True, weight_norm = True)

transform = PadNOneHot(max_length,"pre")
transform_val = PadNOneHot(max_length,"pre", single_read=True)
data = SiameseReads(directory='data/data/hcov_global_2023-11-16_09-28/siamese_reference_set_20k_1k', transform=transform)
val_data = ClusterReads(directory="data/data/validation_data_2/cluster_reference_set", transform=transform_val)
train_count = int(len(data)*0.8)
test_count = len(data) - train_count 
val_count = val_data.length

torch.manual_seed(0)
train_data, test_data = random_split(data, [train_count, test_count])
train_datal = DataLoader(train_data, batch_size=8, shuffle=True, pin_memory=True, num_workers=8,prefetch_factor=8)
val_datal = DataLoader(val_data, batch_size=val_count, shuffle=False, pin_memory=True, num_workers=8,prefetch_factor=8)
test_datal = DataLoader(test_data, batch_size=8, shuffle=False, pin_memory=True, num_workers=8,prefetch_factor=8)

criterion = ContrastiveLoss(margin=0.2)

optimizer = optim.Adam(model.parameters(), lr=1e-3)
scheduler = None

os.environ["WANDB_DIR"] = "/tmp"
wandb.init(project="TCN_siamese_net")
wandb_logger = WandbLogger()
# log loss per epoch
wandb_logger.watch(model, log='all') 

# save the model every 10 epochs
checkpoint_callback = ModelCheckpoint(
    dirpath='checkpoints_bceloss/',
    filename='siamese_net-{epoch:02d}',
    every_n_epochs=10
)
siamese_network = SiameseNetTrainer(model, train_datal, val_datal, test_datal, criterion, optimizer, scheduler)
trainer = pl.Trainer(max_epochs = 60, logger=wandb_logger, callbacks=[checkpoint_callback], devices=5, accelerator="gpu")
trainer.fit(siamese_network)
trainer.save_checkpoint("siamese_net_final_bceloss.ckpt")
wandb_logger.experiment.unwatch(model)
trainer.test()
