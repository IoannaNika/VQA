import torch
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

from tcn_lib import TCN

max_length = 29900

model = TCN(4, -1, [64]*13, 3)

transform = PadNOneHot(max_length,"pre")
transform_val = PadNOneHot(max_length,"pre", single_read=True)
data = SiameseReads(directory='data/data/hcov_global_2023-11-16_09-28/siamese_reference_set', transform=transform)
val_data = ClusterReads(directory="data/data/hcov_global_2022-08-26_05-12/validation_data_1/cluster_reference_set", transform=transform_val)
train_count = int(len(data)*0.8)
test_count = len(data) - train_count 
val_count = val_data.length
print(train_count, val_count, test_count)
torch.manual_seed(0)
train_data, test_data = random_split(data, [train_count, test_count])
train_datal = DataLoader(train_data, batch_size=5, shuffle=True)
val_datal = DataLoader(val_data, batch_size=val_count, shuffle=False)
test_datal = DataLoader(test_data, batch_size=5, shuffle=False)

criterion = ContrastiveLoss(margin=1.0)
optimizer = optim.Adam(model.parameters(), lr=2e-2)
scheduler = None

wandb.init(project="TCN_siamese_net")
wandb_logger = WandbLogger()
# log loss per epoch
wandb_logger.watch(model, log='all', log_freq=80) 
siamese_network = SiameseNetTrainer(model, train_datal, val_datal, test_datal, criterion, optimizer, scheduler)
trainer = pl.Trainer(max_epochs = 5, logger=wandb_logger)
trainer.fit(siamese_network)
trainer.save_checkpoint("siamese_net.ckpt")
wandb_logger.experiment.unwatch(model)
trainer.test()


# print the model weights
# for name, param in model.named_parameters():
#     if param.requires_grad:
#         print(name, param.data)