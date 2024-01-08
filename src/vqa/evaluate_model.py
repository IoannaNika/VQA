from vqa.data.datasets.LumcClusterReads import ClusterReads as LUMCClusterReads
from vqa.data.transforms.PadNOneHot import PadNOneHot
from torch.utils.data import DataLoader, random_split
from vqa.models.LSTM import LSTM
from vqa.lightning.SiameseNetTrainer import SiameseNetTrainer
import torch
from vqa.loss.ContrastiveSiameseLoss import ContrastiveLoss
from vqa.data.datasets.SimulatedAmpliconSiameseReads import SiameseReads
import torch.optim as optim
from collections import OrderedDict
import wandb
from pytorch_lightning.loggers import WandbLogger
import os
import pytorch_lightning as pl
from pytorch_lightning.callbacks import ModelCheckpoint, StochasticWeightAveraging
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from vqa.data.datasets.AmpliconClusterReads import ClusterReads as AmpliconClusterReads

max_length = 1525
lstm = LSTM(4,32)

os.environ["WANDB_DIR"] = "/tmp"
wandb.init(project="TCN_Siamese_net")
wandb_logger = WandbLogger()

wandb_logger.watch(lstm, log='all') 


transform = PadNOneHot(max_length,"pre", single_read=True)
genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                    (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                    (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                    (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                    (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                    (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                    (26766, 27872), (27808, 28985), (28699, 29768)]
checkpoint_path = "checkpoints_p7/siamese_net-epoch=99.ckpt"
checkpoint = torch.load(checkpoint_path)

    
model = SiameseNetTrainer.load_from_checkpoint(checkpoint_path, model = lstm, train_datal = None, val_datal = None, test_datal = None, criterion = None, optimizer = None, scheduler = None)
model.eval()
trainer =  trainer = pl.Trainer(logger=wandb_logger, devices=1, accelerator='gpu')

for gr in genomic_regions:
    gr = str(gr[0]) + "_" + str(gr[1])   
    # data = LUMCClusterReads(directory = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data", transform=transform, genomic_region = gr)
    data = AmpliconClusterReads(directory="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/2022_val_dataset", transform=transform,  genomic_region = gr)

    if data.length <= 0:
        continue
    datal = DataLoader(data, batch_size=data.length, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)
    trainer.test(model, dataloaders = datal)



# model.eval()
# # predict with the model
# y_hat = net.test(datal)
# print(y_hat)