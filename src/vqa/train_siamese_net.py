import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
import torch.optim as optim
from pytorch_lightning.loggers import WandbLogger
import pytorch_lightning as pl
import wandb
from vqa.lightning.SiameseNetTrainer import SiameseNetTrainer
from vqa.data.datasets.SimulatedAmpliconSiameseReads import SiameseReads
# from vqa.data.datasets.SiameseReads import SiameseReads
from vqa.data.datasets.ClusterReads import ClusterReads
from vqa.data.transforms.PadNOneHot import PadNOneHot
from vqa.loss.ContrastiveSiameseLoss import ContrastiveLoss
from vqa.models.LSTM import LSTM
import os
from tcn_lib import TCN
from pytorch_lightning.callbacks import ModelCheckpoint, StochasticWeightAveraging

def main():
    
    max_length = 1525
    # max_length = 29900

    # model = LSTM(4, 64)
    # model = TCN(4, -1, [46]*13, 3, batch_norm = True, weight_norm = True)
    model = TCN(4, -1, [64]*7, 7, batch_norm = True, weight_norm = True)

    # for name, param in model.named_parameters():
    #     # initialize orthogonal weights for recurrent layers
    #     if  "weight_hh" in name:
    #         print("initializing orthogonal weights for recurrent layer")
    #         nn.init.orthogonal_(param)
    #     print(f"Layer: {name}, Mean: {param.data.mean()}, Std: {param.data.std()}")

    transform = PadNOneHot(max_length,"pre")
    # transform_val = PadNOneHot(max_length,"pre", single_read=True)
    data = SiameseReads(directory='/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/dataset_100000_hard', transform=transform)
    # data = SiameseReads(directory='data/data/hcov_global_2023-11-16_09-28/siamese_reference_set_easy', transform=transform)
    # val_data = ClusterReads(directory="data/data/validation_data_2/cluster_reference_set", transform=transform_val)
    train_count = int(len(data)*0.8)
    val_count = int(len(data)*0.1)
    test_count = len(data) - train_count - val_count
    # val_count = val_data.length

    torch.manual_seed(20)
    train_data, val_data, test_data = random_split(data, [train_count, val_count, test_count])
    train_datal = DataLoader(train_data, batch_size=20, shuffle=True, pin_memory=True, num_workers=4, prefetch_factor=8)
    val_datal = DataLoader(val_data, batch_size=20, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)
    test_datal = DataLoader(test_data, batch_size=20, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)
    
    criterion = ContrastiveLoss(2)

    optimizer = optim.Adam(model.parameters(), lr=7e-4)
    scheduler = None #optim.lr_scheduler.ExponentialLR(optimizer, 1e-2)

    os.environ["WANDB_DIR"] = "/tmp"
    wandb.init(project="TCN_Siamese_net")
    wandb_logger = WandbLogger()
    # log loss per epoch
    wandb_logger.watch(model, log='all') 

    # save the model every 10 epochs
    checkpoint_callback = ModelCheckpoint(
        dirpath='checkpoints/',
        filename='siamese_net-{epoch:02d}',
        every_n_epochs=10
    )
    siamese_network = SiameseNetTrainer(model, train_datal, val_datal, test_datal, criterion, optimizer, scheduler)
    # , StochasticWeightAveraging(swa_lrs=1e-2)
    trainer = pl.Trainer(max_epochs=100, logger=wandb_logger, accumulate_grad_batches=50, callbacks=[checkpoint_callback], devices=1, accelerator='gpu')
    trainer.fit(siamese_network)
    trainer.save_checkpoint("checkpoints/siamese_net_final.ckpt")
    wandb_logger.experiment.unwatch(model)
    trainer.test()

if __name__ == "__main__":
    main()
