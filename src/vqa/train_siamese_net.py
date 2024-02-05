import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
import torch.optim as optim
from pytorch_lightning.loggers import WandbLogger
import pytorch_lightning as pl
import wandb
from vqa.lightning.SiameseNetTrainer import SiameseNetTrainer
from vqa.data.datasets.SimulatedAmpliconSiameseReads import SiameseReads
from vqa.data.datasets.AmpliconSiameseReads import SiameseReads as LUMCReads
# from vqa.data.datasets.SiameseReads import SiameseReads
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from vqa.data.datasets.ClusterReads import ClusterReads
from vqa.data.datasets.LumcClusterReads import ClusterReads as LUMCClusterReads
from vqa.data.datasets.AmpliconClusterReads import ClusterReads as AmpliconClusterReads
from vqa.data.transforms.PadNOneHot import PadNOneHot
from vqa.loss.ContrastiveSiameseLoss import ContrastiveLoss
from vqa.models.LSTM import LSTM
import os
from tcn_lib import TCN
from pytorch_lightning.callbacks import ModelCheckpoint, StochasticWeightAveraging

def main():
    
    max_length = 1525
    # max_length = 29900

    model = LSTM(4, 64)
    # model = TCN(4, -1, [32]*7, 7, batch_norm = True, weight_norm = True)
    # model = TCN(4, -1, [56]*7, 7, batch_norm = True, weight_norm = True)

    transform = PadNOneHot(max_length,"pre")
    transform_val = PadNOneHot(max_length,"pre", single_read=True)
    data = SiameseReads(directory='/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/triplet_dataset_test_100000', transform=transform)
    # val_data = SiameseReads(directory="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/2022_val_dataset", transform=transform)
    # val_data = LUMCReads(directory= "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data", transform=transform)
    # val_data = AmpliconClusterReads(directory="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/2022_val_dataset", transform=transform_val)
    # val_data = LUMCClusterReads(directory = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data",transform=transform_val)
    train_count = int(len(data)*0.8)
    val_count = int(len(data)*0.1)
    test_count = len(data) - train_count - val_count
    # val_count = val_data.length

    torch.manual_seed(20)
    train_data, val_data, test_data = random_split(data, [train_count, val_count, test_count])
    # train_data, test_data = random_split(data, [train_count, test_count])
    train_datal = DataLoader(train_data, batch_size=20, shuffle=True, pin_memory=True, num_workers=4, prefetch_factor=8)
    val_datal = DataLoader(val_data, batch_size=20 , shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)
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
        dirpath='checkpoints_triplets_primers/',
        filename='siamese_net-{epoch:02d}',
        every_n_epochs=1
    )
    early_stop_callback = EarlyStopping(monitor="val_loss", patience=5, verbose=False, mode="min")

    siamese_network = SiameseNetTrainer(model, train_datal, val_datal, test_datal, criterion, optimizer, scheduler)
    trainer = pl.Trainer(max_epochs=100, logger=wandb_logger, accumulate_grad_batches=50, callbacks=[checkpoint_callback, early_stop_callback], devices=1, accelerator='gpu')
    trainer.fit(siamese_network)
    print(checkpoint_callback.best_model_path)
    trainer.save_checkpoint("checkpoints_triplets_primers/siamese_net_final.ckpt")
    wandb_logger.experiment.unwatch(model)
    trainer.test()

if __name__ == "__main__":
    main()
