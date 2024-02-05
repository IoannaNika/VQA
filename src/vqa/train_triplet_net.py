import os
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
import torch.optim as optim
from pytorch_lightning.loggers import WandbLogger
import pytorch_lightning as pl
import wandb
from vqa.data.datasets.TripletReads import TripletReads
from vqa.lightning.TripletNetTrainer import TripletNetTrainer
from vqa.data.transforms.PadNOneHot import PadNOneHot
from vqa.models.LSTM import LSTM
from tcn_lib import TCN
from pytorch_lightning.callbacks import ModelCheckpoint, StochasticWeightAveraging
from pytorch_lightning.callbacks.early_stopping import EarlyStopping

def main():
    
    max_length = 2045
    # model = LSTM(4, 256)
    model = TCN(4, -1, [256]*3, 9, batch_norm = True, weight_norm = True)

    transform = PadNOneHot(max_length,"pre", single_read=True)
    data = TripletReads(directory='/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/triplet_dataset_primers_template_excl_neg_identicals_RSII_100000', transform=transform)
    train_count = int(len(data)*0.8)
    val_count = int(len(data)*0.1)
    test_count = len(data) - train_count - val_count

    train_data, val_data, test_data = random_split(data, [train_count, val_count, test_count])
    train_datal = DataLoader(train_data, batch_size=20, shuffle=True, pin_memory=True, num_workers=4, prefetch_factor=8)
    val_datal = DataLoader(val_data, batch_size=20 , shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)
    test_datal = DataLoader(test_data, batch_size=20, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)
    
    optimizer = optim.Adam(model.parameters(), lr=7e-3)

    os.environ["WANDB_DIR"] = "/tmp"
    wandb.init(project="TCN_Siamese_net")
    wandb_logger = WandbLogger()
    # log loss per epoch
    wandb_logger.watch(model, log='all') 

    # save the model every 10 epochs
    # checkpoint_callbacks = ModelCheckpoint(
    #     dirpath='checkpoints_triplet_primers_v4/',
    #     filename='siamese_net-{epoch:02d}',
    #     every_n_epochs=1
    # )
    
    check_point_dir = "checkpoints_triplet_primers_fixed_max_val_e0.3m0.3/best_rsii_3x9/"
    early_stop_callback = EarlyStopping(monitor="val_loss", patience=5, verbose=False, mode="min")
    checkpoint_callback = ModelCheckpoint(dirpath=check_point_dir, save_top_k=3, monitor="val_acc", mode="max")

    triplet_network = TripletNetTrainer(model, train_datal, val_datal, test_datal,optimizer)
    trainer = pl.Trainer(max_epochs=100, logger=wandb_logger, accumulate_grad_batches=50, callbacks=[checkpoint_callback], devices=3, accelerator='gpu')
    trainer.fit(triplet_network)
    print(checkpoint_callback.best_model_path)
    trainer.save_checkpoint(check_point_dir + "/triplet_network_final.ckpt")
    wandb_logger.experiment.unwatch(model)
    trainer.test()

if __name__ == "__main__":
    main()
