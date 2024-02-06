import os
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
import torch.optim as optim
from pytorch_lightning.loggers import WandbLogger
import pytorch_lightning as pl
import wandb
from vqa.data.datasets.TripletReads import TripletReads
from vqa.lightning.TransformerTripletNetTrainer import TripletNetTrainer
from pytorch_lightning.callbacks import ModelCheckpoint, StochasticWeightAveraging
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from transformers import TrainingArguments, Trainer, AutoModelForSequenceClassification, AutoModelForMaskedLM


def main():
    
    model = AutoModelForMaskedLM.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-50m-multi-species", trust_remote_code=True)

    data = TripletReads(directory='/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/triplet_dataset_primers_template_excl_neg_identicals_RSII_100000')

    train_count = int(len(data)*0.8)
    val_count = int(len(data)*0.1)
    test_count = len(data) - train_count - val_count

    train_data, val_data, test_data = random_split(data, [train_count, val_count, test_count])
    train_datal = DataLoader(train_data, batch_size=1, shuffle=True, pin_memory=True, num_workers=4, prefetch_factor=1)
    val_datal = DataLoader(val_data, batch_size=1 , shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=1)
    test_datal = DataLoader(test_data, batch_size=1, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=1)
    
    optimizer = optim.Adam(model.parameters(), lr=7e-3)

    os.environ["WANDB_DIR"] = "/tmp"
    wandb.init(project="NT_siamese")
    wandb_logger = WandbLogger()
    # log loss per epoch
    wandb_logger.watch(model, log='all') 

    check_point_dir = "transformer_checkpoints/"
    early_stop_callback = EarlyStopping(monitor="val_loss", patience=5, verbose=False, mode="min")
    checkpoint_callback = ModelCheckpoint(dirpath=check_point_dir, save_top_k=3, monitor="val_acc", mode="max")
    margin = 0.3
    triplet_network = TripletNetTrainer(model, train_datal, val_datal, test_datal,optimizer, margin)
    trainer = pl.Trainer(max_epochs=100, logger=wandb_logger, accumulate_grad_batches=100, strategy='ddp_find_unused_parameters_true', callbacks=[checkpoint_callback], devices=1, accelerator='gpu')
    trainer.fit(triplet_network)
    print(checkpoint_callback.best_model_path)
    trainer.save_checkpoint(check_point_dir + "/triplet_network_final.ckpt")
    wandb_logger.experiment.unwatch(model)
    trainer.test()

if __name__ == "__main__":
    main()
