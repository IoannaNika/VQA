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
from peft import LoraConfig, get_peft_model, TaskType, IA3Config
from vqa.data.datasets.AmpliconSiameseReads import SiameseReads as LUMCReads


def main():
    model = AutoModelForMaskedLM.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", trust_remote_code=True, cache_dir = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/cache")
    
    data = TripletReads(directory='/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/triplet_dataset_primers_template_excl_neg_identicals_327896')

    train_count = int(len(data)*0.8)
    val_count = int(len(data)*0.1)
    test_count = len(data) - train_count - val_count

    # peft_config = IA3Config(
    # peft_type="IA3", target_modules=["value", "key", "intermediate.dense", "output.dense"], feedforward_modules=["intermediate.dense", "output.dense"])
    
    # model = get_peft_model(model, peft_config)

    for name, param in model.named_parameters():
        param.requires_grad = True


    train_data, val_data, test_data = random_split(data, [train_count, val_count, test_count])
    # train_data, test_data = random_split(data, [train_count, test_count])
    # val_data = LUMCReads(directory= "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data")


    train_datal = DataLoader(train_data, batch_size=3, shuffle=True, pin_memory=True, num_workers=4, prefetch_factor=1)
    val_datal = DataLoader(val_data, batch_size=3 , shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=1)
    test_datal = DataLoader(test_data, batch_size=3, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=1)
    
    adam = optim.Adam(model.parameters(), lr=7e-3)
    optimizer = {"optimizer": adam, "lr_scheduler": optim.lr_scheduler.StepLR(adam, step_size=1, gamma=0.1)}

    os.environ["WANDB_CONFIG_DIR"]= "/tmp"
    os.environ["WANDB_DIR"] ="/tmp"
    os.environ["WANDB_CACHE_DIR"] = "/tmp"

    wandb.init(project="NT_siamese", entity="inika")
    wandb_logger = WandbLogger()
    # log loss per epoch
    wandb_logger.watch(model) 

    check_point_dir = "checkpoint_triplet_transformer/full_fine_tuning_more_data_500"
    early_stop_callback = EarlyStopping(monitor="val_acc", patience=3, verbose=False, mode="max")
    checkpoint_callback = ModelCheckpoint(dirpath=check_point_dir, save_top_k=3, monitor="val_acc", mode="max")
    margin = 0.2
    triplet_network = TripletNetTrainer(model, train_datal, val_datal, test_datal,optimizer, margin)
    trainer = pl.Trainer(max_epochs=100, logger = wandb_logger, accumulate_grad_batches=333, strategy='ddp_find_unused_parameters_true', callbacks=[checkpoint_callback, early_stop_callback], devices=3, accelerator='gpu', deterministic = True,  enable_progress_bar=False)
    trainer.fit(triplet_network)
    trainer.save_checkpoint(check_point_dir + "/triplet_network_final.ckpt")
    wandb_logger.experiment.unwatch(model)
    trainer.test()

if __name__ == "__main__":
    main()
