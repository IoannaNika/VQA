import os
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
import torch.optim as optim
from pytorch_lightning.loggers import WandbLogger
import pytorch_lightning as pl
import wandb
from vqa.data.datasets.SimulatedAmpliconSiameseReads import SiameseReads
from vqa.lightning.TransformerBinaryNetTrainer import TransformerBinaryNetTrainer
from pytorch_lightning.callbacks import ModelCheckpoint, StochasticWeightAveraging
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from transformers import TrainingArguments, Trainer, AutoModelForSequenceClassification, AutoModelForMaskedLM, AutoModel
from peft import LoraConfig, get_peft_model, TaskType, IA3Config
from vqa.models.CustomBinaryTransformer import CustomBinaryTransformer
from vqa.data.datasets.AmpliconSiameseReads import SiameseReads as LUMCReads

def main():
    device = "cuda:0"
    n_devices = 2
    check_point_dir = "checkpoint_binary_transformer/IA3_gpu_dense_no_query_more_data_500m_168h"
    max_length = 800
    batch = 10
    model = AutoModelForSequenceClassification.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", trust_remote_code=True, num_labels=2, cache_dir = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/cache")

    peft_config = IA3Config(
    peft_type="IA3", target_modules=["value", "key", "intermediate.dense", "output.dense"], feedforward_modules=["intermediate.dense", "output.dense"], modules_to_save=["classifier"])
    
    peft_model = get_peft_model(model, peft_config)


    for name, param in peft_model.named_parameters():
        print(name)
        if "classifier" in name:
            param.requires_grad = True
            
    data = SiameseReads(directory='/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/triplet_dataset_primers_template_excl_neg_identicals_327896')

    train_count = int(len(data)*0.8)
    val_count = int(len(data)*0.1)
    test_count =  len(data) - train_count - val_count
    
    # val_data = LUMCReads(directory= "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data")

    train_data, val_data, test_data = random_split(data, [train_count, val_count, test_count])
    # train_data, test_data = random_split(data, [train_count, test_count])


    train_datal = DataLoader(train_data, batch_size=10, shuffle=True, pin_memory=True, num_workers=4, prefetch_factor=1)
    val_datal = DataLoader(val_data, batch_size=10 , shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=1)
    test_datal = DataLoader(test_data, batch_size=10, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=1)
    
    optimizer = optim.Adam(peft_model.parameters(), lr=7e-3)

    os.environ["WANDB_DIR"] = "/tmp"
    wandb.init(project="NT_binary")
    wandb_logger = WandbLogger()
    # log loss per epoch
    wandb_logger.watch(peft_model) 
    
    early_stop_callback = EarlyStopping(monitor="val_acc", patience=3, verbose=False, mode="max")
    # checkpoint_callback = ModelCheckpoint(dirpath=check_point_dir, save_top_k=3, monitor="val_acc", mode="max")
    binary_transformer = TransformerBinaryNetTrainer(peft_model, train_datal, val_datal, test_datal,optimizer, batch, device = device, checkpoint_dir = check_point_dir, max_length=max_length)
    trainer = pl.Trainer(max_epochs=100, logger=wandb_logger,  accumulate_grad_batches=100, strategy='ddp_find_unused_parameters_true', callbacks=[early_stop_callback], devices=n_devices, accelerator="gpu", enable_progress_bar=False)
    trainer.fit(binary_transformer)
    # trainer.save_checkpoint(check_point_dir + "/triplet_network_final.ckpt")
    wandb_logger.experiment.unwatch(peft_model)
    trainer.test()

if __name__ == "__main__":
    main()
