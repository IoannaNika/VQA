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
from vqa.data.datasets.AmpliconSiameseReads import SiameseReads as LUMCReads

def main():
    device = "cuda:0"
    n_devices = 2
    check_point_dir = "checkpoint_binary_transformer_rev_compl/random_init_pb"
    max_length = 800 
    batch = 10
    model = AutoModelForSequenceClassification.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", trust_remote_code=True, num_labels=2, cache_dir = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/cache")
    
    # re-initialize all the weights of the model to random values
    # for name, param in model.named_parameters():
    #     nn.init.xavier_uniform_(param)

    # , "output.dense"  
    peft_config = IA3Config(
    peft_type="IA3", target_modules=["value", "key", "intermediate.dense"], feedforward_modules=["intermediate.dense"], modules_to_save=["classifier"])
    
    peft_model = get_peft_model(model, peft_config)

    for name, param in peft_model.named_parameters():
        print(name)
        if "classifier" in name:
            param.requires_grad = True
            
    data = SiameseReads(directory='/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/tuples_pacbio_sars_cov_2_rev_compl_more/dataset')

    train_count = int(len(data)*0.8)
    val_count = int(len(data)*0.1)
    test_count =  len(data) - train_count - val_count
    
    # val_data = LUMCReads(directory= "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data")

    train_data, val_data, test_data = random_split(data, [train_count, val_count, test_count])
    # train_data, test_data = random_split(data, [train_count, test_count])


    train_datal = DataLoader(train_data, batch_size=batch, shuffle=True, pin_memory=True, num_workers=4, prefetch_factor=1)
    val_datal = DataLoader(val_data, batch_size=batch , shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=1)
    test_datal = DataLoader(test_data, batch_size=batch, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=1)
      
    optimizer = optim.Adam(peft_model.parameters(), lr=7e-3)
    # optimizer = {"optimizer": adam, "lr_scheduler": optim.lr_scheduler.StepLR(adam, step_size=1, gamma=0.1)}

    os.environ["WANDB_DIR"] = "/tmp"
    os.environ["WANDB_START_METHOD"]="thread"
    wandb.init(project="NT_binary")
    wandb_logger = WandbLogger()
    # log loss per epoch
    wandb_logger.watch(peft_model) 
    
    early_stop_callback = EarlyStopping(monitor="val_acc", patience=3, verbose=False, mode="max")
    # checkpoint_callback = ModelCheckpoint(dirpath=check_point_dir, save_top_k=3, monitor="val_acc", mode="max")
    binary_transformer = TransformerBinaryNetTrainer(peft_model, train_datal, val_datal, test_datal,optimizer, batch, device = device, checkpoint_dir = check_point_dir, max_length=max_length, outdir = None)
    trainer = pl.Trainer(max_epochs=15, logger=wandb_logger,  accumulate_grad_batches=100, strategy='ddp_find_unused_parameters_true', callbacks=[early_stop_callback], devices=n_devices, accelerator="gpu", enable_progress_bar=False)
    trainer.fit(binary_transformer)
    wandb_logger.experiment.unwatch(peft_model)
    trainer.test()

if __name__ == "__main__":
    main()
