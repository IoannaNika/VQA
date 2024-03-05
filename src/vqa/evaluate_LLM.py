import os
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
import torch.optim as optim
from pytorch_lightning.loggers import WandbLogger
import pytorch_lightning as pl
from vqa.loss.ContrastiveSiameseLoss import ContrastiveLoss
from vqa.data.datasets.TripletReads import TripletReads
from vqa.lightning.TransformerTripletNetTrainer import TripletNetTrainer
from pytorch_lightning.callbacks import ModelCheckpoint, StochasticWeightAveraging
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from transformers import TrainingArguments, Trainer, AutoModelForSequenceClassification, AutoModelForMaskedLM
from vqa.data.datasets.AmpliconSiameseReads import SiameseReads as LUMCReads


checkpoint_path = "transformer_checkpoints/final_checkpoint/epoch=3-step=1052.ckpt"
nt = AutoModelForMaskedLM.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-50m-multi-species", num_labels =2,  trust_remote_code=True)
criterion = ContrastiveLoss(0.2)
trainer = pl.Trainer(devices=1, accelerator='gpu', enable_progress_bar=False)
model = TripletNetTrainer.load_from_checkpoint(checkpoint_path, model = nt, train_datal = None, val_datal = None, test_datal = None, optimizer = None, margin = 0.3, criterion = criterion)
model.eval()

genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                    (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                    (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                    (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                    (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                    (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                    (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]

for gr in genomic_regions:
    gr = str(gr[0]) + "_" + str(gr[1])  
    data = LUMCReads(directory= "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data",  genomic_region = gr)

    print("Genomic region: ", gr, "\n")
    if data.length <= 0:
        continue 

    datal = DataLoader(data, batch_size = 20, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)
    print("Total amount of data: ", data.length)
    out =  trainer.test(model, dataloaders = datal)

print("Evaluate all data ...")
data = LUMCReads(directory= "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data")
datal = DataLoader(data, batch_size = 20, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)
out =  trainer.test(model, dataloaders = datal)

