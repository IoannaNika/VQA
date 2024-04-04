import os
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
import torch.optim as optim
from pytorch_lightning.loggers import WandbLogger
import pytorch_lightning as pl
from vqa.loss.ContrastiveSiameseLoss import ContrastiveLoss
from vqa.data.datasets.TripletReads import TripletReads
from vqa.lightning.TransformerBinaryNetTrainer import TransformerBinaryNetTrainer
from pytorch_lightning.callbacks import ModelCheckpoint, StochasticWeightAveraging
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from transformers import TrainingArguments, Trainer, AutoModelForSequenceClassification, AutoModelForMaskedLM
from vqa.data.datasets.AmpliconSiameseReads import SiameseReads as LUMCReads
from peft import LoraConfig, get_peft_model, TaskType, IA3Config
from torchsummary import summary
from peft import PeftConfig, PeftModel
from vqa.data.datasets.SimulatedAmpliconSiameseReads import SiameseReads


nt = AutoModelForSequenceClassification.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", num_labels =2,  trust_remote_code=True, cache_dir = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/cache")
adapter_name = "checkpoint_binary_transformer/IA3_gpu_dense_no_query_more_data_500m_168h/11" #"checkpoint_binary_transformer/IA3_gpu_dense_no_query_more_data_500m_long/4" #"checkpoint_binary_transformer/IA3_gpu_dense_no_query_more_data_500m/1" #"checkpoint_binary_transformer/IA3_gpu_dense_no_query/6"

# model = PeftModel.from_pretrained(nt, adapter_name)
# model = model.merge_and_unload()

model = TransformerBinaryNetTrainer(model = nt, train_datal = None, val_datal = None, test_datal = None,optimizer = None, batch_size = 20, checkpoint_dir=None,  device="cuda:0")

model.eval()

trainer = pl.Trainer(devices=1, accelerator='gpu', enable_progress_bar=False)

# sars-cov-2
# genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
#                     (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
#                     (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
#                     (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
#                     (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
#                     (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
#                     (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]

# HIV-1
# genomic_regions =[(140, 1081), (980, 1927), (1824, 2807), (2696, 3679), (3583, 4564), (4429, 5398), (5291, 6249), (6143, 7140), (6968, 7955), (7864, 8844), (8053, 8970)]

# HCV-1b
genomic_regions = [(72, 1065), (985, 1946), (1842, 2800), (2703, 3698), (3495, 4459), (4314, 5279), (5215, 6167), (6068, 7008), (6930, 7899), (7740, 8681), (8300, 9280)]

# for gr in genomic_regions:
#     gr = str(gr[0]) + "_" + str(gr[1])  
#     data = SiameseReads(directory = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/HCV-1b-NCBI/subselection_95/dataset/test_pairs.tsv", genomic_region = gr, test_mode=True)
#     # data = LUMCReads(directory= "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data",  genomic_region = gr)
#     datal = DataLoader(data, batch_size = 20, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)
#     print("Genomic region: ", gr, "\n")
#     if data.length <= 0:
#         continue 
        
#     print("Total amount of data: ", data.length)
#     out =  trainer.test(model, dataloaders = datal)

print("Evaluate all data ...")
# data = LUMCReads(directory= "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data/lumc/dataset", test_mode=True)
data = SiameseReads(directory = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/Experiments_data/HCV-1b/95/cov_100x_mixture_pb_hifi/dataset", test_mode=True)
datal = DataLoader(data, batch_size = 20, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)
out =  trainer.test(model, dataloaders = datal)