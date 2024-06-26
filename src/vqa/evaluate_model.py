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
from tcn_lib import TCN
from vqa.loss.ContrastiveSiameseLoss import ContrastiveLoss
from pytorch_lightning.callbacks import ModelCheckpoint, StochasticWeightAveraging
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from vqa.data.datasets.AmpliconClusterReads import ClusterReads as AmpliconClusterReads
from vqa.data.datasets.AmpliconSiameseReads import SiameseReads as LUMCReads
import vqa.models.cluster_embeddings as cluster_embeddings

max_length = 1525
lstm = LSTM(4,32)
tcn = TCN(4, -1, [256]*7, 7, batch_norm = True, weight_norm = True)


# os.environ["WANDB_DIR"] = "/tmp"
# wandb.init(project="TCN_Siamese_net")
# wandb_logger = WandbLogger()

# wandb_logger.watch(lstm, log='all') 

genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                    (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                    (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                    (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                    (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                    (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                    (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]

# genomic_regions = [(3166, 4240)]
checkpoint_path =  "checkpoints/epoch=19-step=1600.ckpt"
checkpoint = torch.load(checkpoint_path)

criterion = ContrastiveLoss(0.3)


model = SiameseNetTrainer.load_from_checkpoint(checkpoint_path, model = tcn, train_datal = None, val_datal = None, test_datal = None, criterion = criterion, optimizer = None, scheduler = None)
model.eval()
trainer = pl.Trainer(devices=1, accelerator='gpu', enable_progress_bar=False)
outputs = dict()
labels = dict()

for gr in genomic_regions:
    gr = str(gr[0]) + "_" + str(gr[1])  
    outputs[gr] = []
    labels[gr] = []
    # data = LUMCClusterReads(directory = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data/", transform=PadNOneHot(max_length,"pre", single_read=True), genomic_region = gr)
    # data = AmpliconClusterReads(directory="/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/test_datasets/sars_cov2_BA.1.1/dataset/", transform=PadNOneHot(max_length,"pre", single_read=True),  genomic_region = gr)
    data = LUMCReads(directory= "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data", transform= PadNOneHot(max_length,"pre"), genomic_region=gr)
    # data = SiameseReads(directory = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/test_datasets/test_set_BA.1.1/dataset", transform=PadNOneHot(max_length,"pre", single_read=False), genomic_region = gr)
    print("Genomic region: ", gr, "\n")
    if data.length <= 0:
        continue 
        

    datal = DataLoader(data, batch_size = 20, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)
    print("Total amount of data: ", data.length)
    out =  trainer.test(model, dataloaders = datal)
    # print("Accuracy: ", torch.mean(torch.FloatTensor(out)))
    # for batch in out:
    #     for item in batch[0]:
    #         outputs[gr].append(item)
    #     for item in batch[1]:
    #         labels[gr].append(item)
    # predicted_labels , n_clusters_, homogeneity, completeness  = cluster_embeddings.cluster_embeddings_dbscan(outputs[gr], labels[gr], genomic_region=gr, produce_plots=True, eps=0.3, verbose=False)
    # print("Clusters: ", n_clusters_, " Homogeneity: ", homogeneity, " Completeness: ", completeness)
    

print("Evaluate all data ...")
data = LUMCReads(directory= "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/Read_simulators/data/lumc_data", transform=PadNOneHot(max_length,"pre", single_read=False))
datal = DataLoader(data, batch_size = 20, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)
out =  trainer.test(model, dataloaders = datal)

# model.eval()
# # predict with the model
# y_hat = net.test(datal)
# print(y_hat)