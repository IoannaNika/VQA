import os
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, random_split
import torch.optim as optim
from pytorch_lightning.loggers import WandbLogger
import pytorch_lightning as pl
from vqa.lightning.TransformerBinaryNetTrainer import TransformerBinaryNetTrainer
from pytorch_lightning.callbacks import ModelCheckpoint, StochasticWeightAveraging
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
from transformers import TrainingArguments, AutoModelForSequenceClassification, AutoModelForMaskedLM
from vqa.data.datasets.AmpliconSiameseReads import SiameseReads as LUMCReads
from peft import LoraConfig, get_peft_model, TaskType, IA3Config
from torchsummary import summary
from peft import PeftConfig, PeftModel
from vqa.data.datasets.SimulatedAmpliconSiameseReads import SiameseReads
import argparse
import sys

def main(): 
    parser = argparse.ArgumentParser(description="Output predictions for pairs of reads")
    parser.add_argument('--virus_name', dest = 'virus_name', required=True, type=str, help="Virus name: SARS-CoV-2, HIV-1, HCV-1b")
    parser.add_argument('--outdir', dest = 'outdir', required=True, type=str, help="output directory")
    parser.add_argument('--path_to_dataset', dest = 'path_to_dataset', required=True, type=str, help="Path to input dataset")
    parser.add_argument('--append_mode', dest ="append_mode", required=False, default="False", type=str, help="to append to existing prediction file or to create a new one and delete an existing one")
    parser.add_argument('--lumc', dest ="lumc", required=False, default="False", type=str, help="lumc data or not")
    parser.add_argument('--devices', dest="devices", required=False, default=1, type=int, help="number of gpus")
    args = parser.parse_args()
    

    out_file_path = os.path.join(args.outdir, "predictions.tsv")

    
    # if it exists, delete it
    if os.path.exists(out_file_path) and args.append_mode == "False":
        os.remove(out_file_path)
    
        # create the file
        f = open(out_file_path, "x")
        f.write("Genomic_region\tSequence_1_id\tSequence_1\tSequence_2_id\tSequence_2\tPredicted_label\tPredicted_probability\tTrue_label\n")
        f.close()
    
    if os.path.exists(out_file_path)  == False: 
        f = open(out_file_path, "x")
        f.write("Genomic_region\tSequence_1_id\tSequence_1\tSequence_2_id\tSequence_2\tPredicted_label\tPredicted_probability\tTrue_label\n")
        f.close()


    nt = AutoModelForSequenceClassification.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", num_labels =2,  trust_remote_code=True, cache_dir = "/tudelft.net/staff-umbrella/ViralQuasispecies/inika/VQA/src/vqa/cache")
    adapter_name = "checkpoint_binary_transformer/7" #"checkpoint_binary_transformer/IA3_gpu_dense_no_query_more_data_500m_168h/11" #"checkpoint_binary_transformer_rev_compl/500m_pac_bio_with_output_dense/8" #"checkpoint_binary_transformer/7" "#checkpoint_binary_transformer/IA3_gpu_dense_no_query_more_data_500m_168h/11"  #"checkpoint_binary_transformer/IA3_gpu_dense_no_query_more_data_500m_long/4" #"checkpoint_binary_transformer/IA3_gpu_dense_no_query_more_data_500m/1" #"checkpoint_binary_transformer/IA3_gpu_dense_no_query/6"

    model = PeftModel.from_pretrained(nt, adapter_name)
    model = model.merge_and_unload()

    model = TransformerBinaryNetTrainer(model = model, train_datal = None, val_datal = None, test_datal = None,optimizer = None, batch_size = 20, checkpoint_dir=None,  device="cuda:0", outdir = out_file_path)

    model.eval()

    trainer = pl.Trainer(devices=args.devices, accelerator='gpu', enable_progress_bar=False)


    genomic_regions = []
    
    if args.virus_name == "SARS-CoV-2":
        genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                            (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                            (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                            (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                            (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                            (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                            (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]
    
    if args.virus_name == "HIV-1":
        genomic_regions = [(30, 1028), (947, 1907), (1823, 2775), (2687, 3621), (3516, 4474), (4374, 5354), (5272, 6237), (6126, 7120), (7029, 8011), (7931, 8862), (8156, 9137)]

    if args.virus_name == "HCV-1b":
        genomic_regions = [(72, 1065), (985, 1946), (1842, 2800), (2703, 3698), (3495, 4459), (4314, 5279), (5215, 6167), (6068, 7008), (6930, 7899), (7740, 8681), (8300, 9280)]


    for gr in genomic_regions: 

        gr = str(gr[0]) + "_" + str(gr[1])  
        if args.lumc == "True":
            data = LUMCReads(directory = args.path_to_dataset, genomic_region = gr, test_mode = True )
        else:
            data = SiameseReads(directory = args.path_to_dataset, test_mode=True, genomic_region = gr)

        print(data.length)
        if data.length == 0: 
            continue
        
        datal = DataLoader(data, batch_size = 20, shuffle=False, pin_memory=True, num_workers=4, prefetch_factor=8)

        out =  trainer.test(model, dataloaders = datal)
        

if __name__ == "__main__":
    sys.exit(main())