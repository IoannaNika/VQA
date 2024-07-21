import argparse
import pickle
import sys
import pandas as pd
from sklearn.metrics import accuracy_score, f1_score
import os

def main(): 
    parser = argparse.ArgumentParser(description="Output predictions for pairs of reads")
    parser.add_argument('--virus_name', dest = 'virus_name', required=True, type=str, help="Virus name: SARS-CoV-2, HIV-1, HCV-1b")
    parser.add_argument('--model', dest = 'model', required=True, type=str, help="path to model")
    parser.add_argument('--outdir', dest = 'outdir', required=True, type=str, help="path to model")
    parser.add_argument('--dataset', dest = 'dataset', required=True, type=str, help="TSV dataset")
    parser.add_argument('--distance', dest = 'distance', required=True, type=str, help='edit_distance,n_mutations')

    args = parser.parse_args()

    if args.virus_name == "SARS-CoV-2":
        genomic_regions = [(54, 1183), (1128, 2244), (2179, 3235), (3166, 4240), (4189, 5337),
                            (5286, 6358), (6307, 7379), (7328, 8363), (8282, 9378), (9327, 10429),
                            (10370, 11447), (11394, 12538), (12473, 13599), (13532, 14619),
                            (14568, 15713), (15634, 16698), (16647, 17732), (17649, 18684),
                            (18618, 19655), (19604, 20676), (20581, 21620), (21562, 22590),
                            (22537, 23609), (23544, 24714), (24658, 25768), (25712, 26835),
                            (26766, 27872), (27808, 28985), (28699, 29768), (29768, 29790)]
    
    if args.virus_name == "HIV-1":
        genomic_regions =[(140, 1081), (980, 1927), (1824, 2807), (2696, 3679), (3583, 4564), (4429, 5398), (5291, 6249), (6143, 7140), (6968, 7955), (7864, 8844), (8053, 8970)]

    if args.virus_name == "HCV-1b":
        genomic_regions = [(72, 1065), (985, 1946), (1842, 2800), (2703, 3698), (3495, 4459), (4314, 5279), (5215, 6167), (6068, 7008), (6930, 7899), (7740, 8681), (8300, 9280)]


    # # load model
    filename = args.model

    model = pickle.load(open(filename, 'rb'))

    dataset  = pd.read_csv(args.dataset, sep="\t", header=0)
 
    outfile = os.path.join(args.outdir, "predictions.tsv")

    with open(outfile, "w") as f:
        f.write("Genomic_region\tAccuracy\tf1_score\n")
    for gr in genomic_regions:
        start = gr[0]
        end = gr[1]
        if args.distance == "edit_distance":
            samples = dataset[(dataset["start"] == start) & (dataset["end"] == end)]["edit_distance"].tolist()
        if args.distance == "n_mutations":
            samples = dataset[(dataset["start"] == start) & (dataset["end"] == end)]["n_mutations"].tolist()

        samples = [[s] for s in samples]
        if len(samples) <=0: 
            with open(outfile, "a") as f:
                f.write(str(gr[0]) + "-" + str(gr[1]) + "\tNA\tNA\n")
            f.close()
            continue
        y_pred = model.predict(samples)
        y_true = dataset[(dataset["start"] == start) & (dataset["end"] == end)]["label"].tolist()
        # where label is positive replace with 1 and otherwise 0
        y_true = [1 if i == "positive" else 0 for i in y_true]
        accuracy = accuracy_score(y_true, y_pred)
        f1 = f1_score(y_true, y_pred)

        with open(outfile, "a") as f:
            f.write(str(gr[0]) + "-" + str(gr[1]) + "\t" + str(accuracy) + "\t" + str(f1) + "\n")
        f.close()

if __name__ == "__main__":
    sys.exit(main())