import os
import pandas as pd
import random
import editdistance

def check_if_different(read1, read2, threshold= 5):
    length1 = len(read1)
    length2 = len(read2)

    min_length = min(length1, length2)
    read1 = read1[:min_length]
    read2 = read2[:min_length]

    dist = editdistance.eval(read1, read2)

    dist = dist / min_length 
    dist = dist * 100

    if dist > threshold:
        return True
    
    return False

def sufficient_overlap(n_bases_overlap: int, bed1_row, bed2_row) -> bool:
    s1 = bed1_row.iloc[1]
    e1 = bed1_row.iloc[2]
    s2 = bed2_row.iloc[1]
    e2 = bed2_row.iloc[2]
    # check if the reads overlap
    if s1 < e2 and e1 > s2:
        # check if the reads overlap by at least n_bases_overlap
        # 1. Partial overlap
        if s1 < e2 - n_bases_overlap or e1 > s2 + n_bases_overlap:
            return True
        # 2. Full overlap
        if e1 - s1 > n_bases_overlap or e2 - s2 > n_bases_overlap:
            return True
    return False

def find_overlapping_reads(bed_file, bed_file2, fastq1, fastq2, label, outdir = "data/data/amplicon_lumc/reads"):
    try:
        bed1 = pd.read_csv(bed_file, sep="\t", header=None)
    except:
        print("bed_file: ", bed_file)
        raise Exception("bed_file not found")
    bed1.columns = ["ref", "start", "end", "read_name", "score", "strand"]

    bed2 = pd.read_csv(bed_file2, sep="\t", header=None)
    bed2.columns = ["ref", "start", "end", "read_name", "score", "strand"]
    fastq1_parsed = parse_fastq(fastq1)
    fastq2_parsed = parse_fastq(fastq2)
    # group reads by start position
    bed1_grouped = bed1.groupby("start")
    bed2_grouped = bed2.groupby("start")

    # get all unique start positions
    start_positions1 = bed2_grouped.groups.keys()
    start_positions2 = bed1_grouped.groups.keys()

    # find common start positions
    common_start_positions = set(start_positions1).intersection(set(start_positions2))
    
    result = set()
    # find reads that overlap
    for start_pos in common_start_positions: 
        check = True
        patience = 0
        while check and patience < 3:
            patience += 1
            id1 = bed1.groupby('start', group_keys=False).apply(pd.DataFrame.sample, frac=1).where(lambda x : x["start"] == start_pos).dropna().sample(1)["read_name"].values[0]
            id2 = bed2.groupby('start', group_keys=False).apply(pd.DataFrame.sample, frac=1).where(lambda x : x["start"] == start_pos).dropna().sample(1)["read_name"].values[0]
            if id1 != id2:
                id1 = id1.replace("/", "_")
                read1 = fastq1_parsed[id1]
                id2 = id2.replace("/", "_")
                read2 = fastq2_parsed[id2]
                if check_if_different(read1, read2):
                    check = False
                    result.add(frozenset([id1, id2, label]))

                    # write reads to file
                    with open(os.path.join(outdir, id1 + ".fastq"), "w") as f:
                        f.write("@" + id1 + "\n")
                        f.write(read1)
                    
                    with open(os.path.join(outdir, id2 + ".fastq"), "w") as f:
                        f.write("@" + id2 + "\n")
                        f.write(read2)

                    
                
                     
             
    return result

def parse_fastq(fastq_file):

    reads = {}

    with open(fastq_file, "r") as f:
        lines = f.readlines()

    lines = [x.strip() for x in lines]

    l = 0
    while(l != len(lines)):
        read_name = lines[l][1:].replace("/", "_")
        l+=1
        read = lines[l]
        l+=3
        reads[read_name] = read
      
    return reads

##############################################################################################################
n = 1
produce_alignment = False

wuhan_like_1 = ("data/data/amplicon_lumc/trimmed_01_100.fastq", "wuhan_like_1")
wuhan_like_2 = ("data/data/amplicon_lumc/trimmed_02_100.fastq", "wuhan_like_2")

omicron_1 = ("data/data/amplicon_lumc/trimmed_08_0.fastq", "omicron_1")
omicron_2 = ("data/data/amplicon_lumc/trimmed_09_0.fastq", "omicron_2")

data = [wuhan_like_1, wuhan_like_2, omicron_1, omicron_2]

if produce_alignment:
    for i in range(len(data)):
        os.system("minimap2 -ax map-hifi data/data/SARS-CoV-2-NC_045513.fa " + data[i][0] + " > data/data/amplicon_lumc/aln_" + data[i][1] + ".sam")
        os.system("samtools view -@ n -Sb -o data/data/amplicon_lumc/aln_" + data[i][1] + ".bam data/data/amplicon_lumc/aln_" + data[i][1] + ".sam")
        os.system("samtools sort -@ n -o data/data/amplicon_lumc/aln_" + data[i][1] + ".sorted.bam data/data/amplicon_lumc/aln_" + data[i][1] + ".bam")
        os.system("samtools index data/data/amplicon_lumc/aln_" + data[i][1] + ".sorted.bam")
        os.system("bedtools bamtobed -i data/data/amplicon_lumc/aln_" + data[i][1] + ".sorted.bam > data/data/amplicon_lumc/aln_" + data[i][1] + ".sorted.bam.bed")

bedfile_wuhan1 = "data/data/amplicon_lumc/aln_wuhan_like_1.sorted.bam.bed"
bedfile_wuhan2 = "data/data/amplicon_lumc/aln_wuhan_like_2.sorted.bam.bed"
bedfile_omicron1 = "data/data/amplicon_lumc/aln_omicron_1.sorted.bam.bed"
bedfile_omicron2 = "data/data/amplicon_lumc/aln_omicron_2.sorted.bam.bed"


count = 0
all_pairs = set()

#clear outdir 
for f in os.listdir("data/data/amplicon_lumc/reads"):
    os.remove(os.path.join("data/data/amplicon_lumc/reads", f))

while len(all_pairs) < n:

    pairs_wuhan = find_overlapping_reads(bedfile_wuhan1, bedfile_wuhan2, wuhan_like_1[0], wuhan_like_2[0], label="positive")
    all_pairs.update(pairs_wuhan)
    print("wuhan12: ", len(pairs_wuhan))
    pairs_wuhan1 = find_overlapping_reads(bedfile_wuhan1, bedfile_wuhan1, wuhan_like_1[0], wuhan_like_1[0], label="positive")
    all_pairs.update(pairs_wuhan1)
    print("wuhan11: ", len(pairs_wuhan1))
    pairs_wuhan2 = find_overlapping_reads(bedfile_wuhan2, bedfile_wuhan2, wuhan_like_2[0], wuhan_like_2[0], label="positive")
    all_pairs.update(pairs_wuhan2)
    print("wuhan22: ", len(pairs_wuhan2))
    pairs_omicron = find_overlapping_reads(bedfile_omicron1, bedfile_omicron2, omicron_1[0], omicron_2[0], label="positive")
    all_pairs.update(pairs_omicron)
    print("omicron12: ", len(pairs_omicron))
    pairs_omicron1 = find_overlapping_reads(bedfile_omicron1, bedfile_omicron1, omicron_1[0], omicron_1[0], label="positive")
    all_pairs.update(pairs_omicron1)
    print("omicron11: ", len(pairs_omicron1))
    pairs_omicron2 = find_overlapping_reads(bedfile_omicron2, bedfile_omicron2, omicron_2[0], omicron_2[0], label="positive")
    all_pairs.update(pairs_omicron2)
    print("omicron22: ", len(pairs_omicron2))
    pairs_wuhan_omicron = find_overlapping_reads(bedfile_wuhan1, bedfile_omicron1, wuhan_like_1[0], omicron_1[0], label="negative")
    all_pairs.update(pairs_wuhan_omicron)
    print("wuhan1_omicron1: ", len(pairs_wuhan_omicron))
    pairs_wuhan_omicron2 = find_overlapping_reads(bedfile_wuhan2, bedfile_omicron2, wuhan_like_2[0], omicron_2[0], label="negative")
    all_pairs.update(pairs_wuhan_omicron2)
    print("wuhan2_omicron2: ", len(pairs_wuhan_omicron2))
    pairs_wuhan_omicron12 = find_overlapping_reads(bedfile_wuhan1, bedfile_omicron2, wuhan_like_1[0], omicron_2[0], label="negative")
    all_pairs.update(pairs_wuhan_omicron12)
    print("wuhan1_omicron2: ", len(pairs_wuhan_omicron12))
    pairs_wuhan_omicron21 = find_overlapping_reads(bedfile_wuhan2, bedfile_omicron1, wuhan_like_2[0], omicron_1[0], label="negative")
    all_pairs.update(pairs_wuhan_omicron21)
    print("wuhan2_omicron1: ", len(pairs_wuhan_omicron21))


all_pairs = [list(x) for x in all_pairs]
print ("Number of pairs: ", len(all_pairs))
# remove file if already exists
if os.path.exists("data/data/amplicon_lumc/reads/amplicon_lumc_pairs.tsv"):
    os.remove("data/data/amplicon_lumc/reads/amplicon_lumc_pairs.tsv")
# save pairs to file
with open("data/data/amplicon_lumc/reads/amplicon_lumc_pairs.tsv", "w") as f:
    for pair in all_pairs:
        f.write(pair[0] + "\t" + pair[1] + "\t" + pair[2] + "\n")
    






   




