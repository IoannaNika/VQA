import os
import sys
import random 
import pandas as pd



def hamming_distance(seq1: str, seq2: str) -> int:
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

############################################################################################################################################################################
# Input: path to fastq file
# Returns: a dictionary with the fastq reads. Read names are the keys,
#         and the values are a list of the read and the quality score
############################################################################################################################################################################
def parse_fastq(fastq_file):

    reads = {}

    with open(fastq_file, "r") as f:
        lines = f.readlines()

    lines = [x.strip() for x in lines]

    l = 0
    while(l!= len(lines)):
        read_name = lines[l][1:]
        l+=1
        read = lines[l]
        l+=2
        quality_score = lines[l]
        l+=1
        reads[read_name] = [read, quality_score]

    return reads
    

############################################################################################################################################################################
# Input: path to fastq file, and read type (ONT only for now) and virus (SARS-COV-2 only for now)
# Returns: a bed file for the fastq file
############################################################################################################################################################################
def produce_bed_file(fastqid, fastq_file_path, read_type = "ONT", Virus = "SARS-COV-2"):
    # enter the directory of the fastq file
    fastq_file_dir = fastq_file_path.split("/")[:-1]
    fastq_file_dir = "/".join(fastq_file_dir)
    os.chdir(fastq_file_dir)
    os.system("minimap2 -ax map-ont ../../../../../../../data/data/SARS-CoV-2-NC_045513.fa {} > {}.sam".format(fastq_file_path.split("/")[-1], fastqid))
    os.system("samtools view -@ n -Sb -o {}.bam {}.sam".format(fastqid, fastqid))

    os.system("samtools sort -@ n -o {}.sorted.bam {}.bam".format(fastqid, fastqid))
    os.system("samtools index {}.sorted.bam".format(fastqid))
    os.system("bedtools bamtobed -i {}.sorted.bam > {}.bed".format(fastqid, fastqid)) 

    # return to original directory
    os.chdir("../../../../../../../")
    return fastq_file_dir + "/" + fastqid + ".bed"


############################################################################################################################################################################
# Input: directory of lineage to sample from
# directory given must include the lineage directory
# e.g. data/data/hcov_global_2023-11-16_09-28/split_fasta/AY.3/simulated_reads
# expected full directory structure: 
# src/vqa/data/data/hcov_global_2023-11-16_09-28/split_fasta/AY.3/simulated_reads/EPI_ISL_5846086/EPI_ISL_5846086_0001.fastq
# Returns: a bed file (path) and a fastq file (path)
############################################################################################################################################################################
def sample_fastq_bed_file_from_lineage(directory):
    # find all fastq files in directory
    subdirs = [x for x in os.listdir(directory)]
    # randomly sample a subdir 
    # shuffle subdirs
    random.shuffle(subdirs)
    subdir = subdirs[0] # also the genome id
    # find all fastq files in subdir
    fastq_file = [x.split("/")[-1] for x in os.listdir(directory + "/" + subdir) if x.endswith(".fastq")]
    fastq_file =  directory + subdir + "/" + fastq_file[0]
    # check if there is a .bed file created for the fastq file
    # if not, create the .bed file
    bed_file = [x for x in os.listdir(directory + "/" + subdir) if x.endswith(".bed")]

    if len(bed_file) == 0:
        # create bed file
       bed_file = produce_bed_file(subdir, fastq_file)
    else:
        bed_file = directory + "/" + subdir + "/" + bed_file[0]

    return subdir, fastq_file, bed_file

############################################################################################################################################################################
# Input: n_bases_overlap (int), .bed file 1 (pd.Dataframe), .bed file 2 (pd.Dataframe) for each fastq file
# Returns: True if the reads overlap by at least n_bases_overlap, False otherwise
############################################################################################################################################################################
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

############################################################################################################################################################################
# Input: two bed file for a fastq file. If you want to sample from the same 
#       genome, then pass the same bed file twice.
# Returns: two reads that overlap
############################################################################################################################################################################
def find_overlapping_reads(bed_file, bed_file2, fastq_file, fastq_file2, n_bases_overlap = 100):
    # find all reads in bed_file that overlap with reads in bed_file2
    # return the number of overlapping reads
    try:
        bed1 = pd.read_csv(bed_file, sep="\t", header=None)
    except:
        print("bed_file: ", bed_file)
    bed1.columns = ["ref", "start", "end", "read_name", "score", "strand"]

    bed2 = pd.read_csv(bed_file2, sep="\t", header=None)
    bed2.columns = ["ref", "start", "end", "read_name", "score", "strand"]

    # create a dictionary of read names and the reads they sufficiently overlap with
    read_dict = {}
    for i in range(len(bed1)):
        for j in range(len(bed2)):
            # in case reads are from the same genome, don't compare the same read
            if bed1.iloc[i, 3] != bed2.iloc[j, 3]:
                if sufficient_overlap(n_bases_overlap, bed1.iloc[i], bed2.iloc[j]): 
                    if bed1.iloc[i, 3] in read_dict:
                        read_dict[bed1.iloc[i, 3]].add(bed2.iloc[j, 3])
                    else:
                        read_dict[bed1.iloc[i, 3]] = {bed2.iloc[j, 3]}


    # randomly sample a read from the dictionary
    read1 = random.choice(list(read_dict.keys())) # id of read1
    # randomly sample a read that overlaps with read1
    read2 = random.choice(list(read_dict[read1])) # id of read2


    # find the read1 in the first fastq file
    read1_seq = parse_fastq(fastq_file)[read1][0] # read1 sequence
    read2_seq = parse_fastq(fastq_file2)[read2][0] # read2 sequence
    return read1_seq, read2_seq, read1, read2

############################################################################################################################################################################
# Input: maf file
# Returns: original read sequence
############################################################################################################################################################################
def parse_maf_file(maf_file, read_id):
    maf_dict = {}
    with open(maf_file, "r") as f:
        lines = f.readlines()

    cnt = 0
    while cnt < len(lines):
        if lines[cnt].startswith("a"):
            next_line = lines[cnt+1].strip().split(" ")
            original_read_seq = next_line[6].replace("-", "")
            next_line = lines[cnt+2].strip().split(" ")
            sim_read_id = next_line[1]
            maf_dict[sim_read_id] = original_read_seq
            cnt += 4
    try: 
        maf_dict[read_id]
    except:
        print(maf_dict.keys())
    return maf_dict[read_id]

############################################################################################################################################################################
# Input: maf file 1, maf file 2, read1_id, read2_id, bed file 1, bed file 2, distance_threshold
# Returns: True if the reads are sufficiently different, False otherwise
############################################################################################################################################################################
def check_negative_sample(maf_file_1:str, maf_file_2:str, read1_id:str, read2_id:str, bed_file_1:str, bed_file_2:str, distance_threshold:float = 0.0):
    # get the original read sequences
    read1_seq = parse_maf_file(maf_file_1, read1_id)
    read2_seq = parse_maf_file(maf_file_2, read2_id)
    # get the read sequences from the bed files
    bed1 = pd.read_csv(bed_file_1, sep="\t", header=None)
    bed1.columns = ["ref", "start", "end", "read_name", "score", "strand"]
    start_1 = bed1[bed1["read_name"] == read1_id].iloc[0, 1]
    end_1 = bed1[bed1["read_name"] == read1_id].iloc[0, 2]

    bed2 = pd.read_csv(bed_file_2, sep="\t", header=None)
    bed2.columns = ["ref", "start", "end", "read_name", "score", "strand"]
    start_2 = bed2[bed2["read_name"] == read2_id].iloc[0, 1]
    end_2 = bed2[bed2["read_name"] == read2_id].iloc[0, 2]

   
    # make indices 0-based
    # find earlier start
    if start_1 < start_2:
        start = start_1
    else:
        start = start_2
    
    start_1 -= start
    end_1 -= start
    start_2 -= start
    end_2 -= start

    
    overlap_start = max(start_1, start_2)
    overlap_end = min(end_1, end_2)
    
    # calculate hamming distance between the two reads in the overlap region
    hamming_dist = hamming_distance(read1_seq[overlap_start:overlap_end], read2_seq[overlap_start:overlap_end])

    # calculate dissimilarity percentage
    dissimilarity_percentage = hamming_dist / (overlap_end - overlap_start)

    # check if the dissimilarity percentage is less than the threshold
    if dissimilarity_percentage > distance_threshold:
        return True
    else:
        return False



############################################################################################################################################################################
# Input: same_lineage (bool), lineage1 (str), lineage2 (str)
# Returns: two reads that overlap and their respective read and genome ids
############################################################################################################################################################################
def sample_negative(same_lineage: bool, lineage1: str, lineage2: str = "", data_dir: str = "data/data/hcov_global_2023-11-16_09-28/split_fasta"):
    # take subdirectories 
    lineages = [x for x in os.listdir(data_dir)]

  
    if lineage1 not in lineages:
        raise Exception ("Lineage1 not found in data directory")
    
    if lineage2 == "" and same_lineage == False:
        raise Exception ("Cannot sample from two different lineages. Please provide lineage2")
    
    if  lineage2 != "" and lineage2 == True and lineage2 not in lineages:
        raise Exception ("Lineage2 not found in data directory")
            
    # sample from specific lineage
    data_dir1 = data_dir + "/" + lineage1 + "/simulated_reads/"

    if same_lineage == True: 
        # check that there is more than one fastq file in the directory
        ids = [x for x in os.listdir(data_dir1)]


        if len(ids) == 1:
            raise Exception ("Only one genome sequence for this lineage. Cannot sample from the same lineage")
        
        id, fastq, bed = sample_fastq_bed_file_from_lineage(data_dir1)
        id2, fastq2, bed2 = sample_fastq_bed_file_from_lineage(data_dir1)

        maf_file_1 = data_dir1 + id + "/" + id + "_0001.maf"
        
        while(id == id2):
            id2, fastq2, bed2 = sample_fastq_bed_file_from_lineage(data_dir1)

        maf_file_2 = data_dir1 + id2 + "/" + id2 + "_0001.maf"
        read1_seq, read2_seq, read1_id, read_2_id  = find_overlapping_reads(bed, bed2, fastq, fastq2)
        
        # check if the reads are sufficiently different
        while(check_negative_sample(maf_file_1, maf_file_2, read1_id, read_2_id, bed, bed2) == False):
            read1_seq, read2_seq, read1_id, read_2_id  = find_overlapping_reads(bed, bed2, fastq, fastq2)

        return read1_seq, read2_seq, read1_id, read_2_id, id, id2
    
    # sample from two different lineages
    data_dir2 = data_dir + "/"+ lineage2 + "/simulated_reads/"

    id, fastq, bed = sample_fastq_bed_file_from_lineage(data_dir1)
    id2, fastq2, bed2 = sample_fastq_bed_file_from_lineage(data_dir2)

    read1_seq, read2_seq, read1_id, read_2_id  = find_overlapping_reads(bed, bed2, fastq, fastq2)
    # check if the reads are sufficiently different 
    maf_file_1 = data_dir1 + id + "/" + id + "_0001.maf"
    maf_file_2 = data_dir2 + id2 + "/" + id2 + "_0001.maf"
    
    while(check_negative_sample(maf_file_1, maf_file_2, read1_id, read_2_id, bed, bed2) == False):
        read1_seq, read2_seq, read1_id, read_2_id  = find_overlapping_reads(bed, bed2, fastq, fastq2)

    return read1_seq, read2_seq, read1_id, read_2_id, id, id2

############################################################################################################################################################################
# Input: lineage (str). Optional. If not provided, sample from any lineage
# Returns: two reads that overlap and their respective read and genome ids
############################################################################################################################################################################
def sample_positive(lineage: str = "",  data_dir: str = "data/data/hcov_global_2023-11-16_09-28/split_fasta"):

    if lineage not in os.listdir(data_dir):
        raise Exception ("Lineage not found in data directory", lineage)
    
    if lineage == "":
        # take subdirectories
        lineages = [x for x in os.listdir(data_dir)]
        random.shuffle(lineages)
        lineage = lineages[0]
        data_dir += "/" + lineage + "/simulated_reads/"
    else:
        data_dir += "/" + lineage + "/simulated_reads/"

    # randomly sample a fastq file
    id, fastq, bed = sample_fastq_bed_file_from_lineage(data_dir)

    # find two reads that overlap
    read1_seq, read2_seq, read1_id, read_2_id  = find_overlapping_reads(bed, bed, fastq, fastq)

    return read1_seq, read2_seq, read1_id, read_2_id, id, id

############################################################################################################################################################################
# Input: max_len (int), read (str)
# Returns: padded read. If read is longer than max_len, then truncate the read.
############################################################################################################################################################################
def pad_read(max_len: int, read: str):
    if len(read) < max_len:
        read += "N" * (max_len - len(read))

    if len(read) > max_len:
        read = read[:max_len]    
    return read

############################################################################################################################################################################
# Input: read (str)
# Returns: one hot encoded read
############################################################################################################################################################################
def one_hot_encode(read: str):
    # A = [1, 0, 0, 0]
    # C = [0, 1, 0, 0]
    # G = [0, 0, 1, 0]
    # T = [0, 0, 0, 1]
    # N = [0, 0, 0, 0]

    base_to_encoding_dict = {"A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T": [0, 0, 0, 1], "N": [0, 0, 0, 0]}
    encoded_read = []
    for base in read:
        encoded_read.append(base_to_encoding_dict[base])
    
    # to dataframe
    encoded_read = pd.DataFrame(encoded_read)
    return encoded_read.T



# def main():
    # example of sample positive
    # read1, read2, read1_id, read2_id, id, id2 = sample_positive("XBL.2")
    # print("read1: ", read1)
    # print("read2: ", read2)
    # print("read1_id: ", read1_id)
    # print("read2_id: ", read2_id)
    # print("id: ", id)
    # print("id2: ", id2)

    # example of sample negative
    # read1, read2, read1_id, read2_id, id, id2 = sample_negative(True, "GN.4")
    # print("read1: ", read1)
    # print("read2: ", read2)
    # print("read1_id: ", read1_id)
    # print("read2_id: ", read2_id)
    # print("id: ", id)
    # print("id2: ", id2)

    # example of sample negative
    # read1, read2, read1_id, read2_id = sample_negative(False, "BE.1.1.1", "BE.1.1")
    # print("read1: ", read1)
    # print("read2: ", read2)
    # print("read1_id: ", read1_id)
    # print("read2_id: ", read2_id)

    # # example of pad read
    # read = "ATCG"
    # padded_read = pad_read(10, read)
    # print(padded_read)

    # # example of one hot encode
    # encoded_read = one_hot_encode(padded_read)
    # print(encoded_read.T)

    # parse maf file
    # 

# if __name__ == "__main__":
#     sys.exit(main())
