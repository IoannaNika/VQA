import os
import sys
import random 
import pandas as pd



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
# Returns: a bed file and a fastq file
############################################################################################################################################################################
def sample_fastq_bed_file_from_lineage(directory):
    #find all fastq files in directory
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
# Input: two bed file for a fastq file. If you want to sample from the same 
#       genome, then pass the same bed file twice.
# Returns: two reads that overlap
############################################################################################################################################################################
def find_overlapping_reads(bed_file, bed_file2, fastq_file, fastq_file2):
    # find all reads in bed_file that overlap with reads in bed_file2
    # return the number of overlapping reads
    bed1 = pd.read_csv(bed_file, sep="\t", header=None)
    bed1.columns = ["ref", "start", "end", "read_name", "score", "strand"]

    bed2 = pd.read_csv(bed_file2, sep="\t", header=None)
    bed2.columns = ["ref", "start", "end", "read_name", "score", "strand"]

    # create a dictionary of read names and the reads they overlap with
    read_dict = {}

    for i in range(len(bed1)):
        for j in range(len(bed2)):
            if bed1.iloc[i, 3] != bed2.iloc[j, 3]:
                # start_i < end_j  and end_i > start_j
                if bed1.iloc[i, 1] < bed2.iloc[j, 2] and bed1.iloc[i, 2] > bed2.iloc[j, 1]:
                    if bed1.iloc[i, 3] in read_dict:
                        read_dict[bed1.iloc[i, 3]].add(bed2.iloc[j, 3])
                    else:
                        read_dict[bed1.iloc[i, 3]] = {bed2.iloc[j, 3]}

    # randomly sample a read from the dictionary
    read1 = random.choice(list(read_dict.keys()))
    read2 = random.choice(list(read_dict[read1]))


    # find the read1 in the first fastq file
    read1_seq = parse_fastq(fastq_file)[read1][0]
    read2_seq = parse_fastq(fastq_file2)[read2][0]
    return read1_seq, read2_seq, read1, read2

############################################################################################################################################################################
# Input: same_lineage (bool), lineage1 (str), lineage2 (str)
# Returns: two reads that overlap and their respective read and genome ids
############################################################################################################################################################################
def sample_negative(same_lineage: bool, lineage1: str, lineage2: str = ""):
    data_dir =  "data/data/hcov_global_2023-11-16_09-28/split_fasta"
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

        while(id == id2):
            id2, fastq2, bed2 = sample_fastq_bed_file_from_lineage(data_dir1)

        read1_seq, read2_seq, read1_id, read_2_id  = find_overlapping_reads(bed, bed2, fastq, fastq2)

        return read1_seq, read2_seq, read1_id, read_2_id, id, id2
    
    # sample from two different lineages
    data_dir2 = data_dir + "/"+ lineage2 + "/simulated_reads/"

    id, fastq, bed = sample_fastq_bed_file_from_lineage(data_dir1)
    id2, fastq2, bed2 = sample_fastq_bed_file_from_lineage(data_dir2)

    read1_seq, read2_seq, read1_id, read_2_id  = find_overlapping_reads(bed, bed2, fastq, fastq2)

    return read1_seq, read2_seq, read1_id, read_2_id

############################################################################################################################################################################
# Input: lineage (str). Optional. If not provided, sample from any lineage
# Returns: two reads that overlap and their respective read and genome ids
############################################################################################################################################################################
def sample_positive(lineage: str = ""):

    data_dir =  "data/data/hcov_global_2023-11-16_09-28/split_fasta"

    if lineage not in os.listdir(data_dir):
        raise Exception ("Lineage not found in data directory")
    
    if lineage == "":
        # take subdirectories
        lineages = [x[0] for x in os.walk(data_dir)]
        random.shuffle(lineages)
        lineage_to_sample = lineages[0]
        data_dir += "/" + lineage_to_sample + "/simulated_reads/"
    else:
        data_dir += "/" + lineage + "/simulated_reads/"

    # randomly sample a fastq file
    id, fastq, bed = sample_fastq_bed_file_from_lineage(data_dir)

    # find two reads that overlap
    read1_seq, read2_seq, read1_id, read_2_id  = find_overlapping_reads(bed, bed, fastq, fastq)

    return read1_seq, read2_seq, read1_id, read_2_id, id, id


# def main():
    # # example of sample positive
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




# if __name__ == "__main__":
#     sys.exit(main())
