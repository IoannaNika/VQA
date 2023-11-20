import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="simulate long reads from genomes")
    parser.add_argument('--dir', dest = 'dir', required=True, type=str, help="directory with split fasta files. Expected structure: <dir>/<pango_lineage>/<fasta_file>")
    parser.add_argument('--depth', dest = 'depth', default= 20, required=False, type=int, help="depth of coverage")
    parser.add_argument('--error_model', dest = 'error_model', default="R103", required=False, type=str, help="options for error simulation: Nanopore (starts with N), PacBio (starts with P): P4C2, P5C3, P6C4, and R94, R95, R103: Options appear in the order of their release date.")
    args = parser.parse_args()

    # go through all .fasta files in the data directory provided
    for subdir, dirs, files in os.walk(args.dir):
        for file in files:
            if file.endswith(".fasta"):
                # create output directory
                output_dir = args.dir + "/" + subdir.split("/")[-1] + "/simulated_reads/"
                if not os.path.exists(output_dir):
                    os.mkdir(output_dir)
                output_dir += file.split(".")[0]
                if not os.path.exists(output_dir):
                    os.mkdir(output_dir)
                else:
                    # remove all files in output directory
                    os.system("rm " + output_dir + "/*")

                # cd in the output directory
                os.chdir(output_dir)
                    
                if args.error_model.startswith("P"):
                    # simulate pacbio reads
                    os.system("pbsim --depth " + str(args.depth) + " --hmm_model " + "../../../../../simulation_files/" + args.error_model + ".model --prefix " +  file.split(".")[0]  + " --id-prefix " +  file.split(".")[0] + " ../../" + file)
                else: 
                    # nanopore reads
                    os.system("pbsim --depth " + str(args.depth) + " --hmm_model " + "../../../../../simulation_files/" +  args.error_model + ".model --prefix " +  file.split(".")[0]  + " --id-prefix " + file.split(".")[0] + " --difference-ratio 23:31:46 ../../" + file)
                
                # cd back to the original directory
                os.chdir("../../../../../../../")
       

if __name__ == "__main__":
    sys.exit(main())