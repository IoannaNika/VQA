import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--references_file', dest = 'references_file', required=True, type=str, help="Folder with templates")
    parser.add_argument('--data_dir', dest = 'data_dir', required=True, type=str, help="Folder with data, will be used as output folder too")
    args = parser.parse_args()

    outdir = os.path.join(args.data_dir, "quast")
    os.makedirs(outdir, exist_ok=True)
    consensus = os.path.join(args.data_dir, "consensus.fasta")

    # run quast
    os.system("metaquast -o " + outdir + " -r " + args.references_file + "--ambiguity-usage all " + consensus)

if __name__ == "__main__":
    sys.exit(main())