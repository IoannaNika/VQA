import argparse
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--consensus', dest = 'consensus', required=True, type=str, help="File with consensus sequences")
    parser.add_argument('--template_folder', dest = 'template_folder', required=True, type=str, help="Folder with templates")
    args = parser.parse_args()





if __name__ == "__main__":
    sys.exit(main())