import argparse
import sys
import os



parser = argparse.ArgumentParser()
parser.add_argument("--sample_dir", required=True)
parser.add_argument("--genome_fasta", required=True)
parser.add_argument("--annotation_file", required=True)
parser.add_argument("--offset", type=int, default=50)
parser.add_argument("--cores", type=int, default=2)
parser.add_argument("--output_dir", required=True)
args = parser.parse_args()
print(args.sample_dir)


for filename in os.listdir(args.sample_dir):
    if filename.endswith(".md"):
        print (os.path.join(args.sample_dir, filename))