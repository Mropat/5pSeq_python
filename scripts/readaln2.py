import numpy as np
import pysam
import gffutils
import re
from plastid import *
        

def aln_depths_region(filename):
    

    with pysam.AlignmentFile(filename, "r") as bamload:
        for seg in bamload.fetch("I"):
            print(seg)


if __name__ == "__main__":
    coordfile = "data/ribo.txt"
    bam_file = "data/samples/2/post-proc/by47.sam"
    aln_depths_region(bam_file)

