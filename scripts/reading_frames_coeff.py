import numpy as np
from plastid import *
import dill
import matplotlib
import matplotlib.pyplot as plt

def allowed_transcript_ids(filename):
    with open(filename, "r") as gene_list_file:
        return gene_list_file.read().splitlines()


def fetch_vectors(filename):
    gtf_coords_file = list(dill.load(open(gtf_assembly_pickle, "rb")))
    allowed_ids = set(allowed_transcript_ids(gene_set_path))
    alignments = BAMGenomeArray(filename, mapping=FivePrimeMapFactory())

    for transcript in gtf_coords_file:
        if transcript.get_name() in allowed_ids:            
            value_array = transcript.get_counts(alignments)
            value_array = np.reshape(value_array, (int(value_array.size/3), 3))
            value_sums = np.sum(value_array, axis=1)
            value_array = value_array / value_sums [:, np.newaxis]
#            value_array = np.nan_to_num(value_array)
            print(value_array) 
                               



gene_set_path = "/home/maria/Documents/pelechanolab/scripts/coding_list.txt"
genome_fasta_path = "/home/maria/Documents/pelechanolab/data/R64e/genome.fa"
gtf_assembly_pickle = "gtf_assembled.sav"
bam_file_path = "/home/maria/Documents/pelechanolab/data/samples/20160909YP_08-39666715/cutadapt/umitools/star/dedup/BY4741-t5-1_S8_L002_R1_001.fastqAligned.sortedByCoord.out.bam"
fetch_vectors(bam_file_path)