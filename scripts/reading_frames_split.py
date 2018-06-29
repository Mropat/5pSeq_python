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

    count_vectors = []
    filler_array = np.zeros(300, dtype=int)
    for transcript in gtf_coords_file:
        if transcript.get_length() > 400 and transcript.get_name() in allowed_ids:            
            value_array = transcript.get_counts(alignments)
            value_array = np.concatenate((value_array, filler_array))
            count_vectors.append(value_array[:600])

    vector_array = np.vstack(count_vectors)
    outlier = np.percentile(vector_array[vector_array > 2], 99.0)
    vector_array[vector_array > outlier] = outlier

    metagene = vector_array.sum(axis=0)
    metagene_stack = np.reshape(metagene, (200, 3))
    metagene_stack = np.hsplit(metagene_stack, 3)

    frames = []
    for arr in metagene_stack:
        frame_vector = np.concatenate((np.zeros(2), arr.T[0], np.zeros(2)))
        window_iter = (np.sum(frame_vector[i:i+5]) for i in range(len(frame_vector[2:-3])))
        sliding_vector = np.fromiter(window_iter, dtype=float)
        frames.append(sliding_vector)
    return frames


def plot_results_coding():
    frames = fetch_vectors(bam_file_path)
    frame_0, frame_1, frame_2 = frames[0], frames[1], frames[2]

    plt.grid(True, alpha=0.3)
    plt.plot(frame_0, linewidth=0.5, color="blue", label="Frame 0" )
    plt.plot(frame_1, linewidth=1, color="red", label ="Frame +1 (main)")
    plt.plot(frame_2, linewidth=0.5, color="green", label="Frame -1")
    plt.savefig("frames_sliding_t5_outliers.pdf")
    plt.close()


gene_set_path = "/home/maria/Documents/pelechanolab/outlier_ids_t5.txt"
genome_fasta_path = "/home/maria/Documents/pelechanolab/data/R64e/genome.fa"
gtf_assembly_pickle = "gtf_assembled.sav"
bam_file_path = "/home/maria/Documents/pelechanolab/data/samples/20160909YP_08-39666715/cutadapt/umitools/star/dedup/BY4741-t5-1_S8_L002_R1_001.fastqAligned.sortedByCoord.out.bam"
offset = 50

plot_results_coding()
