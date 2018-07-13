import numpy as np
from plastid import *
import dill
import matplotlib
import matplotlib.pyplot as plt

def allowed_transcript_ids(filename):
    with open(filename, "r") as gene_list_file:
        return gene_list_file.read().splitlines()


def fetch_vectors(filenames):
    gtf_coords_file = list(dill.load(open(gtf_assembly_pickle, "rb")))
    allowed_ids = set(allowed_transcript_ids(gene_set_path))
    alignments = BAMGenomeArray(filenames, mapping=FivePrimeMapFactory())
    print("Genomes loaded!")

    count_vectors = []
    filler_array = np.zeros(300, dtype=int)

    for transcript in gtf_coords_file:
        if transcript.get_length() > 400 and transcript.get_name() in allowed_ids:            
            value_array = transcript.get_counts(alignments)
            if np.sum(value_array) < 50:
                continue
            
            value_array = np.concatenate((value_array, filler_array))
            count_vectors.append(value_array[:600])

    print("Vectors retrieved!")
    vector_array = np.vstack(count_vectors)
    outlier = np.percentile(vector_array[vector_array > 1], 99.9)
    vector_array[vector_array > outlier] = outlier
    print("Outliers handled! Values above %f considered outliers!" % outlier)

    metagene = vector_array.sum(axis=0)
    metagene_stack = np.reshape(metagene, (200, 3))
    metagene_stack = np.hsplit(metagene_stack, 3)

    frames = []
    for arr in metagene_stack:
        frame_vector = np.concatenate((np.zeros(2), arr.T[0], np.zeros(2)))
        window_iter = (np.sum(frame_vector[i:i+5]) for i in range(len(frame_vector[2:-3])))
        frames.append(np.fromiter(window_iter, dtype=float))
    
    print("Fames split")
    return frames


def plot_results_coding():
    frames = fetch_vectors(bam_files_merged)
    frame_0, frame_1, frame_2 = frames[0], frames[1], frames[2]

    plt.grid(True, alpha=0.3)
    plt.plot(frame_0, linewidth=1, color="red", label="Frame 0" )
    plt.plot(frame_1, linewidth=1, color="blue", label ="Frame +1 (main)")
    plt.plot(frame_2, linewidth=1, color="green", label="Frame -1")
    plt.savefig("frames_sliding_t0_salt_all_reps_head.pdf")
    plt.close()


gene_set_path = "/home/maria/Documents/pelechanolab/scripts/coding_list.txt"
genome_fasta_path = "/home/maria/Documents/pelechanolab/data/R64e/genome.fa"
gtf_assembly_pickle = "gtf_assembled.sav"

bam_file_path = "/home/maria/Documents/pelechanolab/data/samples/salt_stress/cutadapt/umitools/star/dedup/wt-0-R1_S1.fastqAligned.sortedByCoord.out.bam"
#bam_file_path_2 = "/home/maria/Documents/pelechanolab/data/samples/salt_stress/cutadapt/umitools/star/dedup/wt-30-R1_S2.fastqAligned.sortedByCoord.out.bam"
#bam_file_path_3 = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment_3/cutadapt/umitools/star/dedup/BY4741-t0-3_S15.fastqAligned.sortedByCoord.out.bam"

bam_files_merged = [bam_file_path]

offset = 50

plot_results_coding()
