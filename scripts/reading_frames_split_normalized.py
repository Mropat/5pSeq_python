import numpy as np
from plastid import *
import dill
import matplotlib
import matplotlib.pyplot as plt
from time import time

def timeit(f):
    def timer(*args, **kwargs):
        st = int(time())
        print("Start!")
        f(*args, **kwargs)
        et = int(time())
        print('%s took %i seconds to complete! (args: %s)' % (f.__name__, et - st, ', '.join(map(str, args))))
    return timer


def allowed_transcript_ids(filename):
    with open(filename, "r") as gene_list_file:
        return gene_list_file.read().splitlines()


def fetch_vectors(filenames):
    gtf_coords_file = list(dill.load(open(gtf_assembly_pickle, "rb")))
    allowed_ids = set(allowed_transcript_ids(gene_set_path))
    alignments = BAMGenomeArray(filenames, mapping=FivePrimeMapFactory())
    print("Genomes loaded!")

    count_vectors = []

    for transcript in gtf_coords_file:
        if transcript.get_length() > 1800  and transcript.get_name() in allowed_ids:            
            value_array = transcript.get_counts(alignments)
            if np.sum(value_array[:1800]) > 1:
                count_vectors.append(value_array[:1800])

    print("Vectors retrieved!")
    vector_array = np.vstack(count_vectors)
    vector_normsum = np.sum(vector_array, axis=1)
    vector_array_normalized = vector_array / vector_normsum[:, np.newaxis]

    metagene = vector_array_normalized.sum(axis=0)
    metagene_stack = np.reshape(metagene, (-1, 3))
    metagene_stack = np.hsplit(metagene_stack, 3)

    frames = []
    for arr in metagene_stack:
        frame_vector = np.concatenate((np.zeros(2), arr.T[0], np.zeros(2)))
        window_iter = (np.sum(frame_vector[i:i+5]) for i in range(len(frame_vector[2:-3])))
        frames.append(np.fromiter(window_iter, dtype=float))
    
    print("Frames split")
    return frames

@timeit
def plot_results_coding(bam_files_merged):
    frames = fetch_vectors(bam_files_merged)
    frame_0, frame_1, frame_2 = frames[0], frames[1], frames[2]

    plt.grid(True, alpha=0.3)
    plt.plot(frame_0, linewidth=1, color="red", label="Frame 0" )
    plt.plot(frame_1, linewidth=1, color="blue", label ="Frame +1 (main)")
    plt.plot(frame_2, linewidth=1, color="green", label="Frame -1")
    plt.savefig("frames_sliding_%s_all_reps_longrange_normsum.pdf"% setarg)
    plt.close()


gene_set_path = "/home/maria/Documents/pelechanolab/scripts/coding_list.txt"
genome_fasta_path = "/home/maria/Documents/pelechanolab/data/R64e/genome.fa"
gtf_assembly_pickle = "gtf_assembled.sav"

def exposure_time(arg):
    
    if arg == "t0":           
        bam_file_path = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment/cutadapt/umitools/star/dedup/BY4741-t0-2_S23.fastqAligned.sortedByCoord.out.bam"
        bam_file_path_2 = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment_2/cutadapt/umitools/star/dedup/BY4741-t0-1_S7.fastqAligned.sortedByCoord.out.bam"
        bam_file_path_3 = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment_3/cutadapt/umitools/star/dedup/BY4741-t0-3_S15.fastqAligned.sortedByCoord.out.bam"

    elif arg == "t5":        
        bam_file_path = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment/cutadapt/umitools/star/dedup/BY4741-t5-2_S24.fastqAligned.sortedByCoord.out.bam"
        bam_file_path_2 = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment_2/cutadapt/umitools/star/dedup/BY4741-t5-1_S8.fastqAligned.sortedByCoord.out.bam"
        bam_file_path_3 = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment_3/cutadapt/umitools/star/dedup/BY4741-t5-3_S16.fastqAligned.sortedByCoord.out.bam"

    elif arg == "t15":        
        bam_file_path = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment/cutadapt/umitools/star/dedup/BY4741-t15-2_S25.fastqAligned.sortedByCoord.out.bam"
        bam_file_path_2 = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment_2/cutadapt/umitools/star/dedup/BY4741-t15-1_S9.fastqAligned.sortedByCoord.out.bam"
        bam_file_path_3 = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment_3/cutadapt/umitools/star/dedup/BY4741-t15-3_S17.fastqAligned.sortedByCoord.out.bam"

    elif arg == "t30":        
        bam_file_path = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment/cutadapt/umitools/star/dedup/BY4741-t30-2_S26.fastqAligned.sortedByCoord.out.bam"
        bam_file_path_2 = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment_2/cutadapt/umitools/star/dedup/BY4741-t30-1_S10.fastqAligned.sortedByCoord.out.bam"
        bam_file_path_3 = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment_3/cutadapt/umitools/star/dedup/BY4741-t30-3_S18.fastqAligned.sortedByCoord.out.bam"

    return [bam_file_path, bam_file_path_2, bam_file_path_3]

setarglist = ["t0", "t5", "t15", "t30"]
for setarg in setarglist:
    plot_results_coding(exposure_time(setarg))
