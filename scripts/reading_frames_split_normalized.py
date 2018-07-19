# Script for plotting the ribosome protection frames

import numpy as np
from plastid import GTF2_TranscriptAssembler, GFF3_TranscriptAssembler, Transcript, GenomicSegment, BAMGenomeArray, FivePrimeMapFactory
import dill
import matplotlib
import matplotlib.pyplot as plt
from time import time
from multiprocessing import Process
import os


def allowed_transcript_ids(filename):
    with open(filename, "r") as gene_list_file:
        return gene_list_file.read().splitlines()


#Generates vectors of 5P coverages using mapping coordinates

def fetch_vectors(filenames):

    gtf_coords_file = list(dill.load(open(gtf_assembly_pickle, "rb")))
    allowed_ids = set(allowed_transcript_ids(gene_set_path))


    for transcript in gtf_coords_file:
        span = transcript.spanning_segment
        new_region = GenomicSegment(
            span.chrom, span.start - offset, span.start, span.strand)
        new_region_2 = GenomicSegment(
            span.chrom, span.end, span.end + offset, span.strand)
        transcript.add_segments(new_region, new_region_2)

         
    alignments = BAMGenomeArray(filenames, mapping=FivePrimeMapFactory())
    print("Genomes loaded!")


    except_count = 0
    count_vectors = []
    for transcript in gtf_coords_file:
        if transcript.attr.get('Name') in allowed_ids:
            try: 
                value_array = transcript.get_counts(alignments)
    #            if np.sum(value_array[-1800:]) > 1:
    #                count_vectors.append(np.concatenate((np.zeros(4000, dtype=int), value_array))[-1800:])
                if np.sum(value_array[:1800]) > 1:
                    count_vectors.append(np.concatenate((value_array, np.zeros(4000, dtype=int)))[:1800])
            except ValueError:
                except_count += 1


    vector_array = np.vstack(count_vectors)

    print("Vectors retrieved!")
    print("Removed %i transcripts!"%except_count)

#This normalizes vectors by read (optional) 
#    vector_normsum = np.sum(vector_array, axis=1)
#    vector_array_normalized = vector_array / vector_normsum[:, np.newaxis]



    metagene = vector_array.sum(axis=0)


    metagene_stack = np.reshape(metagene, (-1, 3))
    metagene_stack = np.hsplit(metagene_stack, 3)

    frames = []
    for arr in metagene_stack:
        frame_vector = np.concatenate((np.zeros(2), arr.T[0], np.zeros(2)))
        window_iter = (np.sum(frame_vector[i:i+5])
                       for i in range(len(frame_vector[2:-3])))
        frames.append(np.fromiter(window_iter, dtype=float))

    print("Frames split")
    return frames


def plot_results_coding(bam_files_merged, filename):
    ''' 
    Plots each frame. If sequence boundries are not modified, frame 1 is the one producing
    initiation peak at -17
    '''
    
    frames = fetch_vectors(bam_files_merged)

    plt.grid(True, alpha=0.3)
    plt.plot(frames[0], linewidth=1, color="red", label="Frame 0")
    plt.plot(frames[1], linewidth=1, color="blue", label="Frame +1 (main)")
    plt.plot(frames[2], linewidth=1, color="green", label="Frame -1")
    plt.savefig("count_coding/b_sub_heads/b_sub_frames_sliding_%s_global_head.pdf" % filename)
    plt.close()



gene_set_path = "/home/maria/Documents/pelechanolab/coding_bsub_6633.txt"
gtf_assembly_pickle = "/home/maria/Documents/pelechanolab/gtf_assembled_bsub_6633.sav"
offset = 50
sample_directory = "/home/maria/Documents/pelechanolab/data/samples/ATCC6633/cutadapt/umitools/star/dedup"


def runall_samples(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".bam"):
            yield os.path.join(directory, filename), filename[:5]
    


def gogo():    
    proc = []
    directory = sample_directory
    for bampath, bamname in runall_samples(directory):        
        proc.append(Process(target=plot_results_coding, args=(bampath, bamname)))

    for p in proc:
        p.start()
    for p in proc:
        p.join()

gogo()


