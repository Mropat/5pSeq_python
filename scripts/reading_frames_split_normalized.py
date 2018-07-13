# Script for analyzing and plotting the ribosome protection 

import numpy as np
from plastid import *
import dill
import matplotlib
import matplotlib.pyplot as plt
from time import time
from multiprocessing import Process
import os

def timeit(f):
    '''Gives feedback on how long plotting takes (optional) '''
    def timer(*args, **kwargs):
        st = int(time())
        print("Start!")
        f(*args, **kwargs)
        et = int(time())
        print('%s took %i seconds to complete!' % (f.__name__, et - st))
    return timer


#Reads a file with line separated transcript id. Any group of transcripts can be analyzed
#In this instance allows only protein coding transcripts in the list to be analyzed

def allowed_transcript_ids(filename):
    with open(filename, "r") as gene_list_file:
        return gene_list_file.read().splitlines()


#Generates vectors of 5P coverages using mapping coordinates

def fetch_vectors(filenames):
    ''' 
    gtf_coords_file contains a list of transcript boundries and intron coordinates
    Connected to transcript ID, genomic fasta coordinates, splicing information if needed
    and other information extracted from gtf annotation. 
    This file takes some time to generate, so it is saved with dill(same as pickle but handles complicated classes)
    and later loarded and reused for all the samples mapped with this gtf
    '''
    gtf_coords_file = list(dill.load(open(gtf_assembly_pickle, "rb")))
    allowed_ids = set(allowed_transcript_ids(gene_set_path))

#This part extends the boundries of the transcript by chosen amount of bases. 
#This is the easiest way to do this and retain all transcript information. 
#If the boundries are expanded beyond what is indexed in BAM file, 
#script may produce error when trying to retrieve coverage at later stage

    for transcript in gtf_coords_file:
        span = transcript.spanning_segment
        new_region = GenomicSegment(
            span.chrom, span.start - offset, span.start, span.strand)
        new_region_2 = GenomicSegment(
            span.chrom, span.end, span.end + offset, span.strand)
        transcript.add_segments(new_region, new_region_2)

#This reads in the mapped, indexed bam files.
#For this to work, bai file should be in the same folder
# BAMGenomeArray can also read list of multiple bam files and combine them, 
# In case we want to analyze multiple repeat experiments together  
# Mapping rule can be changed (see documentation), but we retrieve only 5P ends of 
# each read currently
         
    alignments = BAMGenomeArray(filenames, mapping=FivePrimeMapFactory())
    print("Genomes loaded!")

# Fetches vectors corresponding to transcript coordinates from BAM files
# transcripts store information about introns, id, length, genomic coordinates
# and more (see documentation). 
# We can allow only transcripts in our allowed ID list to retrieve vectors
# Or define lengths of transcripts we want to view
# Or even retrieve vectors surrounding a specific codon
# If vectors are going to be analyzed for frame protection,
# It is necessary to make sure the vector length is divisible by 3

    count_vectors = []
    for transcript in gtf_coords_file:
        if transcript.attr.get('Name') in allowed_ids:
            value_array = transcript.get_counts(alignments)
            if np.sum(value_array[-1800:]) > 1:
#            if np.sum(value_array[:1800]) > 1:
#                count_vectors.append(np.concatenate((value_array, np.zeros(4000, dtype=int)))[:1800])
                count_vectors.append(np.concatenate((np.zeros(4000, dtype=int), value_array))[-1800:])


    vector_array = np.vstack(count_vectors)
    print("Vectors retrieved!")

#This normalizes vectors by read (optional) 
    vector_normsum = np.sum(vector_array, axis=1)
    vector_array_normalized = vector_array / vector_normsum[:, np.newaxis]


#This creates a metagene from all retrieved vectors
    metagene = vector_array_normalized.sum(axis=0)

# Vectors are split into 3 frames for ribosome protection analysis
    metagene_stack = np.reshape(metagene, (-1, 3))
    metagene_stack = np.hsplit(metagene_stack, 3)

#This transforms each frame vector with sliding window of 5 for plotting for smoother result
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
    plt.savefig("normalized_coding/b_sub_tails/b_sub_frames_sliding_%s_global_tail.pdf" % filename)
    plt.close()


# Paths to genome coordinate file, allowed id set
gene_set_path = "/home/maria/Documents/pelechanolab/coding_bsub_6633.txt"
gtf_assembly_pickle = "/home/maria/Documents/pelechanolab/gtf_assembled_bsub_6633.sav"
offset = 0

# fetches sample names from directory in which we analyze for plot labeling
def runall_samples(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".bam"):
            yield os.path.join(directory, filename), filename[:5]
    

@timeit
# (Optional) when many samples are analyzed, this will arrange multithreading. 

def gogo():    
    proc = []
    directory = "/home/maria/Documents/pelechanolab/data/samples/ATCC6633/cutadapt/umitools/star/dedup"
    for bampath, bamname in runall_samples(directory):        
        proc.append(Process(target=plot_results_coding, args=(bampath, bamname)))

    for p in proc:
        p.start()
    for p in proc:
        p.join()

# gogo()

print(help(plot_results_coding))
