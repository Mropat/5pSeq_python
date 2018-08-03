# Script for plotting the ribosome protection frames

import numpy as np
from plastid import GTF2_TranscriptAssembler, GFF3_TranscriptAssembler, Transcript, GenomicSegment, BAMGenomeArray, FivePrimeMapFactory
import dill
import matplotlib
import matplotlib.pyplot as plt
import os
import argparse
import sys
import multiprocessing
from multiprocessing import Process, Pool


def create_assembly_dill(annotation_file):
    if annotation_file.endswith("gtf"):
        gtf_file = list(GTF2_TranscriptAssembler(
            annotation_file, return_type=Transcript))
        dill.dump(gtf_file, open(global_args.output_dir +
                                 annotation_file[:-4].split("/")[-1] + ".sav", "wb"))
    elif annotation_file.endswith("gff"):
        gff_file = list(GFF3_TranscriptAssembler(
            annotation_file, return_type=Transcript))
        dill.dump(gff_file, open(global_args.output_dir +
                                 annotation_file[:-4].split("/")[-1] + ".sav", "wb"))


def extend_gtf_frame(pickle):

    pickle_path = global_args.output_dir + \
        global_args.annotation_file[:-4].split("/")[-1] + ".sav"

    if os.path.isfile(pickle_path) == False:
        create_assembly_dill(global_args.annotation_file)
    gtf_coords_file = list(dill.load(open(pickle_path, "rb")))

    for transcript in gtf_coords_file:
        span = transcript.spanning_segment
        new_region = GenomicSegment(
            span.chrom, span.start - global_args.offset, span.start, span.strand)
        new_region_2 = GenomicSegment(
            span.chrom, span.end, span.end + global_args.offset, span.strand)
        transcript.add_segments(new_region, new_region_2)
        yield transcript


def allowed_transcript_ids(filename):
    try:
        with open(filename, "r") as gene_list_file:
            return gene_list_file.read().splitlines()
    except:
        return list()


def fetch_vectors(filenames):

    allowed_ids = set(allowed_transcript_ids(global_args.gene_set))
    alignments = BAMGenomeArray(filenames, mapping=FivePrimeMapFactory())
    print("Genomes loaded!")

    except_count = 0
    count_vectors_start = []
    count_vectors_term = []

    for transcript in extend_gtf_frame(global_args.annotation_file):
        if any([transcript.attr.get('Name') in allowed_ids, transcript.get_name() in allowed_ids]):
            try:
                value_array = transcript.get_counts(alignments)
                if global_args.shortest<(len(value_array)-global_args.offset*2)<global_args.longest:

                    if np.sum(value_array[(-global_args.longest-global_args.offset):]) > 1:
                        count_vectors_term.append(np.concatenate(
                            (np.zeros(global_args.longest, dtype=int), value_array))[-global_args.longest-global_args.offset:])

                    if np.sum(value_array[:global_args.longest+global_args.offset]) > 1:
                        count_vectors_start.append(np.concatenate(
                            (value_array, np.zeros(global_args.longest, dtype=int)))[:global_args.longest+global_args.offset])
            except ValueError:
                except_count += 1

    vector_array_start = np.vstack(count_vectors_start)
    vector_array_term = np.vstack(count_vectors_term)

    print("Vectors retrieved!")
    print("Removed %i transcripts!" % except_count)

    if global_args.normalize == True:
        vector_normsum_start = np.sum(vector_array_start, axis=1)
        vector_array_start = vector_array_start / \
            vector_normsum_start[:, np.newaxis]

        vector_normsum_term = np.sum(vector_array_term, axis=1)
        vector_array_term = vector_array_term / \
            vector_normsum_term[:, np.newaxis]

    metagene_start = vector_array_start.sum(axis=0)
    metagene_stack_start = np.reshape(metagene_start, (-1, 3))
    metagene_stack_start = np.hsplit(metagene_stack_start, 3)

    metagene_term = vector_array_term.sum(axis=0)
    metagene_stack_term = np.reshape(metagene_term, (-1, 3))
    metagene_stack_term = np.hsplit(metagene_stack_term, 3)

    frames_start = []
    for arr in metagene_stack_start:
        frame_vector = np.concatenate((np.zeros(2), arr.T[0], np.zeros(2)))
        window_iter = (np.sum(frame_vector[i:i+5])
                       for i in range(len(frame_vector[2:-3])))
        frames_start.append(np.fromiter(window_iter, dtype=float))

    frames_term = []
    for arr in metagene_stack_term:
        frame_vector = np.concatenate((np.zeros(2), arr.T[0], np.zeros(2)))
        window_iter = (np.sum(frame_vector[i:i+5])
                       for i in range(len(frame_vector[2:-3])))
        frames_term.append(np.fromiter(window_iter, dtype=float))

    print("Frames split")
    return frames_start, frames_term


def plot_results_start(bampath, bamname):
    frames = fetch_vectors(bampath)
    norm = ""
    if global_args.normalize == True:
        norm = "_norm"

    plt.grid(True, alpha=0.3)
    #plt.xticks(range(int((global_args.offset/3)-2 + global_args.longest/3, 10)), labels=([x for x in range((-global_args.offset/3)-2, global_args.longest/3, 100)]))
    plt.plot(frames[0][0][2:], linewidth=1, color="red", label="Frame 0")
    plt.plot(frames[0][1][2:], linewidth=1, color="blue", label="Frame +1 (main)")
    plt.plot(frames[0][2][2:], linewidth=1, color="green", label="Frame -1")
    plt.savefig(global_args.output_dir +
                "frames_coverage_start/ribosome_protection_%s%s.pdf" % (bamname, norm))
    plt.close()

    plt.grid(True, alpha=0.3)
    #plt.xticks(range(int((global_args.offset/3)-2 + global_args.longest/3, 10)), labels=([x for x in range((-global_args.longest/3, global_args.offset/3)-2,  100)]))
    plt.plot(frames[1][0][:-2], linewidth=1, color="red", label="Frame 0")
    plt.plot(frames[1][1][:-2], linewidth=1, color="blue", label="Frame +1 (main)")
    plt.plot(frames[1][2][:-2], linewidth=1, color="green", label="Frame -1")
    plt.savefig(global_args.output_dir +
                "frames_coverage_term/ribosome_protection_%s%s.pdf" % (bamname, norm))
    plt.close()


def runall_samples(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".bam"):
            yield os.path.join(directory, filename), filename.split(".")[0]


def executable():
    pool = Pool(processes=global_args.cores)
    thread_args = []
    for bampath, bamname in runall_samples(global_args.sample_dir):
        try:
            os.mkdir(global_args.output_dir + "frames_coverage_start/")
            os.mkdir(global_args.output_dir + "frames_coverage_term/")
        except:
            pass
        thread_args.append((bampath, bamname))

    pool.starmap_async(plot_results_start, thread_args)
    pool.close()
    pool.join()


def executable_2():
    for bampath, bamname in runall_samples(global_args.sample_dir):
        try:
            os.mkdir(global_args.output_dir + "frames_coverage_start/")
            os.mkdir(global_args.output_dir + "frames_coverage_term/")
        except:
            pass

        plot_results_start(bampath, bamname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_dir", required=True)
    parser.add_argument("--annotation_file", required=True)
    parser.add_argument("--offset", type=int, default=36)
    parser.add_argument("--cores", type=int, default=2)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--gene_set")
    parser.add_argument("--normalize", type=bool, default=False)
    parser.add_argument("--shortest", type=int, default=0)
    parser.add_argument("--longest", type=int, default=1800)
    parser.add_argument("--genome_fasta")

    global_args = parser.parse_args()

    if global_args.cores == 1:
        executable_2()
    else:
        executable()
