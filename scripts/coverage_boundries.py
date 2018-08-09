#!/usr/bin/python3

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
matplotlib.style.use("ggplot")


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


def fetch_vectors(filename):
    allowed_ids = set(allowed_transcript_ids(global_args.gene_set))
    alignments = BAMGenomeArray(filename, mapping=FivePrimeMapFactory())

    count_vectors = []
    count_vectors_term = []

    for transcript in extend_gtf_frame(global_args.annotation_file):
        if any([transcript.attr.get('Name') in allowed_ids, transcript.get_name() in allowed_ids]):
            try:
                value_array=transcript.get_counts(alignments)
                count_vectors.append(value_array[:global_args.offset*2])
                count_vectors_term.append(value_array[-global_args.offset*2:])
            except ValueError:
                pass

    vector_array=np.vstack(count_vectors)
    vector_array_term=np.vstack(count_vectors_term)

    if global_args.normalize == True:
            vector_array = vector_array[~np.all(vector_array == 0, axis=1)]
            vector_array = vector_array / vector_array.sum(1)[:, np.newaxis]
            vector_array_term = vector_array_term[~np.all(vector_array_term == 0, axis=1)]
            vector_array_term = vector_array_term / vector_array_term.sum(1)[:, np.newaxis]

    metagene=vector_array.sum(axis=0)
    metagene_term=vector_array_term.sum(axis=0)

    return metagene, metagene_term

def scale_labels():
    xlabels = []
    for x in range(-global_args.offset+1, global_args.offset+1):
        if x % 10 == 0:
            xlabels.append(x)
        else:
            xlabels.append(" ")
    return xlabels


def plot_results_start(bampath, bamname):
    metagene=fetch_vectors(bampath)
    labels=scale_labels()
    norm = ""
    if global_args.normalize == True:
        norm = "_norm"

    plt.title(("Peak at: " + str(metagene[0].argmax() - (global_args.offset))))
    plt.grid(True, alpha=0.3)
    plt.step(np.linspace(-global_args.offset, global_args.offset, num=global_args.offset*2),
             metagene[0], linewidth=0.5, label=bamname, fillstyle="full", color="blue", aa=True)
    plt.xticks(np.linspace(-global_args.offset, global_args.offset, num=global_args.offset*2),
               labels, size="xx-small")
    plt.legend()
    plt.savefig(global_args.output_dir + "coverage_start/%s%s.png" % (bamname, norm))
    plt.close()

    plt.title(("Peak at: " + str(metagene[1].argmax() - (global_args.offset))))
    plt.grid(True, alpha=0.3)
    plt.step(np.linspace(-global_args.offset, global_args.offset, num=global_args.offset*2),
             metagene[1], linewidth=0.5, label=bamname, fillstyle="full", color="blue", aa=True)
    plt.xticks(np.linspace(-global_args.offset, global_args.offset, num=global_args.offset*2),
               labels, size="xx-small")
    plt.legend()
    plt.savefig(global_args.output_dir + "coverage_term/%s%s.png" % (bamname, norm))
    plt.close()


def runall_samples(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".bam"):
            yield os.path.join(directory, filename), filename.split(".")[0]


def executable():
    pool=Pool(processes=global_args.cores)
    thread_args=[]
    for bampath, bamname in runall_samples(global_args.sample_dir):
        try:
            os.mkdir(global_args.output_dir + "coverage_start/")
            os.mkdir(global_args.output_dir + "coverage_term/")
        except:
            pass
        thread_args.append((bampath, bamname))

    pool.starmap_async(plot_results_start, thread_args)
    pool.close()
    pool.join()


def executable_2():
    for bampath, bamname in runall_samples(global_args.sample_dir):
        try:
            os.mkdir(global_args.output_dir + "coverage_start/")
            os.mkdir(global_args.output_dir + "coverage_term/")
        except:
            pass

        plot_results_start(bampath, bamname)

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--sample_dir", required=True, help="Path to directory with samples")
    parser.add_argument("--genome_fasta")
    parser.add_argument("--annotation_file", required=True)
    parser.add_argument("--offset", type=int, default=50)
    parser.add_argument("--cores", type=int, default=2)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--gene_set")
    parser.add_argument("--normalize", type=bool, default=False)
    global_args=parser.parse_args()
    if global_args.cores == 1:
        executable_2()
    else:
        executable()
