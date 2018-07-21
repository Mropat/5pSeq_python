import numpy as np
from plastid import GTF2_TranscriptAssembler, GFF3_TranscriptAssembler, Transcript, GenomicSegment, BAMGenomeArray, FivePrimeMapFactory
import dill
import matplotlib
import matplotlib.pyplot as plt
from multiprocessing import Process, Pool
import os
import argparse
import sys
import matplotlib.style
matplotlib.style.use("ggplot")
import pandas as pd
from plastid import plotting
from plastid.plotting.plots import *


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


def allowed_transcript_ids(filename):
    try:
        with open(filename, "r") as gene_list_file:
            return gene_list_file.read().splitlines()
    except:
        return list()


def fetch_vectors(filename):
    pickle_path = global_args.output_dir + \
        global_args.annotation_file[:-4].split("/")[-1] + ".sav"
    if os.path.isfile(pickle_path) == False:
        create_assembly_dill(global_args.annotation_file)

    gtf_coords_file = list(dill.load(open(pickle_path, "rb")))
    allowed_ids = set(allowed_transcript_ids(global_args.gene_set))
    alignments = BAMGenomeArray(filename, mapping=FivePrimeMapFactory())
    print("Genomes loaded for %s " % filename)

    vector_array = []

    for transcript in gtf_coords_file:
        if any([transcript.attr.get('Name') in allowed_ids, transcript.attr.get("type") == "mRNA", transcript.attr.get("gene_biotype") == "protein_coding"]):
            readvec = transcript.get_counts(alignments)
            if np.sum(readvec) > 100 and len(readvec)%3 ==0:
                readvec = np.reshape(readvec, (-1, 3))
                vector_array.append(np.sum(readvec, axis=0))

    vector_array = np.vstack(vector_array)
    vector_array = vector_array / vector_array.sum(1)[:, numpy.newaxis]
    return vector_array


def plot_results_start(bampath, bamname):
    metagene = fetch_vectors(bampath)

    fig, ax = triangle_plot(metagene, grid=[
                            0.5, 0.75], marker=".", linewidth=0.1, alpha=0.2, vertex_labels=["A", "B", "C"])

    plt.title("frame preference in %s" % bamname)

    plt.savefig(global_args.output_dir + "trigplot/%s.pdf" % bamname)
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
            os.mkdir(global_args.output_dir + "trigplot/")

        except:
            pass
        thread_args.append((bampath, bamname))

    pool.starmap_async(plot_results_start, thread_args)
    pool.close()
    pool.join()


def executable_2():
    for bampath, bamname in runall_samples(global_args.sample_dir):
        try:
            os.mkdir(global_args.output_dir + "trigplot/")
        except:
            pass

        plot_results_start(bampath, bamname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_dir", required=True)
    parser.add_argument("--annotation_file", required=True)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--gene_set")
    global_args = parser.parse_args()
    if global_args.cores == 1:
        executable_2()
    else:
        executable()
