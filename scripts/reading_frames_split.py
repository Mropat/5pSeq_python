import numpy as np
from plastid import *
import dill
import matplotlib
import matplotlib.pyplot as plt


def extend_gtf_frame(pickle):
    gtf_coords_file = list(dill.load(open(pickle, "rb")))
    for transcript in gtf_coords_file:
        span = transcript.spanning_segment
        new_region = GenomicSegment(
            span.chrom, span.start - offset, span.start, span.strand)
        new_region_2 = GenomicSegment(
            span.chrom, span.end, span.end + offset, span.strand)
        transcript.add_segments(new_region, new_region_2)
    return gtf_coords_file


def match_coords_to_seq(filename):
    fasta_dict = {}
    with open(filename, "r") as genomeseq:
        chrom = ""
        genseq = ""
        genomelines = genomeseq.readlines()
        for line in genomelines:
            if line.startswith(">"):
                fasta_dict[chrom] = genseq
                chrom = line.strip()[1:]
                genseq = ""
            else:
                genseq = genseq + line.strip()
        fasta_dict[chrom] = genseq
    return fasta_dict


def allowed_transcript_ids(filename):
    with open(filename, "r") as gene_list_file:
        return gene_list_file.read().splitlines()


def fetch_vectors(filename):
    fasta_dict = match_coords_to_seq(genome_fasta_path)
    allowed_ids = set(allowed_transcript_ids(gene_set_path))
    alignments = BAMGenomeArray(filename, mapping=FivePrimeMapFactory())

    count_vectors = []
    filler_array = np.zeros(300, dtype=int)
    for transcript in extend_gtf_frame(gtf_assembly_pickle):
        transcript_length = transcript.get_length()
        if transcript.get_name() in allowed_ids:
            if transcript_length in range(400, 800):
                value_array = transcript.get_counts(alignments)
                value_array = np.concatenate((value_array, filler_array))
                count_vectors.append(value_array[:600])

    vector_array = np.vstack(count_vectors)
    metagene = vector_array.sum(axis=0)
    metagene_stack = np.reshape(metagene, (200, 3))
    metagene_stack = np.hsplit(metagene_stack, 3)

    frames = []
    for arr in metagene_stack:
        frame_vector = np.reshape(arr.T, (1, -1))
        frame_vector = np.concatenate(
            (np.zeros(2), frame_vector[0], np.zeros(2)))

        window_iter = (np.sum(frame_vector[i-2:i+2])
                       for i in range(len(frame_vector[2:-3])))
        sliding_vector = np.fromiter(window_iter, dtype=int)
        frames.append(sliding_vector)

    return frames


def scale_labels():
    xlabels = []
    for x in range(-48, 49):
        if x % 10 == 0:
            xlabels.append(x)
        else:
            xlabels.append(" ")
    return xlabels


def plot_results_coding(frames):
    labels = scale_labels()

    frame_0, frame_1, frame_2 = frames[0], frames[1], frames[2]

    plt.grid(True, alpha=0.3)
    plt.plot(frame_0, linewidth=0.5, color="red")
    plt.plot(frame_1, linewidth=0.5, color="blue")
    plt.plot(frame_2, linewidth=0.5, color="green")

    plt.savefig("frames_sliding.pdf")
    plt.close()


gene_set_path = "/home/maria/Documents/pelechanolab/scripts/coding_list.txt"
genome_fasta_path = "/home/maria/Documents/pelechanolab/data/R64e/genome.fa"
gtf_assembly_pickle = "gtf_assembled.sav"
bam_file_path = "/home/maria/Documents/pelechanolab/data/samples/20160909YP_08-39666715/cutadapt/umitools/star/BY4741-t5-1_S8_L002_R1_001.fastqAligned.sortedByCoord.out.bam"
offset = 50


frames = fetch_vectors(bam_file_path)
plot_results_coding(frames)
