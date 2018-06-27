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

def fetch_vectors(filename, filename_dedup):
    fasta_dict = match_coords_to_seq(genome_fasta_path)
    allowed_ids = set(allowed_transcript_ids(gene_set_path))

    alignments = BAMGenomeArray(filename, mapping=FivePrimeMapFactory())
    alignments_dedup = BAMGenomeArray(filename_dedup, mapping=FivePrimeMapFactory())

    count_vectors = []
    count_vectors_dedup = []
    count_vectors_term = []
    count_vectors_dedup_term = []

    for transcript in extend_gtf_frame(gtf_assembly_pickle):
        if transcript.get_name() in allowed_ids:
            value_array = transcript.get_counts(alignments)
            count_vectors.append(value_array[1:100])
            count_vectors_term.append(value_array[-99:])
            value_array = transcript.get_counts(alignments_dedup)
            count_vectors_dedup.append(value_array[1:100])
            count_vectors_dedup_term.append(value_array[-99:])

    vector_array = np.vstack(count_vectors)
    vector_array_dedup = np.vstack(count_vectors_dedup)
    vector_array_term = np.vstack(count_vectors_term)
    vector_array_dedup_term = np.vstack(count_vectors_dedup_term)

    metagene = vector_array.sum(axis=0)
    metagene_dedup = vector_array_dedup.sum(axis=0)
    metagene_term = vector_array_term.sum(axis=0)
    metagene_dedup_term = vector_array_dedup_term.sum(axis=0)

    return metagene, metagene_dedup, metagene_term, metagene_dedup_term

def scale_labels():    
    xlabels = []
    for x in range(-48, 49):
        if x % 10 == 0:
            xlabels.append(x)
        else:
            xlabels.append(" ")
    return xlabels

def plot_results_start(metagene, metagene_dedup):    
    labels = scale_labels()
    plt.grid(True, alpha=0.3)
    plt.step(np.linspace(-49, 51, num=99), metagene, linewidth=0.5, color="red")
    plt.step(np.linspace(-49, 51, num=99), metagene_dedup, linewidth=0.5, color="blue")
    plt.xticks(np.linspace(-49, 51, num=99), labels, size="xx-small")
    plt.savefig("rdn_dedup_head.pdf")
    plt.close()
    print(metagene.argmax() - (offset-1))

def plot_results_stop(metagene_term, metagene_dedup_term):
    labels = scale_labels()
    plt.grid(True, alpha=0.3)
    plt.step(np.linspace(-49, 51, num=99), metagene_term, linewidth=0.5, color="red")
    plt.step(np.linspace(-49, 51, num=99), metagene_dedup_term, linewidth=0.5, color="blue")
    plt.xticks(np.linspace(-49, 51, num=99), labels, size="xx-small")
    plt.savefig("rdn_dedup_tail.pdf")
    plt.close()
    print(metagene_term.argmax() - (offset-1))


gene_set_path = "/home/maria/Documents/pelechanolab/scripts/coding_list.txt"
genome_fasta_path = "/home/maria/Documents/pelechanolab/data/R64e/genome.fa"
gtf_assembly_pickle = "gtf_assembled.sav"
bam_file_path = "/home/maria/Documents/pelechanolab/data/samples/20160909YP_08-39666715/cutadapt/umitools/star/BY4741-t5-1_S8_L002_R1_001.fastqAligned.sortedByCoord.out.bam"
bam_file_dedup_path = "/home/maria/Documents/pelechanolab/data/samples/20160909YP_08-39666715/cutadapt/umitools/star/dedup/BY4741-t5-1_S8_L002_R1_001.fastqAligned.sortedByCoord.out.bam"
offset = 50


metagene, metagene_dedup, metagene_term, metagene_dedup_term = fetch_vectors(bam_file_path, bam_file_dedup_path)
plot_results_stop(metagene, metagene_dedup)
plot_results_start(metagene_term, metagene_dedup_term)