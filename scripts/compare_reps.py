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
        yield transcript

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

def fetch_vectors(filename, filename_rep1, filename_rep2):
    fasta_dict = match_coords_to_seq(genome_fasta_path)
    allowed_ids = set(allowed_transcript_ids(gene_set_path))

    alignments = BAMGenomeArray(filename, mapping=FivePrimeMapFactory())
    print("First genome done!")
    alignments_rep1 = BAMGenomeArray(filename_rep1, mapping=FivePrimeMapFactory())
    print("Second genome done!")
    alignments_rep2 = BAMGenomeArray(filename_rep2, mapping=FivePrimeMapFactory())
    print("Third genome done!")

    count_vectors = []
    count_vectors_rep1 = []
    count_vectors_rep2 = []

    count_vectors_term = []
    count_vectors_rep1_term = []
    count_vectors_rep2_term = []

    for transcript in extend_gtf_frame(gtf_assembly_pickle):
        if transcript.get_name() in allowed_ids:
            
            value_array = transcript.get_counts(alignments)
            count_vectors.append(value_array[1:100])
            count_vectors_term.append(value_array[-99:])

            value_array = transcript.get_counts(alignments_rep1)
            count_vectors_rep1.append(value_array[1:100])
            count_vectors_rep1_term.append(value_array[-99:])

            value_array = transcript.get_counts(alignments_rep2)
            count_vectors_rep2.append(value_array[1:100])
            count_vectors_rep2_term.append(value_array[-99:])

    vector_array = np.vstack(count_vectors)
    vector_array_rep1 = np.vstack(count_vectors_rep1)
    vector_array_rep2 = np.vstack(count_vectors_rep2)

    vector_array_term = np.vstack(count_vectors_term)
    vector_array_rep1_term = np.vstack(count_vectors_rep1_term)
    vector_array_rep2_term = np.vstack(count_vectors_rep2_term)

    metagene = vector_array.sum(axis=0)
    metagene_rep1 = vector_array_rep1.sum(axis=0)
    metagene_rep2 = vector_array_rep2.sum(axis=0)

    metagene_term = vector_array_term.sum(axis=0)
    metagene_rep1_term = vector_array_rep1_term.sum(axis=0)
    metagene_rep2_term = vector_array_rep2_term.sum(axis=0)

    print("Metagenes done!")

    return metagene, metagene_rep1, metagene_rep2, metagene_term, metagene_rep1_term, metagene_rep2_term

def scale_labels():    
    xlabels = []
    for x in range(-48, 49):
        if x % 10 == 0:
            xlabels.append(x)
        else:
            xlabels.append(" ")
    return xlabels

def plot_results_start(metagene, metagene_rep1, metagene_rep2):    
    labels = scale_labels()
    plt.grid(True, alpha=0.3)
    plt.step(np.linspace(-49, 51, num=99), metagene, linewidth=0.5, color="red")
    plt.step(np.linspace(-49, 51, num=99), metagene_rep1, linewidth=0.5, color="green")
    plt.step(np.linspace(-49, 51, num=99), metagene_rep2, linewidth=0.5, color="blue")
    plt.xticks(np.linspace(-49, 51, num=99), labels, size="xx-small")
    plt.savefig("rdn_rep_t30_head.pdf")
    plt.close()


def plot_results_stop(metagene_term, metagene_rep1_term, metagene_rep2_term):
    labels = scale_labels()
    plt.grid(True, alpha=0.3)
    plt.step(np.linspace(-49, 51, num=99), metagene_term, linewidth=0.5, color="red")
    plt.step(np.linspace(-49, 51, num=99), metagene_rep1_term, linewidth=0.5, color="green")
    plt.step(np.linspace(-49, 51, num=99), metagene_rep2_term, linewidth=0.5, color="blue")
    plt.xticks(np.linspace(-49, 51, num=99), labels, size="xx-small")
    plt.savefig("rdn_rep_t30_tail.pdf")
    plt.close()


gene_set_path = "/home/maria/Documents/pelechanolab/scripts/coding_list.txt"
genome_fasta_path = "/home/maria/Documents/pelechanolab/data/R64e/genome.fa"
gtf_assembly_pickle = "gtf_assembled.sav"
bam_file_path = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment/cutadapt/umitools/star/dedup/BY4741-t30-2_S26.fastqAligned.sortedByCoord.out.bam"
bam_file_rep1_path = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment_2/cutadapt/umitools/star/dedup/BY4741-t30-1_S10.fastqAligned.sortedByCoord.out.bam"
bam_file_rep2_path = "/home/maria/Documents/pelechanolab/data/samples/WT_ox_stress_experiment_3/cutadapt/umitools/star/dedup/BY4741-t30-3_S18.fastqAligned.sortedByCoord.out.bam"
offset = 50


metagene, metagene_rep1, metagene_rep2, metagene_term, metagene_rep1_term, metagene_rep2_term = fetch_vectors(bam_file_path, bam_file_rep1_path, bam_file_rep2_path)
plot_results_stop(metagene, metagene_rep1, metagene_rep2)
plot_results_start(metagene_term, metagene_rep1_term, metagene_rep2_term)