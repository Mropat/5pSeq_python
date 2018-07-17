import numpy as np
from plastid import *
import dill
import matplotlib
import matplotlib.pyplot as plt
import os


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

    
def allowed_transcript_ids(filename):    
    with open(filename, "r") as gene_list_file:
        return gene_list_file.read().splitlines()


def fetch_vectors(filename):
    allowed_ids = set(allowed_transcript_ids(gene_set_path))
    alignments = BAMGenomeArray(filename, mapping=FivePrimeMapFactory())

    count_vectors = []
    count_vectors_term = []

    for transcript in extend_gtf_frame(gtf_assembly_pickle):
        if transcript.attr.get('Name') in allowed_ids:
        #if transcript.attr.get("gene_biotype") == "protein_coding":        
            try:                          
                value_array = transcript.get_counts(alignments)
                count_vectors.append(value_array[:offset*2])
                count_vectors_term.append(value_array[-offset*2:])
            except ValueError:
                pass

    vector_array = np.vstack(count_vectors)
    vector_array_term = np.vstack(count_vectors_term)
    metagene = vector_array.sum(axis=0)
    metagene_term = vector_array_term.sum(axis=0)

    return metagene, metagene_term

def scale_labels():    
    xlabels = []
    for x in range(-offset, offset):
        if x % 10 == 0:
            xlabels.append(x)
        else:
            xlabels.append(" ")
    return xlabels

def plot_results_start(bampath, bamname): 
    metagene = fetch_vectors(bampath)   
    labels = scale_labels()

    plt.grid(True, alpha=0.3)
    plt.step(np.linspace(-offset, offset, num=offset*2), metagene[0], linewidth=0.5, label=bamname)
    plt.xticks(np.linspace(-offset, offset, num=offset*2), labels, size="xx-small")
    plt.legend()
    plt.savefig("plots/bsub_start/start_codon_coverage_%s.pdf" % bamname)
    plt.close()

    plt.grid(True, alpha=0.3)
    plt.step(np.linspace(-offset, offset, num=offset*2), metagene[1], linewidth=0.5, label=bamname)
    plt.xticks(np.linspace(-offset, offset, num=offset*2), labels, size="xx-small")
    plt.legend()
    plt.savefig("plots/bsub_term/term_codon_coverage_%s.pdf" % bamname)
    plt.close()

gene_set_path = "/home/maria/Documents/pelechanolab/coding_bsub_6633.txt"
gtf_assembly_pickle = "/home/maria/Documents/pelechanolab/gtf_assembled_bsub_6633.sav"
bam_path = "/home/maria/Documents/pelechanolab/data/samples/ATCC6633/cutadapt/umitools/star/dedup"
offset = 50

def runall_samples(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".bam"):
            yield os.path.join(directory, filename), filename[:5]

for bampath, bamname in runall_samples(bam_path):
    plot_results_start(bampath, bamname)


