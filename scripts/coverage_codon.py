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

def fetch_vectors(filename, filename_cont):
    fasta_dict = match_coords_to_seq(genome_fasta_path)
    allowed_ids = set(allowed_transcript_ids(gene_set_path))

    alignments = BAMGenomeArray(filename, mapping=FivePrimeMapFactory())
    alignments_cont = BAMGenomeArray(filename_cont, mapping=FivePrimeMapFactory())

    count_vectors = []
    count_vectors_cont = []


    for transcript in extend_gtf_frame(gtf_assembly_pickle):
        if transcript.get_name() in allowed_ids:
            readseq = transcript.get_sequence(fasta_dict)
            for ind in range(50, len(readseq)-50, 3):
                codon = readseq[ind] + readseq[ind+1] + readseq[ind+2]
                if codon == "CCG":
                    count_vectors.append(transcript.get_counts(alignments)[ind-50:ind+50])
                    count_vectors_cont.append(transcript.get_counts(alignments_cont)[ind-50:ind+50])

    vector_array = np.vstack(count_vectors)
    vector_array_cont = np.vstack(count_vectors_cont)


    metagene = vector_array.sum(axis=0)
    metagene_cont = vector_array_cont.sum(axis=0)

    return metagene, metagene_cont

def scale_labels():    
    xlabels = []
    for x in range(-49, 49):
        if x % 10 == 0:
            xlabels.append(x)
        else:
            xlabels.append(" ")
    return xlabels

def plot_results_start(metagene, metagene_cont):    
    labels = scale_labels()
    plt.grid(True, alpha=0.3)
    plt.step(np.linspace(-49, 51, num=100), metagene, linewidth=0.5, color="red")
    plt.step(np.linspace(-49, 51, num=100), metagene_cont, linewidth=0.5, color="blue")
    plt.xticks(np.linspace(-49, 51, num=100), labels, size="xx-small")
    plt.savefig("codon_coverage.pdf")
    plt.close()
    print("Peak at: " + str(metagene.argmax() - (offset)))


offset = 50
genome_fasta_path = "/home/maria/Documents/pelechanolab/data/R64e/genome.fa"
gtf_assembly_pickle = "gtf_assembled.sav"
bam_file_path = "/home/maria/Documents/pelechanolab/data/samples/2/post-proc/by14_clu_dedup_nomm.bam"
bam_file_cont_path = "/home/maria/Documents/pelechanolab/data/samples/20160909YP_08-39666715/cutadapt/umitools/star/dedup/BY4741-t5-1_S8_L002_R1_001.fastqAligned.sortedByCoord.out.bam"
gene_set_path = "/home/maria/Documents/pelechanolab/scripts/coding_list.txt"


metagene, metagene_cont = fetch_vectors(bam_file_path, bam_file_cont_path)
plot_results_start(metagene, metagene_cont)