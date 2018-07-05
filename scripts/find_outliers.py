import numpy as np
from plastid import *
import dill
import matplotlib
import matplotlib.pyplot as plt


def extend_gtf_frame(pickle):
    gtf_coords_file = list(dill.load(open(pickle, "rb")))
    allowed_ids = set(allowed_transcript_ids(gene_set_path))
    for transcript in gtf_coords_file:
        if transcript.get_name() in allowed_ids:
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
#    fasta_dict = match_coords_to_seq(genome_fasta_path)
    alignments = BAMGenomeArray(filename, mapping=FivePrimeMapFactory())
    count_vectors = []
    filler_array = np.zeros(300, dtype=int)
    
    for transcript in extend_gtf_frame(gtf_assembly_pickle):
        if transcript.get_length() > 400:            
            value_array = transcript.get_counts(alignments)
            value_array = np.concatenate((value_array, filler_array))
            count_vectors.append(value_array[50:650])

    vector_array = np.vstack(count_vectors)
    outlier = np.percentile(vector_array[vector_array > 2], 99.0)

    for transcript in extend_gtf_frame(gtf_assembly_pickle):
        value_array = transcript.get_counts(alignments)
        if np.max(value_array) > outlier:
            with open("outliers_t5.txt", "a") as wh:
                wh.write("Outlier" + "\n")
                wh.write(str(transcript.get_name()) + "\n")
                wh.write("Max value: " + str(np.max(value_array)) +"\n")
                argpos = np.argmax(value_array) - 50
                wh.write("Position: " + str(argpos) + "\n")
                for r in range(0, 3):
                    if (argpos + r)%3 == 0:
                        wh.write("Frame: " + str(r) + "\n" + "\n")
                        if r !=0 and argpos > 1:
                            with open ("outlier_ids_t5.txt", "a") as whl:
                                whl.write(str(transcript.get_name()) + "\n")


gene_set_path = "/home/maria/Documents/pelechanolab/scripts/coding_list.txt"
genome_fasta_path = "/home/maria/Documents/pelechanolab/data/R64e/genome.fa"
gtf_assembly_pickle = "gtf_assembled.sav"
bam_file_path = "/home/maria/Documents/pelechanolab/data/samples/20160909YP_08-39666715/cutadapt/umitools/star/dedup/BY4741-t5-1_S8_L002_R1_001.fastqAligned.sortedByCoord.out.bam"
offset = 50
fetch_vectors(bam_file_path)