# Script for plotting metagene coverage around a specific codon of interest, like Proline CCG
# Or a specific position of interest, such as start/stop codon

import numpy as np
from plastid import *
import dill
import matplotlib
import matplotlib.pyplot as plt
from multiprocessing import Process, Pool
import os
import argparse
import sys


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
                chrom = line.strip().split(" ")[0][1:]
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
    codon_dict = {}
    print("Genomes loaded")

    for transcript in extend_gtf_frame(gtf_assembly_pickle):
        if transcript.attr.get('Name') in allowed_ids:
            try:
                # Necessary! Exception if out of bounds
                readvec = transcript.get_counts(alignments)
                readseq = transcript.get_sequence(fasta_dict)

                for ind in range(offset+3, len(readseq)-offset, 3):
                    codon = readseq[ind] + readseq[ind+1] + readseq[ind+2]
                    if codon not in codon_dict:
                        codon_dict[codon] = np.atleast_2d(
                            readvec[ind-offset:ind+offset])
                    else:
                        codon_dict[codon] = np.concatenate((codon_dict[codon], np.atleast_2d(
                            readvec[ind-offset:ind+offset])), axis=0)
            except ValueError:
                pass

    for key in codon_dict:
        codon_dict[key] = codon_dict[key].sum(axis=0)
    print("file processed")

    return codon_dict


def scale_labels():
    xlabels = []
    for x in range(-offset, offset):
        if x % 10 == 0:
            xlabels.append(x)
        else:
            xlabels.append(" ")
    return xlabels


def plot_results_start(bam_file_path, bamname):
    codon_dict = fetch_vectors(bam_file_path)
    labels = scale_labels()

    for key in codon_dict:

        plt.title(("Peak at: " + str(codon_dict[key].argmax() - (offset))))
        plt.grid(True, alpha=0.2)
        plt.step(np.linspace(-offset, offset, num=offset*2),
                 codon_dict[key], linewidth=0.5, color="red")
        plt.xticks(np.linspace(-offset, offset, num=offset*2),
                   labels, size="xx-small")
        plt.savefig("plots/l_plan_codons/%s/codon_coverage_%s.pdf" %
                    (bamname, key))
        plt.close()

        print("Peak at: " + str(codon_dict[key].argmax() - (offset)))

    for key in aa_dict:
        collapsed_codons = []

        for x in aa_dict[key]:
            collapsed_codons.append(np.atleast_2d(codon_dict[x]))

        collapsed_aa = np.vstack(collapsed_codons).sum(axis=0)

        plt.title(("Peak at: " + str(collapsed_aa.argmax() - (offset))))
        plt.grid(True, alpha=0.2)
        plt.step(np.linspace(-offset, offset, num=offset*2),
                 collapsed_aa, linewidth=0.5, color="red")
        plt.xticks(np.linspace(-offset, offset, num=offset*2),
                   labels, size="xx-small")
        plt.savefig("plots/l_plan_aa/%s/codon_coverage_%s.pdf" %
                    (bamname, key))
        plt.close()


# Paths and args to files used
offset = 50
cores = 4
sample_directory = "/home/maria/Documents/pelechanolab/data/samples/ATCC8014/umitools/star/dedup"
genome_fasta_path = "/home/maria/Documents/pelechanolab/data/l_plan/GCA_002631775.1_ASM263177v1_genomic.fna"
gtf_assembly_pickle = "/home/maria/Documents/pelechanolab/gtf_assembled_l_plan8014.sav"
gene_set_path = "/home/maria/Documents/pelechanolab/coding_lplan8014.txt"


aa_dict = {
    "ALA": ["GCT", "GCC", "GCA", "GCG"],
    "ARG": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "ASN": ["AAT", "AAC"],
    "ASP": ["GAT", "GAC"],
    "CYS": ["TGT", "TGC"],
    "GLN": ["TGT", "TGC"],
    "GLU": ["TGT", "TGC"],
    "GLY": ["GGT", "GGC", "GGA", "GGG"],
    "HIS": ["CAT", "CAC"],
    "ILE": ["ATT", "ATC", "ATA"],
    "LEU": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "LYS": ["AAA", "AAG"],
    "PHE": ["TTT", "TTC"],
    "PRO": ["CCT", "CCC", "CCA", "CCG"],
    "SER": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "THR": ["ACT", "ACC", "ACA", "ACG"],
    "TRP": ["TGG"],
    "TYR": ["TAT", "TAC"],
    "VAL": ["GTT", "GTC", "GTA", "GTG"],
    "TERM": ["TAA", "TGA", "TAG"],
    "MET": ["ATG"]}


def runall_samples(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".bam"):
            yield os.path.join(directory, filename), filename[:5]


def gogo():

    pool = Pool(processes=cores)
    args = []
    for bampath, bamname in runall_samples(sample_directory):
        try:
            os.mkdir("plots/l_plan_aa/" + bamname)
            os.mkdir("plots/l_plan_codons/" + bamname)
        except:
            pass
        args.append((bampath, bamname))

    pool.starmap_async(plot_results_start, args)
    pool.close()
    pool.join()


gogo()



