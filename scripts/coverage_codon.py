# Script for plotting metagene coverage around a specific codon of interest, like Proline CCG

import numpy as np
from plastid import GTF2_TranscriptAssembler, GFF3_TranscriptAssembler, Transcript, GenomicSegment, BAMGenomeArray, FivePrimeMapFactory
import dill
import matplotlib
import matplotlib.pyplot as plt
from multiprocessing import Process, Pool
import os
import argparse
import sys


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
    try:
        with open(filename, "r") as gene_list_file:
            return gene_list_file.read().splitlines()
    except:
        return list()


def fetch_vectors(filename):
    fasta_dict = match_coords_to_seq(global_args.genome_fasta)
    allowed_ids = set(allowed_transcript_ids(global_args.gene_set))
    alignments = BAMGenomeArray(filename, mapping=FivePrimeMapFactory())
    codon_dict = {}
    out_range_count = 0

    print("Genomes loaded for %s " % filename)

    for transcript in extend_gtf_frame(global_args.annotation_file):
        if any([transcript.attr.get('Name') in allowed_ids, transcript.get_name() in allowed_ids]):
            try:
                # Necessary! Exception if out of bounds
                readvec = transcript.get_counts(alignments)
                readseq = transcript.get_sequence(fasta_dict)

                for ind in range(global_args.offset+3, len(readseq)-global_args.offset, 3):
                    codon = readseq[ind:ind+3]
                    if codon not in codon_dict:
                        codon_dict[codon] = np.atleast_2d(
                            readvec[ind-global_args.offset:ind+global_args.offset])
                    else:
                        codon_dict[codon] = np.concatenate((codon_dict[codon], np.atleast_2d(
                            readvec[ind-global_args.offset:ind+global_args.offset])), axis=0)
            except ValueError:
                out_range_count += 1
    print ("Transctipts out of bounds for %s: %i" % (filename, out_range_count))
    for key in codon_dict:
        if global_args.normalize == True:
            codon_dict[key] = codon_dict[key][~np.all(codon_dict[key] == 0, axis=1)]
            codon_dict[key] = codon_dict[key] / codon_dict[key].sum(1)[:, np.newaxis]
        codon_dict[key] = codon_dict[key].sum(axis=0)
    print("Codons gathered for %s" % filename)

    return codon_dict


def scale_labels():
    xlabels = []
    for x in range(-global_args.offset+1, global_args.offset+1):
        if x % 10 == 0:
            xlabels.append(x)
        else:
            xlabels.append(" ")
    return xlabels


def plot_results_start(bam_file_path, bamname):
    codon_dict = fetch_vectors(bam_file_path)
    labels = scale_labels()

    norm = ""
    if global_args.normalize == True:
        norm = "_norm"

    for key in codon_dict:

        plt.title(
            ("Peak at: " + str(codon_dict[key].argmax() - (global_args.offset))))
        plt.grid(True, alpha=0.2)
        plt.step(np.linspace(-global_args.offset, global_args.offset, num=global_args.offset*2),
                 codon_dict[key], linewidth=0.5, color="red")
        plt.xticks(np.linspace(-global_args.offset, global_args.offset, num=global_args.offset*2),
                   labels, size="xx-small")
        plt.savefig(global_args.output_dir +
                    "coverage_codon/%s/%s%s.pdf" % (bamname, key, norm))
        plt.close()

    for key in aa_dict:
        collapsed_codons = []

        for x in aa_dict[key]:
            collapsed_codons.append(np.atleast_2d(codon_dict[x]))

        collapsed_aa = np.vstack(collapsed_codons).sum(axis=0)

        plt.title(("Peak at: " + str(collapsed_aa.argmax() - (global_args.offset))))
        plt.grid(True, alpha=0.2)
        plt.step(np.linspace(-global_args.offset, global_args.offset, num=global_args.offset*2),
                 collapsed_aa, linewidth=0.5, color="red")
        plt.xticks(np.linspace(-global_args.offset, global_args.offset, num=global_args.offset*2),
                   labels, size="xx-small")
        plt.savefig(global_args.output_dir +
                    "coverage_amino_acid/%s/%s%s.pdf" % (bamname, key, norm))
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
            os.mkdir(global_args.output_dir + "coverage_amino_acid/%s" % bamname)
            os.mkdir(global_args.output_dir + "coverage_codon/%s" % bamname)
        except:
            pass
        thread_args.append((bampath, bamname))

    pool.starmap_async(plot_results_start, thread_args)
    pool.close()
    pool.join()


def executable_2():

    for bampath, bamname in runall_samples(global_args.sample_dir):
        try:
            os.mkdir(global_args.output_dir + "coverage_amino_acid/%s" % bamname)
            os.mkdir(global_args.output_dir + "coverage_codon/%s" % bamname)
        except:
            pass

        plot_results_start(bampath, bamname)


if __name__ == "__main__":

    aa_dict = {
        "ALA": ["GCT", "GCC", "GCA", "GCG"],
        "ARG": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "ASN": ["AAT", "AAC"],
        "ASP": ["GAT", "GAC"],
        "CYS": ["TGT", "TGC"],
        "GLN": ["CAA", "CAG"],
        "GLU": ["GAA", "GAG"],
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

    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_dir", required=True)
    parser.add_argument("--genome_fasta", required=True)
    parser.add_argument("--annotation_file", required=True)
    parser.add_argument("--offset", type=int, default=50)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--gene_set")
    parser.add_argument("--normalize", type=bool, default=False)
    global_args = parser.parse_args()

    try:
        os.mkdir(global_args.output_dir + "coverage_amino_acid")
        os.mkdir(global_args.output_dir +"coverage_codon")

    except:
        pass

    if global_args.cores == 1:
        executable_2()
    else:
        executable()
