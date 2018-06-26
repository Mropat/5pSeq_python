import numpy as np
from plastid import *
import dill
import matplotlib
import matplotlib.pyplot as plt

gtffile = list(dill.load(open("gtf_assembled.sav", "rb")))
offset = 50
for transcript in gtffile:
    span = transcript.spanning_segment
    new_region = GenomicSegment(
        span.chrom, span.start - offset, span.start, span.strand)
    new_region_2 = GenomicSegment(
        span.chrom, span.end, span.end + offset, span.strand)
    transcript.add_segments(new_region, new_region_2)

fasta_dict = {}
with open("/home/maria/Documents/pelechanolab/data/R64e/genome.fa", "r") as genomeseq:
    chrom = ""
    genseq = ""
    genomelines = genomeseq.readlines()
    for ind, line in enumerate(genomelines):
        if line.startswith(">"):
            fasta_dict[chrom] = genseq
            chrom = line.strip()[1:]
            genseq = ""
        else:
            genseq = genseq + line.strip()
    fasta_dict[chrom] = genseq            


alignments = BAMGenomeArray("/home/maria/Documents/pelechanolab/data/samples/20160909YP_08-39666715/cutadapt/umitools/star/dedup/BY4741-t5-1_S8_L002_R1_001.fastqAligned.sortedByCoord.out.bam", mapping=FivePrimeMapFactory())
count_vectors = []

for transcript in gtffile:     
    if transcript.get_sequence(fasta_dict)[50:53] != "ATG":
        continue
    else:
        readseq = transcript.get_sequence(fasta_dict)
        for ind in range(50, len(readseq)-50, 3):
            codon = readseq[ind] + readseq[ind+1] + readseq[ind+2]
            if codon == "CCG":
                count_vectors.append(transcript.get_counts(alignments)[ind-50:ind+50])
            

vector_array = np.vstack(count_vectors)
metagene = vector_array.sum(axis=0)

xlabels = []
for x in range(-48,49):
    if x%10 == 0:
        xlabels.append(x)
    else:
        xlabels.append(" ")  

plt.grid(True, alpha=0.3)
plt.step(np.linspace(-48, 49, num=100), metagene, linewidth=0.5, color="red", fillstyle="full")
plt.xticks(np.linspace(-48, 49, num=100), xlabels, size="xx-small")
plt.savefig("codon.pdf")
