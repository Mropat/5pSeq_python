import numpy as np
from plastid import *
import dill
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# gtffile = list(GTF2_TranscriptAssembler(
#    "/home/maria/Documents/pelechanolab/data/R64e/genes.gtf", return_type=Transcript))
#dill.dump(gtffile, open("gtf_assembled.sav", "wb"))
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


alignments = BAMGenomeArray(
    "/home/maria/Documents/pelechanolab/data/samples/20160909YP_08-39666715/cutadapt/umitools/star/BY4741-t5-1_S8_L003_R1_001.fastqAligned.sortedByCoord.out.bam", mapping=FivePrimeMapFactory())
alignments_dedup = BAMGenomeArray(
    "/home/maria/Documents/pelechanolab/data/samples/20160909YP_08-39666715/cutadapt/umitools/star/dedup/BY4741-t5-1_S8_L003_R1_001.fastqAligned.sortedByCoord.out.bam", mapping=FivePrimeMapFactory())
count_vectors = []
count_vectors_dedup = []

for transcript in gtffile:
    if transcript.get_sequence(fasta_dict)[50:53] != "ATG":
        continue
    else:
        value_array = transcript.get_counts(alignments)
#        count_vectors.append(value_array[1:100])
        count_vectors.append(value_array[-99:])
        value_array = transcript.get_counts(alignments_dedup)
#        count_vectors_dedup.append(value_array[1:100])
        count_vectors_dedup.append(value_array[-99:])
        print(transcript.get_sequence(fasta_dict)[-53:-50])

vector_array = np.vstack(count_vectors)
metagene = vector_array.sum(axis=0)

vector_array_dedup = np.vstack(count_vectors_dedup)
metagene_dedup = vector_array_dedup.sum(axis=0)


xlabels = []
for x in range(-48, 49):
    if x % 10 == 0:
        xlabels.append(x)
    else:
        xlabels.append(" ")

plt.grid(True, alpha=0.3)
plt.step(np.linspace(-49, 51, num=99), metagene, linewidth=0.5, color="red")
plt.step(np.linspace(-49, 51, num=99),
         metagene_dedup, linewidth=0.5, color="blue")
plt.xticks(np.linspace(-49, 51, num=99), xlabels, size="xx-small")
plt.savefig("rdn_dedup_tail.pdf")

print(metagene.argmax() - (offset-1))
