import numpy as np
from plastid import *
import dill
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#gtffile = list(GTF2_TranscriptAssembler(
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


alignments = BAMGenomeArray("/home/maria/Documents/pelechanolab/data/samples/2/post-proc/by14_clu_dedup_nomm.bam", mapping=FivePrimeMapFactory())
alignments.set_normalize(value=True)
count_vectors = []

for transcript in gtffile:
    if "RDN" in transcript.get_gene():
        continue
    if transcript.get_gene().startswith("t"):
        continue
    else:
        value_array = transcript.get_counts(alignments)
#        count_vectors.append(value_array[1:100])
        count_vectors.append(value_array[-99:])        

vector_array = np.vstack(count_vectors)
metagene = vector_array.sum(axis=0)


xlabels = []
for x in range(-49,50):
    if x%10 == 0:
        xlabels.append(x)
    else:
        xlabels.append(" ")  
plt.grid(True, alpha=0.3)
plt.step(np.linspace(-49, 50, num=99), metagene, linewidth=0.3, color="red")
plt.xticks(np.linspace(-49, 50, num=99), xlabels, size="xx-small")
plt.savefig("metagene.pdf")

print(metagene.argmax() - (offset-1))


