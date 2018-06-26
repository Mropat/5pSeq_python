import numpy as np 
from plastid import *
import matplotlib
import matplotlib.pyplot as plt 

gtffile = list(GTF2_TranscriptAssembler("/home/maria/Documents/pelechanolab/yal003w.gtf"))
extend_utr_length = 200
for my_orf in gtffile:
    span = my_orf.spanning_segment
    new_region = GenomicSegment(span.chrom,span.start - extend_utr_length, span.start,span.strand)
    new_region_2 = GenomicSegment(span.chrom,span.end, span.end + extend_utr_length,span.strand)
    my_orf.add_segments(new_region, new_region_2)

alignments = BAMGenomeArray("/home/maria/Documents/pelechanolab/data/samples/2/post-proc/by14_clu_dedup_nomm.bam", mapping=CenterMapFactoryDenorm())
count_vectors = []

for ind, transcript in enumerate(gtffile):
    count_vectors.append(transcript.get_counts(alignments))
        
    print(transcript)
    value_array=count_vectors[ind]
    plt.plot(value_array, linewidth=0.3, drawstyle="steps", fillstyle="bottom")
    plt.axvline(transcript.cds_start + extend_utr_length, linewidth = 0.1, color="red", ls="dotted")
    plt.axvline(transcript.cds_end + extend_utr_length, linewidth = 0.1, color="red", ls="dotted")
    plt.axvline(transcript.cds_start + extend_utr_length-17, linewidth = 0.1, color="blue", ls="dotted")
    plt.axhline(0, linewidth = 0.2, color="black", ls="dashed")
    plt.xlabel("Position in transcript (5' to 3')")
    plt.ylabel("Ribosome counts")
    plt.grid(True, lw=0.1, )
    plt.tight_layout(pad=0.1)
    plt.savefig("plt_%s.pdf"%transcript.get_name())

    alignments.set_mapping(FivePrimeMapFactory)
    
    value_array = transcript.get_counts(alignments)

    five_prime_max = max(value_array)
    print(five_prime_max)
    five_prime_max_pos= value_array.argmax()-200
    print(five_prime_max_pos)

