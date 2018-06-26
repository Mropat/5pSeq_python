import numpy as np 
from plastid import *
import dill 

#gtffile = list(GTF2_TranscriptAssembler("/home/maria/Documents/pelechanolab/data/R64e/genes.gtf"))
#dill.dump(gtffile, open("gtf_assembled.sav", "wb"))
gtffile = list(dill.load(open("gtf_assembled.sav", "rb")))

extend_utr_length_up = 50
extend_utr_length_down = 100
for transcript in gtffile:
    span = transcript.spanning_segment
    new_region = GenomicSegment(span.chrom,span.start - extend_utr_length_up, span.start,span.strand)
    new_region_2 = GenomicSegment(span.chrom,span.end, span.end + extend_utr_length_down,span.strand)
    transcript.add_segments(new_region, new_region_2)


alignments = BAMGenomeArray("/home/maria/Documents/pelechanolab/data/samples/2/post-proc/by14_clu_dedup_nomm.bam", mapping=FivePrimeMapFactory())
count_vectors = []
pos_vectors = []

for transcript in gtffile:

    
    value_array = transcript.get_counts(alignments)
    count_vectors.append(value_array)

    five_prime_max = max(value_array)
    five_prime_max_pos= value_array.argmax()-extend_utr_length_up

    if five_prime_max > 5:
        pos_vectors.append(five_prime_max_pos)

print(pos_vectors)
print(np.mean(pos_vectors))

