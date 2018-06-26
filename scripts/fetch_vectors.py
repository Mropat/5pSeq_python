from plastid import BAMGenomeArray, FivePrimeMapFactory, BED_Reader, Transcript, CenterMapFactory
import numpy as np 
import matplotlib.pyplot as plt


roi_file = open("/home/maria/Documents/pelechanolab/bed_trans_dev.bed", "r", encoding="UTF-8")
transcripts = list(BED_Reader(roi_file))

alignments = BAMGenomeArray("/home/maria/Documents/pelechanolab/data/samples/2/post-proc/by14_clu_dedup_nomm.bam", mapping=FivePrimeMapFactory())
count_vectors = []

for transcript in transcripts:
    count_vectors.append(transcript.get_counts(alignments))


my_transcript = transcripts[156]
my_vector = count_vectors[156]

print(my_vector)
print(my_vector.shape)
print(my_vector.sum())

print(alignments[my_transcript])


