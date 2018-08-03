from plastid import GTF2_TranscriptAssembler, GFF3_TranscriptAssembler, Transcript, GenomicSegment, BAMGenomeArray, FivePrimeMapFactory

filename = "/home/maria/Documents/pelechanolab/data/samples/DH52/umitools/star/dedup/S1_14_S14_R1_001.fastqAligned.sortedByCoord.out.bam"
alignments = BAMGenomeArray(filename, mapping=FivePrimeMapFactory())

with open("ecoli.bedGraph", "w") as wh:
    alignments.to_bedgraph(wh, ecoli_all, "+")