import pysam

wh = open("chr1reads.txt", "w")
with pysam.AlignmentFile("data/r1.bam", "rb") as r1_reads:
    for read in r1_reads.fetch():
        wh.write(read)
        wh.write("\n")
