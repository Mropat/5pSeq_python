import Bio
from Bio import SeqIO

misses = SeqIO.parse("/home/maria/Documents/pelechanolab/misses/misses_clust.fasta", "fasta")
for m in misses:
    sequence = m._seq
    if sequence.count("G") / len(sequence) > 0.33:
        continue
    elif len(sequence) < 100:
        continue
    else:
        with open("misses_clust_clean.txt", "a") as wh:
            SeqIO.write(m, wh, "fasta")
