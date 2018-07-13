from Bio.SeqIO.QualityIO import FastqGeneralIterator

fastqset = set()
counter = 0

wh = open("unmapped.fasta", "w")

with open ("/home/maria/Documents/pelechanolab/data/samples/ATCC6633/umitools/star/S1_11_S11_R1_001.fastqUnmapped.out.mate1", "r") as fh:
    for title, seq, qual in FastqGeneralIterator(fh):
        counter += 1
        if seq not in fastqset:
            wh.write(">" + str(counter) + "\n" + seq + "\n")
        fastqset.add(seq)


print("Sequences: " + str(counter))
print("Unique: " + str(len(fastqset)))
        
