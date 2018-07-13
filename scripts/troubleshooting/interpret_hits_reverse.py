#Troubleshooting script 2 to interpret BLAST output:
#Data from BLAST query of coding unique unmapped reads against unique unmapped sequences from alignment
#This script writes a fasta file of sequences which could not be mapped nor could provide a BLAST hit within the coding regions

from Bio.Blast import NCBIXML
from Bio import SeqIO


q_dict = SeqIO.index("/home/maria/Documents/pelechanolab/BLAST/unmapped.fasta", "fasta")

hits = []
for record in NCBIXML.parse(open("/home/maria/Documents/pelechanolab/BLAST/seqs_to_genome")):
    if record.alignments:
        hits.append(record.query.split()[0])

misses = set(q_dict.keys()) - set(hits)
orphan_records = [q_dict[name] for name in misses]

SeqIO.write(orphan_records, "misses.txt", "fasta")