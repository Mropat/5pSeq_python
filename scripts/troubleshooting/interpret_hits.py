#Troubleshooting to interpret BLAST output:
#Data from BLAST query of coding nucleotide sequences against unique unmapped reads from alignment
#This script shows how many additional sequences could be mapped using more permissive mapping strategy

from Bio.Blast import NCBIXML
result_handle= open("/home/maria/Documents/pelechanolab/BLAST/unmapped_hits_reference", "r")
save_results = open("parsed_blast_rtrna.txt", "w")
blast_records = NCBIXML.parse(result_handle)

accounted_seqs = 0
unique_hits = 0
for record in blast_records:
    if len(record.alignments) > 0:
        accounted_seqs += len(record.alignments)
        unique_hits += 1
        save_results.write(record.query + "\n" + "Hits: " + str(len(record.alignments)) + "\n")

print("Accounted for:" + str(accounted_seqs))
print("Unique Hits:" + str(unique_hits))

result_handle.close
save_results.close
