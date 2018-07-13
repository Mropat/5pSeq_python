from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

gbank = SeqIO.parse(open("/home/maria/Documents/pelechanolab/data/l_plan/GCA_002631775.1_ASM263177v1_protein.gpff", "r"), "genbank")
output_handle=open("t_rRNAs.fa","w")
rRNAs = []
for genome in gbank :
    for feature in genome.features:
        if feature.type == "rRNA":
            ID = feature.qualifiers['db_xref'][0]
            desc = feature.qualifiers['locus_tag'][0]
            seq = feature.extract(genome.seq)
            record = SeqRecord(seq, id=ID, description=desc)
            rRNAs.append(record)

        elif feature.type == "tRNA":
            ID = feature.qualifiers['db_xref'][0]
            desc = feature.qualifiers['locus_tag'][0]
            seq = feature.extract(genome.seq)
            record = SeqRecord(seq, id=ID, description=desc)
            rRNAs.append(record)
            

SeqIO.write(rRNAs, output_handle, "fasta")
output_handle.close()

