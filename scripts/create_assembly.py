# Script which reads and saves GTF hierarchy for reuse (optional with a lot of RAM)

from plastid import GTF2_TranscriptAssembler, Transcript, GenomicSegment, GFF3_TranscriptAssembler
import dill

#gtffile = list(GTF2_TranscriptAssembler("/home/maria/Documents/pelechanolab/data/e_coli/genes.gtf", return_type=Transcript))
gtffile = list(GFF3_TranscriptAssembler("/home/maria/Documents/pelechanolab/data/b_subtilis/6633/GCA_000177595.1_ASM17759v1_genomic.gff", return_type=Transcript))
dill.dump(gtffile, open("gtf_assembled_bsub_6633.sav", "wb"))