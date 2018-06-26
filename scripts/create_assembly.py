from plastid import GTF2_TranscriptAssembler
import dill

gtffile = list(GTF2_TranscriptAssembler("/home/maria/Documents/pelechanolab/data/R64e/genes.gtf", return_type=Transcript))
dill.dump(gtffile, open("gtf_assembled.sav", "wb"))