

from plastid import *

new_gtf = open("expanded_gtf.gtf","w")
extend_utr_length = 100 # change to whatever you like

ORFS = GTF2_TranscriptAssembler(open("/home/maria/Documents/pelechanolab/data/R64e/genes.gtf"))
for my_orf in ORFS:
    span = my_orf.spanning_segment
    if my_orf.strand == "+":
        new_region = GenomicSegment(span.chrom,span.start - extend_utr_length,span.end + extend_utr_length ,span.strand)

    else:
        new_region = GenomicSegment(span.chrom,span.start - extend_utr_length,span.end+extend_utr_length,span.strand)  
    new_transcript = SegmentChain(new_region, *my_orf.segments, **my_orf.attr) # copy metadata attributes from old ORF
    new_gtf.write(new_transcript.as_gtf())

new_gtf.close()            
            
