import numpy as np
import pysam
from plastid import *


#Write the indexing generator


wh = open("chr1reads_dev.txt", "w")

def get_aligns(filename):
    
    with pysam.AlignmentFile(filename, "rb") as ribotest:
        no_umi_list = set()
        segment_map = [[0, 0]]

        for pileupcolumn in ribotest.pileup("II"):
            
            unique_umi = set()
            base = set()

            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del == 1:
                    continue

                
                #if pileupread.alignment.qstart == 0:
                   #no_umi_list.add(pileupread.alignment.query_alignment_sequence)
                #will find alignments without 8 unmapped bases at the start

#                else:
#                    unique_umi.add(pileupread.alignment.query_name[-8:])
#                    base.add(
#                    pileupread.alignment.query_alignment_sequence[pileupread.query_position_or_next - pileupread.alignment.qstart])

            if len(unique_umi) > 0:
                if len(base) > 1:
                    base = "N"
                else:
                    base = "".join(base)
                mapped_pos = ([pileupcolumn.pos, len(unique_umi), base])
                if pileupcolumn.pos-1 not in segment_map[-1]:
                    yield ((segment_map[1:]))
                    segment_map = [[0, 0]]
                    segment_map.append((mapped_pos))

                else:
                    segment_map.append((mapped_pos))

            wh.write("\n" + "%s:%s" % (pileupcolumn.pos, len(unique_umi)))
        with open("no_umi_seqs.txt", "w") as exc_wh:
            exc_wh.write("\n".join(sorted(no_umi_list)))
        yield ((segment_map[1:]))




if __name__ == "__main__":
    aligns = {}
    filename = "/home/maria/Documents/pelechanolab/data/samples/2/post-proc/by14_clu_dedup_nomm.bam"
    for seg in get_aligns(filename):
        fullread = ("".join(x[-1] for x in seg))
        for ind, x in enumerate(seg):
            seg[ind] = x[:-1]
        aligns[fullread] = seg

    with open("dump_aligns.txt", "w") as segsave:
        for key in aligns:
            segsave.write(">" + key + "\n" + str(aligns[key]) + "\n" + "\n")

    print(len(aligns))