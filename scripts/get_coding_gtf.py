from plastid import GTF2_TranscriptAssembler, Transcript
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--infile", required=True, help="Provide the *.gtf file")
parser.add_argument("--outfile", required=True)
global_args = parser.parse_args()

gtf_file = list(GTF2_TranscriptAssembler(global_args.infile, return_type=Transcript))

with open (global_args.outfile, "w") as wh:
    for transcript in gtf_file:
        if transcript.attr.get('type') == "mRNA":
            wh.write(transcript.get_name() + "\n")
