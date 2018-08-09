import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--infile", required=True, help="Provide file with the feature table from assembly")
parser.add_argument("--outfile", required=True)
global_args = parser.parse_args()


with open (global_args.infile, "r") as fh:
    lines = fh.read().split("\n")
    coding_ids = []
    for line in lines:
        if line.startswith("CDS"):
            contents = line.split("\t")
            if contents[1] == "with_protein":
                coding_ids.append(contents[10])

    with open (global_args.outfile, "w") as wh:
        wh.write("\n".join(coding_ids))


