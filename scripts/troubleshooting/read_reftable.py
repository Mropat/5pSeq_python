

with open ("/home/maria/Documents/pelechanolab/data/b_subtilis/6633/GCA_000177595.1_ASM17759v1_feature_table.txt", "r") as fh:
    lines = fh.read().split("\n")
    coding_ids = []
    for line in lines:
        if line.startswith("CDS"):
            contents = line.split("\t")
            if contents[1] == "with_protein":
                coding_ids.append(contents[10])

    with open ("coding_bsub_6633.txt", "w") as wh:
        wh.write("\n".join(coding_ids))
