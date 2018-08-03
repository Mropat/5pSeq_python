

with open ("/home/maria/Documents/pelechanolab/data/l_reuteri/GCF_000159455.2_ASM15945v2_feature_table.txt", "r") as fh:
    lines = fh.read().split("\n")
    coding_ids = []
    for line in lines:
        if line.startswith("CDS"):
            contents = line.split("\t")
            if contents[1] == "with_protein":
                coding_ids.append(contents[10])

    with open ("coding_l_reuteri.txt", "w") as wh:
        wh.write("\n".join(coding_ids))
