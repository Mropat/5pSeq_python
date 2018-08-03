from plastid import * 
import dill

pickle_path = "plots/e_coli/genes.sav"
gtf_coords_file = list(dill.load(open(pickle_path, "rb")))
with open ("coding_list.txt", "w") as wh:
    for transcript in gtf_coords_file:
        if transcript.attr.get('type') == "mRNA":
            wh.write(transcript.get_name() + "\n")
