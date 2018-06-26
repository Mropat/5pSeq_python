#Creates tab delimited ROI + Gene ID for plastid

bed_ind = open("bed_trans_dev.gtf", "w")

with open ("data/R64e/genes.gtf", "r") as bed_handle:
    for line in bed_handle:
        line = line.split("\t")
        if line[7] == "transcript":
            line[1] = str(int(line[1]) - 65)
            line[2] = str(int(line[2]) + 65)
            line[4] = "0"
            line[6] = line[1]
            line[7] = line[2]
            bed_ind.write("\t".join(line[:8]) + "\n")

