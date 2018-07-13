import os

directory = "/home/maria/Documents/pelechanolab/data/samples/ATCC6633/cutadapt/umitools/star/"
summarypath = "/home/maria/Documents/pelechanolab/sample_reports/"


for filename in os.listdir(directory):
    if filename.endswith("final.out"):
        filepath = directory + filename
        with open (filepath, "r") as reportfile:
            reportlines = reportfile.read().split("\n")

            with open(summarypath + filename[:5], "w") as summaryfile:

                summaryfile.write(reportlines[5].strip().replace("|", " ") + "\n")
                summaryfile.write(reportlines[8].strip().replace("|", " ") + "\n")
                summaryfile.write(reportlines[9].strip().replace("|", " ") + "\n")
                summaryfile.write(reportlines[26].strip().replace("|", " ") + "\n")
                summaryfile.write(reportlines[29].strip().replace("|", " ") + "\n" + "\n")



