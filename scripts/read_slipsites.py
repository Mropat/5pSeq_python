
slip_sites = {}

with open ("test.csv", "r") as fh:
    for line in fh:
        items = line.rstrip().split("\t")
        if items[0] not in slip_sites:
            slip_sites[items[0]] = set()
            slip_sites[items[0]].add(items[1])
        else:
            slip_sites[items[0]].add(items[1])

print(len(slip_sites))

