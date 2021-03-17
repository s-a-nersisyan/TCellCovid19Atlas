import sys
import os

netMHCpan_outdir = sys.argv[1]
print("HLA,peptide,affinity")
for fname in os.listdir(netMHCpan_outdir):
    netMHCpan_file = open(fname, "r")
    iterator = iter(netMHCpan_file)

    line = next(iterator)
    while not ("Pos" in line.split()):
        line = line(next)
    line = next(iterator)
    line = next(iterator)
    while not (line.startswith("-")):
        descript = line.split()
        print("{},{},{}".format(
            descript[1],
            descript[2],
            descript[14]
        ))
        line = next(iterator)
    netMHCpan_file.close()
