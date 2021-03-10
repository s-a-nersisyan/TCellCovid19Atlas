import sys
import os

netMHCpan_outdir = sys.argv[1]
print("Allele,Peptide,Affinity")
for fname in sorted(os.listdir(netMHCpan_outdir)):
    netMHCpan_file = open(netMHCpan_outdir + "/" + fname, "r")
    table_started = False
    for line in netMHCpan_file:
        if "Pos" in line:
            table_started = True
            netMHCpan_file.readline()
            continue
        
        if table_started and line[0] == "-":
            break
        
        if not table_started:
            continue

        line = line.split()
        allele = line[1]
        peptide = line[2]
        affinity = int(float(line[15]))
        print("{},{},{}".format(allele, peptide, affinity))
    
    netMHCpan_file.close()
