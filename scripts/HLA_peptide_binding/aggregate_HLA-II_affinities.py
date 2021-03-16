import sys
import os

netMHCIIpan_outdir = sys.argv[1]
print("Allele,Peptide,Affinity")
for fname in sorted(os.listdir(netMHCIIpan_outdir)):
    netMHCIIpan_file = open(netMHCIIpan_outdir + "/" + fname, "r")
    table_started = False
    for line in netMHCIIpan_file:
        if "Pos" in line:
            table_started = True
            netMHCIIpan_file.readline()
            continue
        
        if table_started and line[0] == "-":
            break
        
        if not table_started:
            continue

        line = line.split()
        allele = line[1]
        peptide = line[2]
        affinity = int(float(line[11]))
        print("{},{},{}".format(allele, peptide, affinity))
    
    netMHCIIpan_file.close()
