import os
import pandas as pd
os.getcwd()
import re


import re

ids_to_haplotypes = {}
current_id = ""
best1 = ""
best2 = ""
extracting_best1 = False
extracting_best2 = False

with open("haplotype_phase.out", "r") as f:
    for line in f:
        line = line.strip()
        
        if "BEGIN BESTPAIRS1" in line:
            extracting_best1 = True
        elif "END BESTPAIRS1" in line:
            extracting_best1 = False
        elif extracting_best1:
            if line.startswith("0 #"):
                current_id = line.strip()[3:]
            else:
                # If not starting with "0 #", assume it's the next line after the sample ID
                if not best1:
                    # Remove non-numeric characters from Best1 haplotype
                    best1 = re.sub(r'[^0-9]', '', line.strip())
                elif not best2:
                    # Remove non-numeric characters from Best2 haplotype
                    best2 = re.sub(r'[^0-9]', '', line.strip())
                    ids_to_haplotypes[current_id] = {"Best1": best1, "Best2": best2}
                    # Reset variables for next sample ID
                    best1 = ""
                    best2 = ""

# Write the extracted data into a file
with open("haplotypes_ttn_r2_0.8.txt", "w") as f:
    f.write("sample\tBest1\tBest2\n")
    for id, haplotypes in ids_to_haplotypes.items():
        f.write(f"{id}\t{haplotypes['Best1']}\t{haplotypes['Best2']}\n")




