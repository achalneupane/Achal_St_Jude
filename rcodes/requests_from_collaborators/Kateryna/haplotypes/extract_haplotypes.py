import os
import pandas as pd
os.getcwd()

        
ids_to_haplotypes = {}


# r2 > 0.8
with open("haplotype_input_edited_eur.txt.out", "r") as f:
    extract = False
    current_id = ""

    for line in f:
        if "BEGIN BESTPAIRS1" in line:
            extract = True
        elif "END BESTPAIRS1" in line:
            extract = False
        elif extract:
            if line.startswith("0 #"):
                current_id = line.strip()[3:]
            elif current_id:
                # Assume the haplotype is on the second line after the ID line,
                # and that line always starts with "0".
                haplotype_line = next(f)
                haplotype = haplotype_line.strip().split()[0:]
                ids_to_haplotypes[current_id] = " ".join(haplotype)

with open("haplotype_input_edited_eur_extracted.txt", "w") as f:
    for id, haplotype in ids_to_haplotypes.items():
        f.write(f"{id}\t{haplotype}\n")


