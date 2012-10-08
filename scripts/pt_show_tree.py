#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""Replace abbreviated name with long versions and display the tree."""
import sys

def main():
    """Main"""
    
    # mapping from short names to long names
    names = {
        "arath": "Arabidopsis thaliana",
        "carpa": "Carica papaya",
        "glyma": "Glycine max",
        "medtr": "Medicago truncatula",
        "poptr": "Populus trichocarpa",
        "ricco": "Ricinus communis",
        "vitvi": "Vitis vinifera",
        "orysa": "Oryza sativa",
        "sorbi": "Sorghum bicolor",
        "zeama": "Zea mays",
        "selmo": "Selaginella moellendorfii",
        "phypa": "Phscomitrella patens",
        "volca": "Volvox carteri",
        "chlre": "Chlamydomonas reinhardtii",
        "chlsp": "Chlorella sp.",
        "micp1": "Micromonas pusilla",
        "micp2": "Micromonas pusilla NOUM 17",
        "ostlu": "Ostreococcus lucimarinus",
        "ostta": "Ostreococcus tauri",
        "cyame": "Cyanidioschyzon merolae"  
    }

    # read in newick tree
    fp = open("../trees/outtree", "r")
    newick = "".join(fp.readlines())
    fp.close()
    
    # replace short names
    for k,v in names.items():
        newick = newick.replace(k, v)
        
    # write new version to a file
    fp = open("../trees/outtree.long", "w")
    fp.write(newick)
    fp.close()
    
if __name__ == '__main__':
    sys.exit(main())


