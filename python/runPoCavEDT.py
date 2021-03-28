#!/usr/bin/env python
'''
Simple script that runs PoCavEDT on all proteins in the dataset. 
Please remember to set the paths to the executable and to the 
dataset folder according to their actual location on your computer.
The program output will be recorder in "results.txt". 
Use this file as input for "analyse_results.py" and "make_results_table.py".
'''
import os
from glob import glob
from subprocess import Popen, PIPE

with open("results.txt", "w") as out_results : 
    for protein_b in sorted(glob("../../Protein-Ligand_dataset_48_structures_bu/*_b.pdb")) : 
        pdb_id = os.path.basename(protein_b)[:4]
        protein_u = protein_b[:-5] + "u.pdb"
        ligand = protein_b[:-5] + "l.pdb"
        out_results.write("\n#########################\n")
        out_results.write("%s - bound\n" % pdb_id)
        out_results.write("#########################\n")
        bound = "../bin/PoCavEDT %s %s -r 5 --hetatm --evaluate_prediction --atom_radii ../bin/CHARMM22" %(protein_b, ligand)
        p = Popen(bound, shell=True, stdout=PIPE, stderr=PIPE)
        out_results.write(p.stdout.read())
        out_results.write("\n#########################\n")
        out_results.write("%s - unbound\n" % pdb_id)
        out_results.write("#########################\n")
        unbound = "../bin/PoCavEDT %s %s -r 5 --hetatm --evaluate_prediction --atom_radii ../bin/CHARMM22" %(protein_u, ligand)
        p = Popen(unbound, shell=True, stdout=PIPE, stderr=PIPE)
        out_results.write(p.stdout.read())