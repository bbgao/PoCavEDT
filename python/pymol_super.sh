#!/bin/bash

# Script that performs the superimposition of the unbound protein over the bound one using pymol.

for ids in ./*_b.pdb
do
    pymol -c pymol_super.pml $ids "./$(basename $ids _b.pdb)_u.pdb"
done
