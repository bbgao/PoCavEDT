# -*- coding: utf-8 -*-
#!/usr/bin/env python
'''
Simple script that analyses the output of PoCavEDT and produces the results table for each PDB file.
PoCavEDT only outputs the top 3 pockets and top 3 cavities (if present), the script actually
merges the pockets and cavities into a single list and ranks them according to their distance 
from the protein center.
'''
import re

for f in ["other_results.txt"] :
    bound_results = { }
    unbound_results = { }
    with open(f) as input_results : 
        bound = False
        for l in input_results : 
            m=re.search(r"(....) - bound\n", l)
            if m :
                current_pdb_id = m.group(1)
                bound = True
                bound_results[current_pdb_id] = [[float("inf"), float("inf")], #pocket 1 
                                                 [float("inf"), float("inf")], #pocket 2
                                                 [float("inf"), float("inf")], #pocket 3
                                                 [float("inf"), float("inf")], #cavity 1
                                                 [float("inf"), float("inf")], #cavity 2
                                                 [float("inf"), float("inf")]] #cavity 3
            m=re.search(r"(....) - unbound\n", l)
            if m :
                current_pdb_id = m.group(1)
                bound = False
                unbound_results[current_pdb_id] = [[float("inf"), float("inf")], #pocket 1 
                                                 [float("inf"), float("inf")], #pocket 2
                                                 [float("inf"), float("inf")], #pocket 3
                                                 [float("inf"), float("inf")], #cavity 1
                                                 [float("inf"), float("inf")], #cavity 2
                                                 [float("inf"), float("inf")]] #cavity 3
            p = re.search(r"Pocket\[(\d)\] - weighted centroid distance from ligand: (.+)Å\.\n" % weighted, l)
            if p : 
                ii = int(p.group(1))
                d = float(p.group(2))
                if bound : 
                    bound_results[current_pdb_id][ii][0] = d
                else :
                    unbound_results[current_pdb_id][ii][0] = d
            p = re.search(r"Pocket\[(\d)\] - weighted centroid distance from protein centroid: (.+)Å\.\n" % weighted, l)
            if p : 
                ii = int(p.group(1))
                d = float(p.group(2))
                if bound : 
                    bound_results[current_pdb_id][ii][1] = d
                else :
                    unbound_results[current_pdb_id][ii][1] = d                

            c = re.search(r"Cavity\[(\d)\] - centroid distance from ligand: (.+)Å\.\n", l)
            if c : 
                ii = int(c.group(1)) + 3
                d = float(c.group(2))
                if bound : 
                    bound_results[current_pdb_id][ii][0] = d
                else :
                    unbound_results[current_pdb_id][ii][0] = d
            c = re.search(r"Cavity\[(\d)\] - centroid distance from protein centroid: (.+)Å\.\n", l)
            if c : 
                ii = int(c.group(1)) + 3
                d = float(c.group(2))
                if bound : 
                    bound_results[current_pdb_id][ii][1] = d
                else :
                    unbound_results[current_pdb_id][ii][1] = d                  
                
                
    for pdb_id in bound_results : 
        bound_results[pdb_id] = [y[0] for y in sorted(bound_results[pdb_id], key=lambda x: x[1])[:3]]
    for pdb_id in unbound_results : 
        unbound_results[pdb_id] = [y[0] for y in sorted(unbound_results[pdb_id], key=lambda x: x[1])[:3]]
    
    
    print "PDB_ID\trank\tdistance\trank\tdistance"
    
    for k in sorted(bound_results.keys()) : 
        b_rank = -1
        b_dist = -1
        if bound_results[k][0] <= 4.0 : 
            b_rank = 1
            b_dist = bound_results[k][0]
        elif bound_results[k][1] <= 4.0 : 
            b_rank = 2
            b_dist = bound_results[k][1]
        elif bound_results[k][2] <= 4.0 : 
            b_rank = 3
            b_dist = bound_results[k][2]

        u_rank = -1
        u_dist = -1
        if unbound_results[k][0] <= 4.0 : 
            u_rank = 1
            u_dist = unbound_results[k][0]
        elif unbound_results[k][1] <= 4.0 : 
            u_rank = 2
            u_dist = unbound_results[k][1]
        elif unbound_results[k][2] <= 4.0 : 
            u_rank = 3
            u_dist = unbound_results[k][2]
        print "%s\t%d\t%1.2f\t%d\t%1.2f" %(k, b_rank, b_dist, u_rank, u_dist)