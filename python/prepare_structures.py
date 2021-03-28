# -*- coding: utf-8 -*-

import os
from Bio.PDB import Select
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
pdb_parser = PDBParser(QUIET=True, PERMISSIVE=True)
io = PDBIO()

class LigandSelect(Select):
    def __init__(self, l_name):
        self.l_name = l_name
    def accept_residue(self, residue):
        if residue.get_resname()==self.l_name:
            return 1
        else:
            return 0

class NonHeteroSelect(Select):
    def accept_residue(self, residue):
        if residue.get_id()[0] == " " : 
            return 1
        else:
            return 0

class AllHeteroSelect(Select):
    def accept_residue(self, residue):
        if residue.get_id()[0] != " " : 
            return 1
        else:
            return 0

pdb_list = [("1BID", "3TMS", "UMP"), ("1CDO", "8ADH", "NAD"), ("1DWD", "1HXF", "MID"), ("1FBP", "2FBP", "F6P"), ("1GCA", "1GCG", "GAL"), ("1HEW", "1HEL", "NAG"), ("1HYT", "1NPC", "BZS"), ("1INC", "1ESA", "ICL"), ("1RBP", "1BRQ", "RTL"), ("1ROB", "8RAT", "C2P"), ("1STP", "1SWB", "BTN"), ("1ULB", "1ULA", "GUN"), ("2IFB", "1IFB", "PLM"), ("3PTB", "3PTN", "BEN"), ("2YPI", "1YPI", "PGA"), ("4DFR", "5DFR", "MTX"), ("4PHV", "3PHV", "VAC"), ("5CNA", "2CTV", "MMA"), ("7CPA", "5CPA", "FVF"), ("1A6W", "1A6U", "NIP"), ("1ACJ", "1QIF", "THA"), ("1APU", "3APP", "MAN"), ("1BLH", "1DJB", "FOS"), ("1BYB", "1BYA", "GLC"), ("1HFC", "1CGE", "PLH"), ("1IDA", "1HSI", "PPL"), ("1IGJ", "1A4J", "DGX"), ("1IMB", "1IME", "LIP"), ("1IVD", "1NNA", "ST1"), ("1MRG", "1AHC", "ADN"), ("1MTW", "2TGA", "DX9"), ("1OKM", "4CA2", "SAB"), ("1PDZ", "1PDY", "PGA"), ("1PHD", "1PHC", "PIW"), ("1PSO", "1PSN", "STA"), ("1QPE", "3LCK", "PP2"), ("1RNE", "1BBS", "C60"), ("1SNC", "1STN", "THP"), ("1SRF", "1PTS", "MTB"), ("2CTC", "2CTB", "HFA"), ("2H4N", "2CBA", "AZM"), ("2PK4", "1KRN", "ACA"), ("2SIM", "2SIL", "DAN"), ("2TMN", "1L3F", "0FA"), ("3GCH", "1CHG", "OAC"), ("3MTH", "6INS", "MPB"), ("5P2P", "3P2P", "DHG"), ("6RSA", "7RAT", "UVC")]

if not os.path.exists("./structures/"): os.makedirs("./structures/")

for bound_pdb_id, unbound_pdb_id, ligand_name in pdb_list : 
    bound_structure = pdb_parser.get_structure(bound_pdb_id, "./PDB/pdb%s.ent" % bound_pdb_id.lower())
    unbound_structure = pdb_parser.get_structure(unbound_pdb_id, "./PDB/pdb%s.ent" % unbound_pdb_id.lower())
    
    io.set_structure(bound_structure)
    io.save("./structures/%s_b.pdb" % bound_pdb_id, select=NonHeteroSelect(), preserve_atom_numbering=True)
    io.save("./structures/%s_l.pdb" % bound_pdb_id, select=LigandSelect(ligand_name), preserve_atom_numbering=True)
    io.save("./structures/%s_l2.pdb" % bound_pdb_id, select=AllHeteroSelect(), preserve_atom_numbering=True)

    io.set_structure(unbound_structure)
    io.save("./structures/%s_u.pdb" % bound_pdb_id, select=NonHeteroSelect(), preserve_atom_numbering=True)
