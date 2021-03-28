# -*- coding: utf-8 -*-

import Bio
from Bio.PDB import PDBList

# bound pdb id, unbound pdb id, ligand name
pdb_list = [("1BID", "3TMS", "UMP"), ("1CDO", "8ADH", "NAD"), ("1DWD", "1HXF", "MID"), ("1FBP", "2FBP", "F6P"), ("1GCA", "1GCG", "GAL"), ("1HEW", "1HEL", "NAG"), ("1HYT", "1NPC", "BZS"), ("1INC", "1ESA", "ICL"), ("1RBP", "1BRQ", "RTL"), ("1ROB", "8RAT", "C2P"), ("1STP", "1SWB", "BTN"), ("1ULB", "1ULA", "GUN"), ("2IFB", "1IFB", "PLM"), ("3PTB", "3PTN", "BEN"), ("2YPI", "1YPI", "PGA"), ("4DFR", "5DFR", "MTX"), ("4PHV", "3PHV", "VAC"), ("5CNA", "2CTV", "MMA"), ("7CPA", "5CPA", "FVF"), ("1A6W", "1A6U", "NIP"), ("1ACJ", "1QIF", "THA"), ("1APU", "3APP", "STA"), ("1BLH", "1DJB", "FOS"), ("1BYB", "1BYA", "GLC"), ("1HFC", "1CGE", "PLH"), ("1IDA", "1HSI", "PPL"), ("1IGJ", "1A4J", "DGX"), ("1IMB", "1IME", "LIP"), ("1IVD", "1NNA", "ST1"), ("1MRG", "1AHC", "ADN"), ("1MTW", "2TGA", "DX9"), ("1OKM", "4CA2", "SAB"), ("1PDZ", "1PDY", "PGA"), ("1PHD", "1PHC", "PIM"), ("1PSO", "1PSN", "STA"), ("1QPE", "3LCK", "PP2"), ("1RNE", "1BBS", "C60"), ("1SNC", "1STN", "PTP"), ("1SRF", "1PTS", "MTB"), ("2CTC", "2CTB", "LOF"), ("2H4N", "2CBA", "AZM"), ("2PK4", "1KRN", "ACA"), ("2SIM", "2SIL", "DAN"), ("2TMN", "1L3F", "LEP"), ("3GCH", "1CHG", "CIN"), ("3MTH", "6INS", "MPB"), ("5P2P", "3P2P", "DHG"), ("6RSA", "7RAT", "UVC")]

pdbl = PDBList()

for bound_pdb, unbound_pdb, _ in pdb_list:
    pdbl.retrieve_pdb_file(pdb_code=bound_pdb, pdir='PDB', file_format="pdb")
    pdbl.retrieve_pdb_file(pdb_code=unbound_pdb, pdir='PDB', file_format="pdb")
