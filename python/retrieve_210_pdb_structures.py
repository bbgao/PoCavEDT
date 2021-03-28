# -*- coding: utf-8 -*-

import Bio
from Bio.PDB import PDBList
import os.path

# bound pdb id, unbound pdb id, ligand name
pdb_list = ["1a0q", "1a28", "1a42", "1a4g", "1a6w", "1a9u", "1aaq", "1abe", "1ac0", "1acj", "1aco", "1adb", "1add", "1adf", "1aec", "1aha", "1ai5", "1aj7", "1ake", "1anf", "1aoe", "1apt", "1ase", "1azm", "1b59", "1b6n", "1b9v", "1baf", "1bap", "1bcd", "1bgo", "1bhf", "1bl7", "1blh", "1bma", "1bmq", "1bra", "1byb", "1byg", "1c2t", "1c5c", "1c5x", "1c83", "1cbs", "1cbx", "1cdg", "1ckp", "1cla", "1cle", "1coy", "1cps", "1cqp", "1ctr", "1ctt", "1d0l", "1d3h", "1dbb", "1dd7", "1dg5", "1dhf", "1did", "1dih", "1dmp", "1dog", "1dr1", "1e96", "1eap", "1ebg", "1eed", "1ei1", "1ejn", "1ela", "1eoc", "1epb", "1eta", "1exw", "1f0r", "1fbl", "1fen", "1fgi", "1fkb", "1fki", "1fmo", "1frp", "1glp", "1gpy", "1hak", "1hbv", "1hdy", "1hew", "1hfc", "1hti", "1hyt", "1ibg", "1icn", "1ida", "1imb", "1inc", "1ivb", "1ivc", "1jao", "1l82", "1lah", "1lcp", "1ldm", "1lgr", "1lic", "1lmo", "1lpm", "1mbi", "1mfc", "1mmp", "1mmq", "1mrg", "1mrk", "1mts", "1mup", "1nco", "1nsc", "1okl", "1pbd", "1pdz", "1pgp", "1pha", "1poc", "1ppi", "1ppk", "1pso", "1qbr", "1qcf", "1qh7", "1qpe", "1rbp", "1rds", "1rgk", "1rne", "1rob", "1rpa", "1rt2", "1sln", "1slt", "1snc", "1sre", "1stp", "1tdb", "1thl", "1tlc", "1tng", "1tph", "1ukz", "1ulb", "1uvs", "1vgc", "1xid", "1ydr", "2aad", "2ack", "2ada", "2ak3", "2cmd", "2cpp", "2csc", "2ctc", "2er0", "2fox", "2gbp", "2gpb", "2ifb", "2msb", "2phh", "2pk4", "2qwb", "2sim", "2sns", "2tsc", "2xis", "2yhx", "2ypi", "3cla", "3dfr", "3er3", "3ert", "3fx2", "3gch", "3gpb", "3hvt", "3nos", "3ts1", "4cts", "4dfr", "4est", "4gr1", "4hvp", "4lbd", "4mbp", "4tln", "4xia", "5abp", "5cpp", "5er1", "5p21", "5p2p", "6acn", "6cpa", "6rnt", "6rsa", "7lpr", "7tim", "9aat", "9icd"];
pdbl = PDBList()


for bound_pdb in pdb_list:
    pdbl.retrieve_pdb_file(pdb_code=bound_pdb, obsolete=False, pdir='PDB', file_format="pdb")
    if not os.path.isfile("./PDB/pdb%s.ent" % bound_pdb) : 
        pdbl.retrieve_pdb_file(pdb_code=bound_pdb, obsolete=True, file_format="pdb")
