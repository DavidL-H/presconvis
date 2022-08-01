# Root-mean-square deviation (RMSD) for proteins and/or ligands
# in PDB files.
# from: https://github.com/charnley/rmsd

# Set path to current working directory
import os
import calculate_rmsd

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)


pdb1 = "./457_A8_DTB_refine_23.pdb"
pdb2 = "./AF-A5ZUL4-F1-model_v3.pdb"

print(help(calculate_rmsd))