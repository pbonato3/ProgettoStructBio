import pymol
from pymol import cmd, util

def visualize_features(prot_id, residues, lip_indexes, feature):
    pymol.finish_launching()  # Open Pymol (not necessary from pyMOL 2.1)
    cmd.fetch(prot_id, prot_id)  # Download the PDB

    cmd.hide("lines", prot_id)  # Hide lines
    cmd.show("cartoon", prot_id)  # Show ribbon

    cmd.alter("(all)", "b=0.0")  # Set all B-factor to 0.0

    pymol.finish_launching()  # Open Pymol (not necessary from pyMOL 2.1)


    for i in range(0,len(residues)):
        residue_string = '{}/{}{}/'.format(residues[i].get_full_id()[2], residues[i].id[1], '')  # Residue selection in Pymol syntax
        cmd.alter(residue_string, "b={:2f}".format(feature[i]))  # Substitute the B-score with the LIP score
        if lip_indexes[i] == 1:
            cmd.show_as("spheres", "chain {} and resi {}".format(residues[i].get_full_id()[2], residues[i].id[1]))

    cmd.spectrum("b", palette="rainbow", selection="(all)")  # color by B-factor values