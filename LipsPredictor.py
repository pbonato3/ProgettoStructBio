#!/usr/bin/env python

from os import listdir
import copy
import pandas as pd
import math
import random
import numpy as np
import requests
import json
import time
import zipfile
import pymol
from pymol import cmd, util
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser

### Kyle-Dolittle hydrophobicity
kd = { 
    'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
    'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
    'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
    'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 
}

aa_3to1 = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
}

solvent_accessibility = {
    'S':0.70,'T':0.71,'A':0.48,'G':0.51,'P':0.78,
    'C':0.32,'D':0.81,'E':0.93,'Q':0.81,'N':0.82,
    'L':0.41,'I':0.39,'V':0.40,'M':0.44,'F':0.42,
    'Y':0.67,'W':0.49,'K':0.93,'R':0.84,'H':0.66
}


#######################################################
#################  MAIN PROGRAM   #####################
#######################################################

from Dataset import ProteinDataset

dataset = ProteinDataset()
dataset.parse()
#dataset.download_all_pdb()
#dataset.download_all_ring_files()


prot_list = dataset.get_prot_list()

#long_prot =["1cee","1dev","1dow","1fqj","1g3j","1hrt","1i7w","1j2j","1jsu","1l8c","1p4q","1q68","1rf8","1sc5","1sqq","1tba","1th1","1xtg","1zoq","2a6q","2auh","2c1t","2o8a"]
#short_prot = []
#for prot in prot_list:
#    if prot not in long_prot:
#        short_prot.append(prot)
#
#prot_list_test = short_prot[6:-1] + long_prot[6:-1]
#prot_list_tr = short_prot[0:6] + long_prot[0:6]


residues, X_tr, y_tr= dataset.generate_random_examples(prot_list, 4, 4.5, 10)

print X_tr

dataset.training_set_out(X_tr, y_tr)



#print X_tr
#print y_tr
#print "#############################"
#print "           TRUTH"
#print "#############################"
#print y_test
#
#print "#############################"
#print "         PREDICTIONS"
#print "#############################"
#print predictions
#print "\n"
#print "accuracy:  {}".format(metrics.accuracy_score(y_test, predictions))
#print "precision: {}".format(metrics.precision_score(y_test, predictions))
#print "recall:    {}".format(metrics.recall_score(y_test, predictions))

#print clf.get_params()