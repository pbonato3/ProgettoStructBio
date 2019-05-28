#!/usr/bin/env python

from os import listdir
import sys
import copy
import pandas as pd
import math
import random
import numpy as np
import requests
import json
import time
import zipfile
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser


aa_3to1 = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
}

# polarity index
factor_I = {
	'A': -0.591,'C': -1.343,'D':  1.050,'E':  1.357,'F': -1.006, 
	'G': -0.384,'H':  0.336,'I': -1.239,'K':  1.831,'L': -1.019, 
	'M': -0.663,'N':  0.945,'P':  0.189,'Q':  0.931,'R':  1.538, 
	'S': -0.228,'T': -0.032,'V': -1.337,'W': -0.595,'Y':  0.2
}

# secondary structure
factor_II = {
	'A': -1.302,'C':  0.465,'D':  0.302,'E': -1.453,'F': -0.590,
	'G':  1.652,'H': -0.417,'I': -0.547,'K': -0.561,'L': -0.987,
	'M': -1.524,'N':  0.828,'P':  2.081,'Q': -0.179,'R': -0.055,
	'S':  1.399,'T':  0.326,'V': -0.279,'W':  0.009,'Y':  0.830
}

# molecular size
factor_III = {
	'A': -0.733,'C': -0.862,'D': -3.656,'E':  1.477,'F':  1.891,
	'G':  1.330,'H': -1.673,'I':  2.131,'K':  0.533,'L': -1.505,
	'M':  2.219,'N':  1.299,'P': -1.628,'Q': -3.005,'R':  1.502,
	'S': -4.760,'T':  2.213,'V': -0.544,'W':  0.672,'Y':  3.097
}

# relative aa composition 
factor_IV = {
	'A':  1.570,'C': -1.020,'D': -0.259,'E':  0.113,'F': -0.397,
	'G':  1.045,'H': -1.474,'I':  0.393,'K': -0.277,'L':  1.266,
	'M': -1.005,'N': -0.169,'P':  0.421,'Q': -0.503,'R':  0.440,
	'S':  0.670,'T':  0.908,'V':  1.242,'W': -2.128,'Y': -0.838
}

# electrostatic charge
factor_V = {
	'A': -0.146,'C': -0.255,'D': -3.242,'E': -0.837,'F':  0.412,
	'G':  2.064,'H': -0.078,'I':  0.816,'K':  1.648,'L': -0.912,
	'M':  1.212,'N':  0.933,'P': -1.392,'Q': -1.853,'R':  2.897,
	'S': -2.647,'T':  1.313,'V': -1.262,'W': -0.184,'Y':  1.512
}


#######################################################
#################  RING FUNCTIONS #####################
#######################################################

# Calculate contacts with the web API (if the RING software is not working)
def download_ring_file(pdb_id, ring_folder):
    req = {"pdbName": pdb_id, "chain": "all", "seqSeparation": "5", "networkPolicy": "closest", "nowater": "true", "ringmd": "false", "allEdges": "true",
           "thresholds": '{"hbond": 3.5, "vdw": 0.5, "ionic": 4, "pipi": 6.5, "pication": 5, "disulphide": 2.5}'}
    r = requests.post('http://protein.bio.unipd.it/ringws/submit', data=json.dumps(req), headers={'content-type': 'application/json'})
    jobid = json.loads(r.text)['jobid']
    print jobid
    # Check job status and wait until done
    status = None
    while status != 'complete':
        status = json.loads(requests.get('http://protein.bio.unipd.it/ringws/status/{}'.format(jobid)).text).get("status", None)
        print status
        time.sleep(5)
    # Retrieve and write the web result to a file
    archive_file = "{}_network.zip".format(pdb_id)
    r = requests.get("http://protein.bio.unipd.it/ring_download/{}/{}".format(jobid, archive_file))
    with open(archive_file, "wb") as fout:
        fout.write(r.content)
    # Extract the contact file from the archive
    archive = zipfile.ZipFile(archive_file, 'r')
    with open("{}/{}_edges.txt".format(ring_folder,pdb_id), "w") as fout:
        fout.write(archive.read("{}_edges.txt".format(pdb_id)))

# Generate correct ring id from residue
def residue_to_ring_id(res):
    return "{}:{}:{}:{}".format(res.get_full_id()[2], res.id[1], res.id[2] if res.id[2] != ' ' else '_', res.get_resname())

# Get chain id from ring id
def chain_from_ring_id(res_id):
    return res_id.split(":")[0]

# Generate contact network for 
def generate_contact_network(prot_id, dir_path = "./ring_files", contact_threshold = 4):
    contacts = {}
    # with open(out_contact_file) as f:
    with open("{}/{}_edges.txt".format(dir_path, prot_id) ) as f:

        f.next()
        for line in f:

            # node = chain:residue_number:insertion_code:residue_name
            # edge = localization:contact_type (localization MC = main chain, SC = side chain)
            node_1, edge, node_2, distance, _, _, atom_1, atom_2 = line.split()[0:8]

            distance = float(distance)
            if distance <= contact_threshold:
                # print line
                edge_type, edge_loc = edge.split(":")
                contacts.setdefault(node_1, [])
                contacts.setdefault(node_2, [])
                contacts[node_1].append((node_2, edge_loc, (atom_2, atom_1), edge_type))
                contacts[node_2].append((node_1, edge_loc, (atom_1, atom_2), edge_type))

    return (contacts)



#######################################################
#################  OTHER FUNCTIONS ####################
#######################################################

def extract_features(residues, central_index, win_length, contacts_network, contact_threshold):
    X = []
    win_start = central_index - win_length
    if(win_start < 0):
        win_start = 0
    win_end = central_index + win_length
    if(win_end >= len(residues)):
        win_end = len(residues) - 1

    n_of_res = 0
    dist = compute_distance(residues[win_start],residues[win_end])

    # set to float to avoid approx
    inter_cc = 0.0
    intra_cc = 0.0

    chain_length = len(residues)
    chain_dist = compute_distance(residues[0], residues[-1]) 
    chain_angle = compute_angle(residues[0], residues[chain_length/2], residues[-1])

    f_I   = 0.0
    f_II  = 0.0
    f_III = 0.0
    f_IV  = 0.0
    f_V   = 0.0


    for index in range(win_start, win_end):
        curr_residue = residues[index]
        ring_id = residue_to_ring_id(curr_residue)
        if ring_id in contacts_network: 
            for contact in contacts_network[ring_id]:
                if chain_from_ring_id(contact[0]) == curr_residue.get_full_id()[2]:
                    intra_cc += 1
                else:
                    inter_cc += 1
        f_I   += factor_I  [aa_3to1[curr_residue.get_resname()]]
    	f_II  += factor_II [aa_3to1[curr_residue.get_resname()]]
    	f_III += factor_III[aa_3to1[curr_residue.get_resname()]]
    	f_IV  += factor_IV [aa_3to1[curr_residue.get_resname()]]
    	f_V   += factor_V  [aa_3to1[curr_residue.get_resname()]]
        n_of_res += 1

    X.append(inter_cc/(intra_cc+1))
    X.append(intra_cc/n_of_res)
    X.append(inter_cc/n_of_res)
    X.append(f_I  /n_of_res)
    X.append(f_II /n_of_res)
    X.append(f_III/n_of_res)
    X.append(f_IV /n_of_res)
    X.append(f_V  /n_of_res)
    X.append(dist/n_of_res)
    X.append(compute_angle(residues[win_start],residues[win_start + (win_end-win_start)/2], residues[win_end]))
    X.append(float(chain_length))
    X.append(chain_dist)
    X.append(chain_angle)

    return X


def compute_distance(res1, res2):
    a = 0
    b = 0

    for atom in res1:
        if atom.id == "CA":
            a = atom.get_coord()

    for atom in res2:
        if atom.id == "CA":
            b = atom.get_coord()

    return np.linalg.norm(a-b)


def compute_angle(res1, res2, res3):
    a = 0
    b = 0
    c = 0

    for atom in res1:
        if atom.id == "CA":
            a = atom.get_coord()

    for atom in res2:
        if atom.id == "CA":
            b = atom.get_coord()

    for atom in res3:
        if atom.id == "CA":
            c = atom.get_coord()

    ba = a - b
    bc = c - b

    if (np.linalg.norm(ba) * np.linalg.norm(bc) == 0):
        print("Wrong arguments calculating angle between: {}, {}, {}".format(res1.get_full_id(),res2.get_full_id(),res3.get_full_id()))
        return 90

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)



#######################################################
#################  DATASET CLASS  #####################
#######################################################


class ProteinDataset:
    ds = {}
    features_names = [
    	"Inter/Intra_CC", 
    	"Intra_CC/n_res",
    	"Inter_CC/n_res", 
    	"f_I/n_of_res",
    	"f_II/n_of_res",
    	"f_III/n_of_res",
    	"f_IV/n_of_res",
    	"f_V/n_of_res", 
    	"Dist/n_res", 
    	"Angle", 
    	"Chain_Length", 
    	"Chain_Dist", 
    	"Chain_Angle"
    	]

    def parse(self, path = "lips_dataset.txt"):
        dataset_file = open(path, 'r')
        lines = dataset_file.readlines()

        for line in lines[1:]:
            tokens = line.split()
            self.ds.setdefault(tokens[0], [])
            self.ds[tokens[0]].append((tokens[1], tokens[2], tokens[3]))

        dataset_file.close()

    def check_res(self, residue):
        for lip in self.ds[residue.get_full_id()[0]]:
            chain_id, start, stop = lip
            if chain_id == residue.get_full_id()[2]:
                if start == "neg" or stop == "neg":
                    return False
                if int(start) <= residue.id[1] and int(stop) >= residue.id[1]:
                    return True
        return False

    def get_labelled_chains(self, pdb_id):
    	labelled_chains = []
        for entry in self.ds[pdb_id]:
            labelled_chains.append(entry[0])
        return labelled_chains

    def get_prot_list(self):
        return self.ds.keys()

    def download_all_pdb(self, dir_path = "./pdb_files"):
        pdbl = PDBList()
        for prot_id in self.ds.get_prot_list():
            pdbl.retrieve_pdb_file(prot_id, pdir=dir_path, file_format='pdb')

    def download_all_ring_files(self, ring_folder="./ring_files"):
        for pdb_id in self.get_prot_list():
            print ("Calculating contacts for protein: {}".format(pdb_id))
            if not "{}_edges.txt".format(pdb_id) in listdir(ring_folder):
                download_ring_file(pdb_id, ring_folder)

    def generate_test(self, pdb_id, win_length, contact_threshold):
        eval_res = []
        X = []
        y = []
        contacts = generate_contact_network(pdb_id, contact_threshold = contact_threshold)
        structure = PDBParser(QUIET=True).get_structure(pdb_id, "./pdb_files/pdb{}.ent".format(pdb_id))
        for chain in structure[0]:
            if chain.id in self.get_labelled_chains(pdb_id):
                chain_res = []
                for residue in chain:
                    if residue.id[0] == ' ':
                        chain_res.append(residue)
                
                for res in chain_res:
                    central_index = chain_res.index(res)
                    X.append(extract_features(chain_res, central_index, win_length, contacts, contact_threshold))
    
                    if self.check_res(res):
                        y.append(1)
                    else:
                        y.append(0)
                    eval_res.append(res)

        return (eval_res, X, y)

    def generate_random_examples(self, pdb_ids, win_length, contact_threshold, ex_per_chain):
        eval_res = []
        X = []
        y = []
        done_counter = 0
        for pdb_id in pdb_ids:
            done = done_counter*100/len(pdb_ids)
            sys.stdout.write("Generating random examples: {}%\r".format(done))
            sys.stdout.flush()
            contacts = generate_contact_network(pdb_id, contact_threshold = contact_threshold)
            structure = PDBParser(QUIET=True).get_structure(pdb_id, "./pdb_files/pdb{}.ent".format(pdb_id))
            for chain in structure[0]:
                if chain.id in self.get_labelled_chains(pdb_id):
                    chain_res = []
                    for residue in chain:
                        if residue.id[0] == ' ':
                            chain_res.append(residue)
    
                    selected_res = self.select_random_residues(chain_res, ex_per_chain)
                
                    for res in selected_res:
                        central_index = chain_res.index(res)
                        X.append(extract_features(chain_res, central_index, win_length, contacts, contact_threshold))
    
                        if self.check_res(res):
                            y.append(1)
                        else:
                            y.append(0)
                        eval_res.append(res)
            done_counter +=1
        sys.stdout.write("Generating random examples: 100%")

        return (eval_res, X, y)



    def select_random_residues(self, residues, n):
        if n >= len(residues):
            n = len(residues)-1
        selected = []
        i = 0
        while i < n:
            rand = random.randrange(len(residues))
            if not residues[rand] in selected:
                selected.append(residues[rand])
                i += 1
        return selected


    def as_dataframe(self, X, y):
        pd_data = []
        for i in range(0, len(X)):
            pd_data.append([])
            for item in X[i]:
                pd_data[i].append(item)
            pd_data[i].append(y[i])
        print pd_data
        return pd.DataFrame(data = pd_data, columns = self.features_names + ["Lip_Flag"])


    def training_set_out(self, X, y, path ="./training.txt"):
    	self.as_dataframe(X,y).to_csv(path)

	def training_set_in(self, path ="./training.txt"):
	    df = pd.read_csv(path)
	    X_tr = df[features_names]
	    y_tr = df["Lip_Flag"]
	    return X_tr,y_tr


