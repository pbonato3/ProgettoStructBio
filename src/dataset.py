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
from Bio.PDB import PDBList, DSSP, NeighborSearch
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
#################  DSSP FUNCTIONS #####################
#######################################################

def get_asa(residues):
    pdb_id = residues[0].get_full_id()[0]
    chain_id = residues[0].get_full_id()[2]
    structure = PDBParser(QUIET=True).get_structure(pdb_id, "./pdb_files/pdb{}.ent".format(pdb_id))
    dssp = DSSP(structure[0], "./pdb_files/pdb{}.ent".format(pdb_id), dssp="bin/xssp-master/mkdssp")
    dssp = dict(dssp)  # Convert to dict to access residues

    tot_residues = len(residues)
    ss_content = {}
    surface = 0
    for residue in residues:
        if dssp.get((chain_id, residue.id)):
            tot_residues += 1
            dssp_index, aa, ss, asa, phi, psi, NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy, NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy = dssp.get((chain_id, residue.id))
            surface += asa
            ss_content.setdefault(ss, 0)
            ss_content[ss] += 1
    #for ss in ss_content:
        #print "{} {} {} {:.2f}".format(chain.id, ss, ss_content[ss], float(ss_content[ss]) / tot_residues)
    #print "{} ASA {:.2f} {:.2f}\n".format(chain.id, surface, float(surface) / tot_residues)
    return float(surface) / tot_residues

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
def generate_contact_network(prot_id, dir_path = "../ring_files", contact_threshold = 8):
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
                ls_distance = "short"
                if distance > contact_threshold/2:
                    ls_distance = "long"

                edge_type, edge_loc = edge.split(":")
                contacts.setdefault(node_1, [])
                contacts.setdefault(node_2, [])
                contacts[node_1].append((node_2, edge_loc, (atom_2, atom_1), edge_type, ls_distance))
                contacts[node_2].append((node_1, edge_loc, (atom_1, atom_2), edge_type, ls_distance))

    return (contacts)



#######################################################
#################  OTHER FUNCTIONS ####################
#######################################################

def extract_chain_features(chain_res):
    seq_length = len(chain_res)
    chain_dist = compute_distance(chain_res[0], chain_res[-1]) 
    chain_angle = compute_angle(chain_res[0], chain_res[seq_length/2], chain_res[-1])
    #chain_asa = get_asa(chain_res)
    return [
        seq_length, 
        chain_dist/seq_length, 
        chain_angle, 
        (chain_dist/seq_length)/chain_angle
        ]



def extract_features(residues, central_index, win_length, contacts_network):
    X = []
    win_start = central_index - win_length
    if(win_start < 0):
        win_start = 0
    win_end = central_index + win_length
    if(win_end >= len(residues)):
        win_end = len(residues) - 1

    n_of_res = 0

    # set to float to avoid approx
    inter_cc = 0.0
    intra_cc = 0.0
    long_cc  = 0.0
    short_cc = 0.0
    f_I   = 0.0
    f_II  = 0.0
    f_III = 0.0
    f_IV  = 0.0
    f_V   = 0.0

    dist = compute_distance(residues[win_start],residues[win_end])
    angle = compute_angle(residues[win_start],residues[win_start + (win_end-win_start)/2], residues[win_end])

    for index in range(win_start, win_end):
        curr_residue = residues[index]
           

        ring_id = residue_to_ring_id(curr_residue)
        if ring_id in contacts_network: 

            for contact in contacts_network[ring_id]:
                if contact[4] == "long":
                    long_cc += 1
                if contact[4] == "short":
                    short_cc += 1    
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
    X.append(long_cc/ n_of_res)
    X.append(short_cc/n_of_res)
    X.append(f_I  /n_of_res)
    X.append(f_II /n_of_res)
    X.append(f_III/n_of_res)
    X.append(f_IV /n_of_res)
    X.append(f_V  /n_of_res)
    X.append(dist/n_of_res)
    X.append(angle)
    X.append(angle/(dist + 0.001))

    return X

class ContactCounter:
    contacts = {}

    def __init__(self, structure, chain_id, threshold):
        ns = NeighborSearch(list(structure.get_atoms()))

        for chain in structure[0]:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[0] == ' ':
                        inter = 0
                        intra = 0
                        lon = 0
                        shor = 0
    
    
                        center = get_center(residue)
                        for a in ns.search(center, threshold):  # Iterate over contacts
                            if a.get_full_id()[2] == residue.get_full_id()[2]:
                                #if abs(int(a.get_full_id()[1]) - int(residue.get_full_id()[1])) > 3:
                                intra += 1
                            else:
                                inter += 1
                            if np.linalg.norm(a.get_coord()-center) < threshold/2:
                                shor += 1
                            else:
                                lon += 1
                        self.contacts[residue.get_full_id()] = (inter, intra, lon, shor)


    def get_contacts(self, res):
        return self.contacts[res.get_full_id()]

    def get_inter(self, res):
        return self.contacts[res.get_full_id()][0]

    def get_intra(self, res):
        return self.contacts[res.get_full_id()][1]

    def get_long(self, res):
        return self.contacts[res.get_full_id()][2]

    def get_short(self, res):
        return self.contacts[res.get_full_id()][3]

def get_center(residue):
        cm = []
        for atom in residue:
            cm.append(atom.get_coord())
        cm = np.array(cm)
        return np.sum(cm, axis=0) / len(residue)


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
    	"IrIa_CC", 
    	"Intra",
    	"Inter",
        "L_CC",
        "S_CC",
    	"f_I",
    	"f_II",
    	"f_III",
    	"f_IV",
    	"f_V", 
    	"Dist", 
    	"Ang", 
        "Ang/Dist",
    	"Seq_Len", 
    	"Ch_Dist/Seq_Len", 
    	"Ch_Ang",
        "ChD/SqL/An"
    	]

    def parse(self, path = "../sets/lips_dataset.txt"):
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

    def download_pdb(self, ids, dir_path = "../pdb_files"):
        pdbl = PDBList()
        for prot_id in ids:
            pdbl.retrieve_pdb_file(prot_id, pdir=dir_path, file_format='pdb')

    def download_ring(self, ids, ring_folder="../ring_files"):
        for pdb_id in ids:
            if "{}_edges.txt".format(pdb_id) not in listdir(ring_folder):
                print "Downloading ring file for: {}".format(pdb_id)
                download_ring_file(pdb_id, ring_folder)
            else:
                print "Ring file for: {} exists".format(pdb_id)

    def generate_blind_test(self, pdb_id, path, win_length, contact_threshold):
        eval_res = []
        X = []
        y = []
        contacts = generate_contact_network(pdb_id, contact_threshold = contact_threshold)
        structure = PDBParser(QUIET=True).get_structure(pdb_id, "../pdb_files/pdb{}.ent".format(pdb_id))

        for chain in structure[0]:
            chain_res = []
            for residue in chain:
                if residue.id[0] == ' ':
                    chain_res.append(residue)

            cf = extract_chain_features(chain_res)

            for res in chain_res:
                central_index = chain_res.index(res)
                f = extract_features(chain_res, central_index, win_length, contacts )
                X.append(f+cf)

        return (eval_res, X, y)

    def generate_test(self, pdb_id, win_length, contact_threshold):
        eval_res = []
        X = []
        y = []

        contacts = generate_contact_network(pdb_id, contact_threshold = contact_threshold)
        structure = PDBParser(QUIET=True).get_structure(pdb_id, "../pdb_files/pdb{}.ent".format(pdb_id))

        for chain in structure[0]:
            if chain.id in self.get_labelled_chains(pdb_id):
                chain_res = []
                for residue in chain:
                    if residue.id[0] == ' ':
                        chain_res.append(residue)

                cf = extract_chain_features(chain_res)

                for res in chain_res:
                    central_index = chain_res.index(res)
                    f = extract_features(chain_res, central_index, win_length, contacts)
                    X.append(f+cf)

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
            structure = PDBParser(QUIET=True).get_structure(pdb_id, "../pdb_files/pdb{}.ent".format(pdb_id))
            for chain in structure[0]:
                if chain.id in self.get_labelled_chains(pdb_id):
                    chain_res = []
                    for residue in chain:
                        if residue.id[0] == ' ':
                            chain_res.append(residue)
    
                    selected_res = self.select_random_residues(chain_res, ex_per_chain)
                    # chain features
                    cf = extract_chain_features(chain_res)

                    for res in selected_res:
                        central_index = chain_res.index(res)
                        f = extract_features(chain_res, central_index, win_length, contacts )
                        X.append(f+cf)
    
                        if self.check_res(res):
                            y.append(1)
                        else:
                            y.append(0)
                        eval_res.append(res)
            done_counter +=1
        sys.stdout.write("Generating random examples: 100%")
        print

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
        return pd.DataFrame(data = pd_data, columns = self.features_names + ["y"])

    def training_set_in(self, path ="../sets/training.txt"):
        return pd.read_csv(path)

    def training_set_out(self, X, y, path ="../sets/training.txt"):
    	self.as_dataframe(X,y).to_csv(path)
        return 


    def balance_neg_pos(residues, features, lip, prop = 50):
        if prop < 50:
            prop = 50
        new_res = residues
        new_fea = features
        new_lip = lip
        
        tot = len(new_lip)
        pos_count = new_lip.count(1)
        neg_count = new_lip.count(0)
        
        while(float(neg_count)/float(tot) > prop or float(pos_count)/float(tot) > prop ):
            if (neg_count > pos_count):
                rand = random.randrange(len(new_lip))
                if new_lip[rand] == 0:
                    new_res.pop(rand)
                    new_fea.pop(rand)
                    new_lip.pop(rand)
                    neg_count -= 1
            else:
                rand = random.randrange(len(new_lip))
                if new_lip[rand] == 1:
                    new_res.pop(rand)
                    new_fea.pop(rand)
                    new_lip.pop(rand)
                    pos_count -= 1
    
        return(new_res, new_fea, new_lip)

