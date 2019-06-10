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

# calculates asa using dssp
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

# Generate contact network from ring files
def generate_contact_network(prot_id, dir_path = "../ring_files", contact_threshold = 8):
	# save contacts as a dictionary
	contacts = {}
	# with open(out_contact_file) as f:
	with open("{}/{}_edges.txt".format(dir_path, prot_id) ) as f:
		f.next()
		for line in f:
			# tokenize line
			node_1, edge, node_2, distance, _, _, atom_1, atom_2 = line.split()[0:8]
			# parse the distance as a float
			distance = float(distance)
			# if the distance is less thant contact threshold
			if distance <= contact_threshold:
				# set as a short contanct
				ls_distance = "short"
				# if the distance is more than half of contact threshold, set contact as long
				if distance > contact_threshold/2:
					ls_distance = "long"

				# parse edges
				edge_type, edge_loc = edge.split(":")
				# add contacts to contacts dictionary
				contacts.setdefault(node_1, [])
				contacts.setdefault(node_2, [])
				contacts[node_1].append((node_2, edge_loc, (atom_2, atom_1), edge_type, ls_distance))
				contacts[node_2].append((node_1, edge_loc, (atom_1, atom_2), edge_type, ls_distance))

	return (contacts)



#######################################################
#################  OTHER FUNCTIONS ####################
#######################################################

# extract long range features for a residues, given a window
def extract_long_range_features(residues, central_index, win_length):
	# compute window start and end indexes
	win_start = central_index - win_length
	if(win_start < 0):
		win_start = 0
	win_end = central_index + win_length
	if(win_end >= len(residues)):
		win_end = len(residues) - 1

	# compute the number of residues in the window
	seq_length = win_end - win_start + 1 
	# compute the distance between first and last residues of the window
	chain_dist = compute_distance(residues[win_start], residues[win_end]) 
	# compute the angle between start, middle, and end residues
	chain_angle = compute_angle(residues[win_start],residues[win_start + (win_end-win_start)/2], residues[win_end])

	# return the features as a list
	return [
		seq_length, 
		chain_dist/seq_length, 
		chain_angle,
		# adding a small amount to chain dist to avoid division by zero 
		(chain_angle/10)/(chain_dist + 0.001)
		]


# extract short range features for a residues, given a window and the contacts network
def extract_short_range_features(residues, central_index, win_length, contacts_network):
	# compute window start and end indexes
	win_start = central_index - win_length
	if(win_start < 0):
		win_start = 0
	win_end = central_index + win_length
	if(win_end >= len(residues)):
		win_end = len(residues) - 1

	# compute number of residues in window
	n_of_res = win_end - win_start + 1 

	# set to float to avoid approx when computing the mean
	inter_cc = 0.0
	intra_cc = 0.0
	long_cc  = 0.0
	short_cc = 0.0
	f_I   = 0.0
	f_II  = 0.0
	f_III = 0.0
	f_IV  = 0.0
	f_V   = 0.0

	# compute distance from start and end residues in the window
	dist = compute_distance(residues[win_start],residues[win_end])
	angle = compute_angle(residues[win_start],residues[win_start + (win_end-win_start)/2], residues[win_end])

	# for each residue in the window
	for index in range(win_start, win_end + 1):
		# get current residue
		curr_residue = residues[index]
		# compute it's ring id
		ring_id = residue_to_ring_id(curr_residue)
		# if the id is in contact network
		if ring_id in contacts_network: 
			# for each contact
			for contact in contacts_network[ring_id]:
				# count long and short contacts
				if contact[4] == "long":
					long_cc += 1
				if contact[4] == "short":
					short_cc += 1  
				# count inter and intra chain contacts  
				if chain_from_ring_id(contact[0]) == curr_residue.get_full_id()[2]:
					intra_cc += 1
				else:
					inter_cc += 1
		# compute factors
		f_I   += factor_I  [aa_3to1[curr_residue.get_resname()]]
		f_II  += factor_II [aa_3to1[curr_residue.get_resname()]]
		f_III += factor_III[aa_3to1[curr_residue.get_resname()]]
		f_IV  += factor_IV [aa_3to1[curr_residue.get_resname()]]
		f_V   += factor_V  [aa_3to1[curr_residue.get_resname()]]

	# return the features 
	return [
		# adding a small amount to intra to avoid division by zero 
		inter_cc/(intra_cc+1),
		intra_cc/n_of_res,
		inter_cc/n_of_res,
		long_cc/ n_of_res,
		short_cc/n_of_res,
		f_I  /n_of_res,
		f_II /n_of_res,
		f_III/n_of_res,
		f_IV /n_of_res,
		f_V  /n_of_res,
		dist/n_of_res,
		angle,
		# adding a small amount to dist to avoid division by zero 
		(angle/10)/(dist + 0.001)
	]

# Class that uses NS to find contacts, NOT USED
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


# get the spatial center of the residue, USED IN ContactCounter CLASS ONLY
def get_center(residue):
		cm = []
		for atom in residue:
			cm.append(atom.get_coord())
		cm = np.array(cm)
		return np.sum(cm, axis=0) / len(residue)


# compute distance between two residues
def compute_distance(res1, res2):
	# initialize coordinates variables
	a = 0
	b = 0

	# search for the Carbon Alpha atom and get coordinates
	for atom in res1:
		if atom.id == "CA":
			a = atom.get_coord()

	for atom in res2:
		if atom.id == "CA":
			b = atom.get_coord()

	# return norm as the distance
	return np.linalg.norm(a-b)


# compute angle between three residues
def compute_angle(res1, res2, res3):
	# initialize coordinates variables
	a = 0
	b = 0
	c = 0

	# search for the Carbon Alpha atom and get coordinates
	for atom in res1:
		if atom.id == "CA":
			a = atom.get_coord()
	for atom in res2:
		if atom.id == "CA":
			b = atom.get_coord()
	for atom in res3:
		if atom.id == "CA":
			c = atom.get_coord()

	# compute the edges ba and bc
	ba = a - b
	bc = c - b
	# with a very small window can happens that one of the edges is 0
	if (np.linalg.norm(ba) * np.linalg.norm(bc) == 0):
		# print a warning and return a 90 degrees angle
		print("Wrong arguments calculating angle between: {}, {}, {}".format(res1.get_full_id(),res2.get_full_id(),res3.get_full_id()))
		return 90
	# if not 0 compute and return angle
	cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
	angle = np.arccos(cosine_angle)
	return np.degrees(angle)



#######################################################
#################  DATASET CLASS  #####################
#######################################################

# Class that manages the lips_dataset.txt file, from parsing to features extraction
class ProteinDataset:
	# dictionary used to store ids and chains
	ds = {}
	# names of all the features computed
	# ORDER IS IMPORTANT
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
		"S_Dist", 
		"S_Ang", 
		"S_Ang/Dist",
		"L_Seq_Len", 
		"L_Dist/Seq_Len", 
		"L_Ang",
		"L_Ang/Dist"
		]

	# parse the lips_dataset from file
	def parse(self, path = "../sets/lips_dataset.txt"):
		dataset_file = open(path, 'r')
		lines = dataset_file.readlines()
		# for every line in the file (except first one)
		for line in lines[1:]:
			# split into tokens
			tokens = line.split()
			# set default entry for pdb id as an empty list
			self.ds.setdefault(tokens[0], [])
			# fill the list with a triple (CHAIN_ID, START, STOP)
			# START and STOP will both be 'neg' for counter example chains
			self.ds[tokens[0]].append((tokens[1], tokens[2], tokens[3]))
		# close file
		dataset_file.close()

	# check if a residue is labelled as lip
	def check_res(self, residue):
		# for every triple associated with residue's pdb id
		for lip in self.ds[residue.get_full_id()[0]]:
			# extract chain id, start and stop lip's indexes
			chain_id, start, stop = lip
			# if it is the same chain of the residue
			if chain_id == residue.get_full_id()[2]:
				# if start or stop are 'neg' it is not lip, return false
				if start == "neg" or stop == "neg":
					return False
				# if not returned yet, cast start and stop to integers and check if residue's index is in between 
				if int(start) <= residue.id[1] and int(stop) >= residue.id[1]:
					return True
		# if flow reaches here no entry is present for residue, return false
		return False

	# get a list of the chains of a pdb_structure that are labelled (with lips or 'neg')
	def get_labelled_chains(self, pdb_id):
		labelled_chains = []
		# for each triple associated with the id
		for entry in self.ds[pdb_id]:
			# append the first element (chain id)
			labelled_chains.append(entry[0])
		return labelled_chains

	# get list of the pdb structures loaded in this database
	def get_prot_list(self):
		return self.ds.keys()

	# download the given pdb files
	def download_pdb(self, ids, dir_path = "../pdb_files"):
		pdbl = PDBList()
		for prot_id in ids:
			pdbl.retrieve_pdb_file(prot_id, pdir=dir_path, file_format='pdb')

	# dowload ring files for the given ids, checks if already exists in ring folder
	def download_ring(self, ids, ring_folder="../ring_files"):
		for pdb_id in ids:
			if "{}_edges.txt".format(pdb_id) not in listdir(ring_folder):
				print "Downloading ring file for: {}".format(pdb_id)
				download_ring_file(pdb_id, ring_folder)
			else:
				print "Ring file for: {} exists".format(pdb_id)

	# extracts features for a real test, without lips labels and considering all the chains in the model
	def generate_blind_test(self, pdb_id, short_win, large_win, contact_threshold):
		# residues evaluated
		eval_res = []
		# feautres associated with evaluated residues
		X = []
		# network of the contacts in the file
		contacts = generate_contact_network(pdb_id, contact_threshold = contact_threshold)
		# pdb structure
		structure = PDBParser(QUIET=True).get_structure(pdb_id, "../pdb_files/pdb{}.ent".format(pdb_id))

		# for each chain in the structure
		for chain in structure[0]:
			# initialize a list that will contains the residues of the chain
			chain_res = []
			# for each residue check if it is a good residue an if it is append it
			for residue in chain:
				if residue.id[0] == ' ':
					chain_res.append(residue)

			# for each of the residues in the chain
			for res in chain_res:
				# set it as the center of a window
				central_index = chain_res.index(res)
				# extract long and short range features 
				sf = extract_short_range_features(chain_res, central_index, short_win, contacts )
				lf = extract_long_range_features(chain_res, central_index, large_win)
				# append residue's features to the list of computed features 
				X.append(sf+lf)
				# append the residue to the list of evaluated residues
				eval_res.append(res)

		# resturn evaluated residues, their features and a fake list of lip labels set to false
		return (eval_res, X, np.zeros(len(X)))

	# extracts features for a test used during model selection, with lips labels and considering only the labelled chains in the model
	def generate_test(self, pdb_id, short_win, large_win, contact_threshold):
		# residues evaluated
		eval_res = []
		# feautres associated with evaluated residues
		X = []
		# lips labels associated with evaluated residues
		y = []
		# network of the contacts in the file
		contacts = generate_contact_network(pdb_id, contact_threshold = contact_threshold)
		# pdb structure
		structure = PDBParser(QUIET=True).get_structure(pdb_id, "../pdb_files/pdb{}.ent".format(pdb_id))

		# for each chain in the structure
		for chain in structure[0]:
			# consider the chain only if it is labelled in the dataset
			if chain.id in self.get_labelled_chains(pdb_id):
			# initialize a list that will contains the residues of the chain
				chain_res = []
				# for each residue check if it is a good residue an if it is append it
				for residue in chain:
					if residue.id[0] == ' ':
						chain_res.append(residue)

				# for each of the residues in the chain
				for res in chain_res:
					# set it as the center of a window
					central_index = chain_res.index(res)
					# extract long and short range features 
					sf = extract_short_range_features(chain_res, central_index, short_win, contacts )
					lf = extract_long_range_features(chain_res, central_index, large_win)
					# append residue's features to the list of computed features 
					X.append(sf+lf)
					# generate the lip labels list checking residues
					if self.check_res(res):
						y.append(1)
					else:
						y.append(0)
					# append the residue to the list of evaluated residues
					eval_res.append(res)

		# resturn evaluated residues, their features and the list of lip labels
		return (eval_res, X, y)


	# extracts a given number of random examples from every given labelled chain 
	def generate_random_examples(self, pdb_ids, short_win, large_win, contact_threshold, ex_per_chain):
		# residues evaluated
		eval_res = []
		# feautres associated with evaluated residues
		X = []
		# lips labels associated with evaluated residues
		y = []
		# printing purpose only
		done_counter = 0


		# for each chain in the structure
		for pdb_id in pdb_ids:
			# printing purpose only
			done = done_counter*100/len(pdb_ids)
			sys.stdout.write("Generating random examples: {}%\r".format(done))
			sys.stdout.flush()
			
			# network of the contacts in the file
			contacts = generate_contact_network(pdb_id, contact_threshold = contact_threshold)
			# pdb structure
			structure = PDBParser(QUIET=True).get_structure(pdb_id, "../pdb_files/pdb{}.ent".format(pdb_id))

			# for each chain in the structure
			for chain in structure[0]:
				# consider the chain only if it is labelled in the dataset 
				if chain.id in self.get_labelled_chains(pdb_id):
					# initialize a list that will contains the residues of the chain
					chain_res = []
					# for each residue check if it is a good residue an if it is append it
					for residue in chain:
						if residue.id[0] == ' ':
							chain_res.append(residue)
					# select a certain number of random residues from the chain
					selected_res = self.select_random_residues(chain_res, ex_per_chain)

					# for each selectd residue
					for res in selected_res:
						# set it as the center of a window
						central_index = chain_res.index(res)
						# extract long and short range features 
						sf = extract_short_range_features(chain_res, central_index, short_win, contacts )
						lf = extract_long_range_features(chain_res, central_index, large_win)
						# append residue's features to the list of computed features 
						X.append(sf+lf)
						# generate the lip labels list checking residues
						if self.check_res(res):
							y.append(1)
						else:
							y.append(0)
						eval_res.append(res)
			# printing purpose only
			done_counter +=1
		# printing purpose only
		sys.stdout.write("Generating random examples: 100%")
		sys.stdout.flush()
		print

		# resturn evaluated residues, their features and the list of lip labels
		return (eval_res, X, y)


	# function used to select n residues from a list of residues
	def select_random_residues(self, residues, n):
		# if n is greater or equal to the number of given residues return all of them
		if n >= len(residues):
			return residues
		# else init a list of selected residues
		selected = []
		# while there's less than n residues selected
		while len(selected) < n:
			# select a random residue
			rand = random.randrange(len(residues))
			# if not in selected list yet, append it
			if not residues[rand] in selected:
				selected.append(residues[rand])
		# return selected residues
		return selected

	# transform features and lip labels into a pandas data frame
	def as_dataframe(self, X, y):
		pd_data = []
		# append lip flag at the end of features list for every example
		for i in range(0, len(X)):
			pd_data.append([])
			for item in X[i]:
				pd_data[i].append(item)
			pd_data[i].append(y[i])
		# create and return dataset
		return pd.DataFrame(data = pd_data, columns = self.features_names + ["y"])

	# parse a pandas dataframe from file
	def training_set_in(self, path ="../sets/training.txt"):
		return pd.read_csv(path)

	# write dataset to file using pandas dataframe
	def training_set_out(self, X, y, path ="../sets/training.txt"):
		self.as_dataframe(X,y).to_csv(path)
		return 

	# balance number of positive and negative examples in a dataset
	# prop is lower bound percent of positive cases
	def balance_neg_pos(self, new_res, new_fea, new_lip, prop = 50):
		# check if percent is respencted
		while(float(new_lip.count(1))/float(len(new_lip)) < prop ):
			# if not pop a random negative example from every list
			rand = random.randrange(len(new_lip))
			if new_lip[rand] == 0:
				new_res.pop(rand)
				new_fea.pop(rand)
				new_lip.pop(rand)

		# return new lists
		return(new_res, new_fea, new_lip)

