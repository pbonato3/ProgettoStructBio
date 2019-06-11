#!/usr/bin/env python
from dataset import ProteinDataset, aa_3to1
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn import metrics
import copy
import pickle
import random
import sys
from os.path import isfile
from os import getcwd

# find correct paths
main_folder_path = './'
src_folder_path = './src'

if getcwd().endswith("/src"):
    main_folder_path = '../'
    src_folder_path = './'

pdb_folder_path  = main_folder_path+'pdb_files/'
ring_folder_path = main_folder_path+'ring_files/'
sets_folder_path = main_folder_path+'sets/'
test_folder_path = main_folder_path+'tests/'

#default_training_ids = ['2ivz','1dow','1hrt','1i7w','1j2j','1l8c','1rf8','1sqq','3b71','1a3b','1hv2','1ozs','1i8h','1axc','2gl7','1h2k','1ycq','1fv1','1kdx','1cqt']
default_training_ids = ["1p22","1ozs","2gsi","1fqj","1o9a","1kdx","1i7w","1hv2","1dev","1tba","1sc5","1lm8","1sb0","2phe","1i8h","1fv1","1l8c","2o8a","2gl7","1rf8","1cqt","2nl9","1hrt"]
default_choosen_features = [
			"IrIa_CC", 
			"Intra",
			"Inter",
			"S_Dist", 
			"S_Ang", 
			"S_Ang/Dist",
			"L_Dist/Seq_Len", 
			]

def test(clf, test_name,  sw = 4, lw = 60, ct = 6, epc = 100, choosen_features = None, tr_prots = None):
	print "###############################################################################################"
	print "RUNNING TEST {}".format(test_name)
	print "\n"
	# initilize a protein dataset object and parse ids
	prot_dataset = ProteinDataset()
	prot_dataset.parse()

	######## DEFAULT FEATURES FOR THIS TEST #########
	if not choosen_features:
		choosen_features = default_choosen_features

	######## DEFAULT PDB IDS FOR THIS TEST #########
	if not tr_prots:
		tr_prots = default_training_ids
	######## TEST SET GENERATION #########

	if isfile(test_folder_path+"{}_trs.txt".format(test_name)):
		print "Dataset found"
	else:
		res, X, y = prot_dataset.generate_random_examples(tr_prots, short_win = sw, large_win = lw, contact_threshold = ct, ex_per_chain = epc)
		prot_dataset.training_set_out(X,y, path = test_folder_path+"{}_trs.txt".format(test_name))

	df = prot_dataset.training_set_in(path = test_folder_path+"{}_trs.txt".format(test_name))

	# fit classifier
	clf.fit(df[choosen_features],df['y'])

	# save classifier to file
	model_out(clf, test_folder_path+"{}.sav".format(test_name))

	# select as test ids all ids in protein_dataset that are not in test set
	test_prots = []
	for prot_id in prot_dataset.get_prot_list():
		if prot_id not in tr_prots:
			test_prots.append(prot_id)

	# collect all partial results from 
	avg_pred = []
	avg_y = []

	# clear output file if exists
	clear_result_file(test_folder_path+"{}_res.txt".format(test_name))

	# foreach test pdb id 
	for t_p in test_prots:
		# generate test dataset
		res_test, X_test, y_test = prot_dataset.generate_test(t_p, short_win = sw, large_win = lw, contact_threshold = ct)
		# convert dataset to dataframe
		df_test = prot_dataset.as_dataframe(X_test, y_test)
		# predict 
		predictions = clf.predict_proba(df_test[choosen_features])
		# classifier returns a list of pairs [a, b] where 'a' is the 
		# probability of being class '0' and 'b' the probability of being class '1'
		predictions = [pred[1] for pred in predictions]
		# applu blur over predictions
		predictions = blur_by_chain(res_test, predictions)

		# output predictions to file
		out_predictions(t_p, res_test, predictions, path = test_folder_path+"{}_res.txt".format(test_name))

		# compute binary predictions
		bin_pred = prob_to_binary(predictions)

		# print results for every pdb id alone
		print ">{}".format(t_p)
		print 
		print "accuracy:          {:.5f}".format(metrics.accuracy_score(y_test, bin_pred))
		print "balanced accuracy: {:.5f}".format(metrics.balanced_accuracy_score(y_test, bin_pred))
		print "precision:         {:.5f}".format(metrics.precision_score(y_test, bin_pred))
		print "recall:            {:.5f}".format(metrics.recall_score(y_test, bin_pred))
		print "f1 score:          {:.5f}".format(metrics.f1_score(y_test, bin_pred))
		# this metrics throws an error if the test set is one class only
		if y_test.count(0) == len(y_test) or y_test.count(0) == 0:
			print "ROC AUC:           {}".format("Not defined for one class only")
		else:
			print "ROC AUC:           {:.5f}".format(metrics.roc_auc_score(y_test, bin_pred))
		print


		# merge lists for the average metrics
		avg_y = avg_y + list(y_test)
		avg_pred = avg_pred + list(bin_pred)


	# make a dictionary form metrics
	results = {}

	results["acc"] = metrics.accuracy_score(avg_y, avg_pred)
	results["bac"] = metrics.balanced_accuracy_score(avg_y, avg_pred)
	results["pre"] = metrics.precision_score(avg_y, avg_pred)
	results["rec"] = metrics.recall_score(avg_y, avg_pred)
	results["f1s"] = metrics.f1_score(avg_y, avg_pred)
	results["roc"] = metrics.roc_auc_score(avg_y, avg_pred)

	# print average results
	print "#############################"
	print "       AVERAGE RESULTS       "
	print "#############################"
	print 
	print "accuracy:          {:.5f}".format(results["acc"])
	print "balanced accuracy: {:.5f}".format(results["bac"])
	print "precision:         {:.5f}".format(results["pre"])
	print "recall:            {:.5f}".format(results["rec"])
	print "f1 score:          {:.5f}".format(results["f1s"])
	print "ROC AUC:           {:.5f}".format(results["roc"])
	print

	return results


# not used
def fill_1_gaps(pred, gap_len = 6):
	for idx in range(0, len(pred)):
		start = max(idx - gap_len, 0)
		stop = min(idx + gap_len, len(pred)-1)
		zeros_count = 0
		for label in pred[start:stop+1]:
			if label == 0:
				zeros_count += 1
		if zeros_count <= gap_len:
			pred[idx] = 1

	return pred	

# not used
def fill_0_gaps(pred, gap_len = 6):
	for idx in range(0, len(pred)):
		start = max(idx - gap_len, 0)
		stop = min(idx + gap_len, len(pred)-1)
		zeros_count = 0
		for label in pred[start:stop+1]:
			if label == 0:
				zeros_count += 1
		if zeros_count > gap_len:
			pred[idx] = 0

	return pred	

# apply blur filer dividing predictions by chain
def blur_by_chain(residues, predictions):
	# it is better to perform blur over chain instead of entire result
	# so extract chain labels and keep them ordered with a list
	chains_labels = []
	# make a dictionary  chain-label --> [start-index, end-index]
	chains = {}
	for i in range(0,len(residues)):
		if residues[i].get_full_id()[2] not in chains_labels :
			chains.setdefault(residues[i].get_full_id()[2], [i, i])
			chains_labels.append(residues[i].get_full_id()[2])
		chains[residues[i].get_full_id()[2]][1] =  i
	new_pred = []
	# perform blur and merge all together again
	for label in chains_labels:
		first = chains[label][0]
		last = chains[label][1]
		new_pred = new_pred + prob_blur(predictions[first : last + 1])
	return new_pred

# apply a blur filter on probability output, win_len*2 +1
def prob_blur(pred, win_len = 6):
	# initialize a new list
	new = []
	# loop on predictions
	for idx in range(0, len(pred)):
		# calculate the first and the last index of the window
		start = max(idx - win_len, 0)
		stop = min(idx + win_len, len(pred)-1)
		# sum all values in the window
		prob_count = 0.0
		for label in pred[start:stop+1]:
			prob_count += label
		# append to the new list the mean of the probability
		new.append(prob_count/(stop - start + 1))
	# return new list
	return new	

# convert probability prediction to binary prediction (if prob >= 0.5 then 1 else 0)
def prob_to_binary(pred):
	new = []
	for pr in pred:
		if pr >= 0.5:
			new.append(1)
		else:
			new.append(0)
	return new

# the file is opened in 'append' mode so can be useful to clear it at the beginning
def clear_result_file(path = main_folder_path+"results.txt"):
	if isfile(path):
		file = open(path, 'w')
		file.close()

# output predictions to file
def out_predictions(p_id, residues, predictions, path = main_folder_path+"results.txt"):
	# compute the binary classification
	bin_pred = prob_to_binary(predictions)
	# open file in append mode
	file = open(path, 'a')
	# write the pdb id
	file.write(">{}\n".format(p_id))
	# for each residue, print the id and probability in given format
	for idx in range(0,len(residues)):
		res_id = residues[idx].get_full_id()
		insertion_code = res_id[3][2]
		if insertion_code == ' ':
			insertion_code = ''
		file.write("{}/{}/{}/{}/{} {:.3f} {}\n".format(res_id[1],res_id[2],res_id[3][1],insertion_code,aa_3to1[residues[idx].get_resname()],predictions[idx],bin_pred[idx]))

	file.close()

# make a decision tree from configuration
def make_dt(config):
	return DecisionTreeClassifier(max_depth = config["max-depth"], random_state = config["random-state"])

# make a random forest from configuration
def make_rand_forest(config):
	return RandomForestClassifier(n_estimators = config["n-estimators"], max_depth = config["max-depth"], random_state = config["random-state"])

# make a multi layer perceptron from configuration
def make_mlp(config):
	return MLPClassifier(
		hidden_layer_sizes = config["hidden-layers"], 
		activation = config["activation"], 
		solver = config["solver"],
		alpha = config["alpha"],
		learning_rate_init = config["learning-rate-init"],
		max_iter = config["max-iter"],
		tol = config["tol"],
		random_state = config["random-state"])
		


# make a classifier of given type and configuration, train it ad return it
def make_predictor(model_type, config, training_set, features):
	if model_type == 'decision-tree':
		clf = make_dt(config)
		return clf.fit(training_set[features],training_set['y'])
	elif model_type == 'random-forest':
		clf = make_rand_forest(config)
		return clf.fit(training_set[features],training_set['y'])

	elif model_type == 'multi-layer-perceptron':
		clf = make_mlp(config)
		return clf.fit(training_set[features],training_set['y'])
	return 

# output the trained model to file
def model_out(model, path = main_folder_path+"trained-model.sav"):
	pickle.dump(model, open(path, 'wb'))
	return 

# parse trained model from file
def model_in(path = main_folder_path+"trained-model.sav"):
	return  pickle.load(open(path, 'rb'))


# predict a list of pdb ids using the given classifier and parameters. then output the results to file
def predict(clf, pdb_ids, features, short_win, large_win, contact_threshold, path = main_folder_path+"results.txt", blur = True, blur_w = 6):
	# initialize a proteinDataset object
	prot_dataset = ProteinDataset()

	# clear the restult file if exists
	clear_result_file(path)

	# printing purpose only
	done_counter = 0
	# loop on ids
	for pdb_id in pdb_ids:
		# printing purpose only
		done = done_counter*100/len(pdb_ids)
		sys.stdout.write("Predicting: {}%\r".format(done))
		sys.stdout.flush()

		# extract features of every chain of the given ids, with given parameters
		residues, X , y = prot_dataset.generate_blind_test(
			pdb_id,
			short_win = short_win, 
			large_win = large_win, 
			contact_threshold = contact_threshold
			)

		# build a dataframe
		df = prot_dataset.as_dataframe(X, y)
		# predict using given features
		predictions = clf.predict_proba(df[features])
		# classifier returns a list of pairs [a, b] where 'a' is the 
		# probability of being class '0' and 'b' the probability of being class '1'
		predictions = [pred[1] for pred in predictions]

		# if blur option is enabled perform it
		if blur :
			predictions = blur_by_chain(residues, predictions)

		# print predictions to file
		out_predictions(pdb_id, residues, predictions, path = path)
		
		# printing purpose only
		done_counter += 1
	# printing purpose only
	sys.stdout.write("Predicting: 100%")
	sys.stdout.flush()
	print

# test small window sizes on a range
def find_best_small_window():
	results_list = []
	for sw in range(2,13):
		results_list.append(test(DecisionTreeClassifier(max_depth = 7, random_state = 1), "sw_{}_test".format(sw), sw = sw))

	print "###############################################################################################"
	print "FINAL RESULTS"
	it = 2
	for results in results_list:
		print 
		print "Small Window Size {}".format(it)
		print results
		it += 1

# test large window sizes on a range
def find_best_large_window():
	results_list = []
	for lw in [15,20,25,30,35,40,45,50]:
		results_list.append(test(RandomForestClassifier(max_depth = 10, random_state = 1), "lw_{}_test".format(lw), lw = lw))

	print "###############################################################################################"
	print "FINAL RESULTS"
	it = 15
	for results in results_list:
		print 
		print "Large Window Size {}".format(it)
		print results
		it += 5

# test constact threshold length on a range
def find_best_ct_threshold():
	results_list = []
	for ct in range(2,11):
		results_list.append(test(DecisionTreeClassifier(max_depth = 7, random_state = 1), "ct_{}_test".format(ct), ct = ct))

	print "###############################################################################################"
	print "FINAL RESULTS"
	it = 2
	for results in results_list:
		print 
		print "Contact Threshold Length {}".format(it)
		print results
		it += 1

def find_best_max_depth():
	results_list = []
	for md in range(4,11):
		results_list.append(test(DecisionTreeClassifier(max_depth = md, random_state = 1), "md_{}_test".format(md)))

	print "###############################################################################################"
	print "FINAL RESULTS"
	it = 4
	for results in results_list:
		print 
		print "Max Depth {}".format(it)
		print results
		it += 1

def test_features():
	features = [
	"IrIa_CC", 
	"Intra",
	"Inter",
	"S_Dist", 
	"S_Ang", 
	"S_Ang/Dist",
	"L_Seq_Len", 
	"L_Dist/Seq_Len", 
	"L_Ang*Dist"
	]
	print test(RandomForestClassifier(max_depth = 10, random_state = 1), "fea_test", choosen_features = features)

def test_pdb_ids():
	prot_dataset = ProteinDataset()
	prot_dataset.parse()
	ids_list = prot_dataset.get_prot_list()
	results_list = []
	random.seed(1)
	for test_n in range(1,11):
		sel_list = []
		while len(sel_list) < 20:
			index = random.randrange(len(ids_list))
			if ids_list[index] not in sel_list:
				sel_list.append(ids_list[index])

		results_list.append((str(sel_list), test(RandomForestClassifier(max_depth = 10, random_state = 1), "id_test_{}".format(test_n), tr_prots = sel_list)))

	print "###############################################################################################"
	print "FINAL RESULTS"
	for results in results_list:
		ids, resu = results
		print 
		print " training_set: {}".format(ids)
		print resu
	return 

	
#test(RandomForestClassifier(max_depth = 10, random_state = 1), "default_training_set", tr_prots =["1p22", "1ozs", "2gsi", "1fqj", "1o9a", "1kdx", "1i7w", "1hv2", "1dev", "1tba", "1sc5", "1lm8", "1sb0", "2phe", "1i8h", "1fv1", "1l8c", "2o8a", "2gl7", "1rf8", "1cqt", "2nl9", "1hrt"])
#test_features()																						