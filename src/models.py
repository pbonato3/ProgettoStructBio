#!/usr/bin/env python
from dataset import ProteinDataset
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn import metrics
import matplotlib.pyplot as plt
import copy
from os.path import isfile

aa_3to1 = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
}

def test_1():
	prot_dataset = ProteinDataset()
	prot_dataset.parse()
	choosen_features = [
		"IrIa_CC", 
		"Intra",
		"Inter",
		"Dist", 
		"Ang", 
		"Ang/Dist",
		"Ch_Dist/Seq_Len", 
		]

	tr_prots = ['2ivz','1dow','1hrt','1i7w','1j2j','1l8c','1rf8','1sqq','3b71','1a3b','1hv2','1ozs','1i8h','1axc','2gl7','1h2k','1ycq','1fv1','1kdx','1cqt']

	#res, X, y = prot_dataset.generate_random_examples(tr_prots, 4, 5, 100)
	#prot_dataset.training_set_out(X,y)
	#df = prot_dataset.as_dataframe(X,y)

	df = prot_dataset.training_set_in()


	clf = DecisionTreeClassifier(max_depth = 7)

	clf.fit(df[choosen_features],df['y'])

	export_graphviz(clf, out_file="../tree.txt")

	test_prots = []
	for prot_id in prot_dataset.get_prot_list():
		if prot_id not in tr_prots:
			test_prots.append(prot_id)

	avg_pred = []
	avg_y = []

	clear_result_file()

	for t_p in test_prots:
		print "\nPredicting: {}".format(t_p)

		res_test, X_test, y_test = prot_dataset.generate_test(t_p,  4, 5)

		df_test = prot_dataset.as_dataframe(X_test, y_test)
		#predictions = fill_0_gaps(fill_1_gaps(clf.predict(df_test[choosen_features])))
		predictions = clf.predict_proba(df_test[choosen_features])
		#predictions = clf.predict(df_test[choosen_features])
		predictions = [pred[1] for pred in predictions]
		bin_pred = prob_to_binary(prob_blur(predictions))

		out_predictions(t_p, res_test, prob_blur(predictions))

		avg_y = avg_y + list(y_test)
		avg_pred = avg_pred + list(bin_pred)


		#print "#############################"
		#print "           TRUTH"
		#print "#############################"
		#print y_test
		#
		#print "#############################"
		#print "         PREDICTIONS"
		#print "#############################"
		#print predictions
		print "\n"
		print "accuracy:  {}".format(metrics.accuracy_score(y_test, bin_pred))
		print "precision: {}".format(metrics.precision_score(y_test, bin_pred))
		print "recall:    {}".format(metrics.recall_score(y_test, bin_pred))

	print "\n"
	print "accuracy:  {}".format(metrics.accuracy_score(avg_y, avg_pred))
	print "precision: {}".format(metrics.precision_score(avg_y, avg_pred))
	print "recall:    {}".format(metrics.recall_score(avg_y, avg_pred))



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


def prob_blur(pred, win_len = 6):
	new = []
	for idx in range(0, len(pred)):
		start = max(idx - win_len, 0)
		stop = min(idx + win_len, len(pred)-1)
		prob_count = 0.0
		for label in pred[start:stop+1]:
			prob_count += label
		new.append(prob_count/(stop - start + 1))

	return new	

def prob_to_binary(pred):
	new = []
	for pr in pred:
		if pr >= 0.5:
			new.append(1)
		else:
			new.append(0)
	return new

def clear_result_file(path = "../results.txt"):
	if isfile(path):
		file = open(path, 'w')
		file.close()


def out_predictions(p_id, residues, predictions, path = "../results.txt"):
	bin_pred = prob_to_binary(predictions)
	file = open(path, 'a')
	file.write(">{}\n".format(p_id))
	for idx in range(0,len(residues)):
		res_id = residues[idx].get_full_id()
		insertion_code = res_id[3][2]
		if insertion_code == ' ':
			insertion_code = ''
		file.write("{}/{}/{}/{}/{} {:.3f} {}\n".format(res_id[1],res_id[2],res_id[3][1],insertion_code,aa_3to1[residues[idx].get_resname()],predictions[idx],bin_pred[idx]))

	file.close()


test_1()