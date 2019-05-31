#!/usr/bin/env python
from dataset import ProteinDataset
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn import metrics
import matplotlib.pyplot as plt
import copy


def do_predictor():
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

	tr_prots = ['1cee','1dow','1hrt','1i7w','1j2j','1l8c','1rf8','1sqq','3b71','1a3b','1hv2','1ozs','1sb0','1axc','2gl7','1h2k','1ycq','1fv1','1kdx','1cqt']

	res, X, y = prot_dataset.generate_random_examples(tr_prots, 4, 6, 10)
	prot_dataset.training_set_out(X,y)
	df = prot_dataset.as_dataframe(X,y)

	df = prot_dataset.training_set_in()


	clf = DecisionTreeClassifier()

	clf.fit(df[choosen_features],df['y'])

	export_graphviz(clf, out_file="../tree.txt")

	test_prots = []
	for prot_id in prot_dataset.get_prot_list():
		if prot_id not in tr_prots:
			test_prots.append(prot_id)

	avg_pred = []
	avg_y = []
	for t_p in test_prots:
		print "\nPredicting: {}".format(t_p)

		res_test, X_test, y_test = prot_dataset.generate_test(t_p,  6, 4.5)

		df_test = prot_dataset.as_dataframe(X_test, y_test)
		predictions = fill_0_gaps(fill_1_gaps(clf.predict(df_test[choosen_features])))

		avg_y = avg_y + list(y_test)
		avg_pred = avg_pred + list(predictions)

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
		print "\n"
		print "accuracy:  {}".format(metrics.accuracy_score(y_test, predictions))
		print "precision: {}".format(metrics.precision_score(y_test, predictions))
		print "recall:    {}".format(metrics.recall_score(y_test, predictions))

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


do_predictor()