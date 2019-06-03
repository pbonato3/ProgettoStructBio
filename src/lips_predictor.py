#!/usr/bin/env python
from dataset import ProteinDataset
import models
import matplotlib.pyplot as plt
import argparse
from os.path import isfile
import json
import re

def parse_config_file(path = "../config.json"):
    if isfile(path):
        json_obj = {}
        with open(path, 'r') as file:
            json_obj = json.load(file)
            file.close()
        return json_obj

    else:
        print "can't find configuration file."
        return 


def run_from_config(args):
    # parse config file
    params = parse_config_file(args.path)

    #### training set ####

    # if flag is set to true generate a training set to given path
    if params["make-training-set"]:
        # get pdb ids from file
        pdb_ids = params["training-ids"]
        # if the first id is equal to 'all' then use all ids in dataset
        if pdb_ids[0] == "all":
            pdb_ids = prot_dataset.get_prot_list()

        # generate examples
        res, X, y = prot_dataset.generate_random_examples(
            pdb_ids,
            short_win = params["short-window"], 
            large_win = params["large-window"], 
            contact_threshold = params["contact-threshold"], 
            ex_per_chain = params["examples-per-chain"] 
            )
        # if balance flag is set to true
        if params["balance"]:
            # balance number of positive and negative examples
            res, X, y = prot_dataset.balance_neg_pos(res, X, y )
        # output dataset
        prot_dataset.training_set_out(X,y,params["training-set-path"])


    #### model fitting ####

    # if flag is set to true fit the model

    if params["fit-model"]:
        # parse training set
        training_set = prot_dataset.training_set_in(params["training-set-path"])

        # initialize the standard classifier
        predictor = models.make_predictor(
            model_type = params["model-type"], 
            config = params[params["model-type"]], 
            training_set = training_set, 
            features = params["features"]
            )

        models.model_out(predictor, params["trained-model-path"])


    #### prediction ####

    # if flag is set to true predict given ids
    if params["predict"]:

        # load model from file
        predictor = models.model_in(params["trained-model-path"])
        # check if pdb and ring file exists
        prot_dataset.download_pdb(params["predict-ids"])
        prot_dataset.download_ring(params["predict-ids"])

        models.predict(
            clf = predictor, 
            pdb_ids = params["predict-ids"], 
            features = params["features"],
            short_win = params["short-window"], 
            large_win = params["large-window"], 
            contact_threshold = params["contact-threshold"],
            path = params["result-file"],
            blur = params["probability-blur"]
            )

    return



def show_plots(args):
    path = args.path

    df = prot_dataset.training_set_in(path)
    
    for feature in prot_dataset.features_names:
        df.boxplot(column=feature, by="y")
        plt.show()
    
    df.plot.scatter(x="Intra", y="Inter", c="y", colormap='viridis')
    plt.show()
    
    df.plot.scatter(x="L_Dist/Seq_Len", y="L_Ang", c="y", colormap='viridis')
    plt.show()

    df.plot.scatter(x="S_Dist", y="S_Ang", c="y", colormap='viridis')
    plt.show()


def download_pdb_files(args):
    if args.all:
        prot_dataset.download_pdb(prot_dataset.get_prot_list())
    elif args.i :
        prot_dataset.download_pdb(args.i)

    else:
        print "No id given."

def download_ring_files(args):
    if args.all:
        prot_dataset.download_ring(prot_dataset.get_prot_list())
    elif args.i :
        prot_dataset.download_ring(args.i)

    else:
        print "No id given."


def random_dataset(args):
    path = args.path
    win_length = args.win_length
    pdb_ids = args.pdb_ids
    contact_threshold = args.contact_threshold
    n_examples = args.n_examples
    balance = args.balance

    if not path:
        path = "../sets/rand_dataset_w{}c{}.txt".format(win_length, contact_threshold)
    if pdb_ids[0] == "all":
        pdb_ids = prot_dataset.get_prot_list()
    res, X, y = prot_dataset.generate_random_examples(pdb_ids,short_win= win_length, large_win = 60, contact_threshold = contact_threshold, ex_per_chain=n_examples )
    if balance:
        res, X, y = prot_dataset.balance_neg_pos(res, X, y )
    prot_dataset.training_set_out(X,y,path)

def results3D(args):
    import pymol
    from pymol import cmd, util

    path = args.path
    pdb_id = args.pdb_id

    file = open(path,'r')

    pymol.finish_launching()  # Open Pymol (not necessary from pyMOL 2.1)

    cmd.fetch(pdb_id, pdb_id)  # Download the PDB

    cmd.hide("lines", pdb_id)  # Hide lines
    cmd.show("cartoon", pdb_id)  # Show ribbon

    cmd.alter("(all)", "b=0.0")  # Set all B-factor to 0.0

    disp_flag = False
    for line in file:
        if disp_flag:
            if re.match(">", line):
                disp_flag = False
            else:
                identifier, prob, binary_prob = line.split()
                model, chain, index, insertion_code, name = identifier.split('/')

                residue_string = '{}/{}{}/'.format(chain, index, '')  # Residue selection in Pymol syntax
                cmd.alter(residue_string, "b={:2f}".format(float(prob)))  # Substitute the B-score with the LIP score

        if re.match(">{}".format(pdb_id), line):
            print "GOT IT"
            disp_flag = True
            

    cmd.spectrum("b", palette="rainbow", selection="(all)")  # color by B-factor values
    file.close()





def feature3D(args):
    import pymol
    from pymol import cmd, util

    pdb_id = args.pdb_id
    win_length = args.win_length
    contact_threshold = args.contact_threshold
    feature_name = args.feature_name

    residues, features, lip_indexes = prot_dataset.generate_test(pdb_id, short_win= win_length, large_win = 60, contact_threshold= contact_threshold)
    df = prot_dataset.as_dataframe(features,lip_indexes)
    feature = df[feature_name]

    pymol.finish_launching()  # Open Pymol (not necessary from pyMOL 2.1)
    cmd.fetch(pdb_id, pdb_id)  # Download the PDB

    cmd.hide("lines", pdb_id)  # Hide lines
    cmd.show("cartoon", pdb_id)  # Show ribbon

    cmd.alter("(all)", "b=0.0")  # Set all B-factor to 0.0


    for i in range(0,len(residues)):
        residue_string = '{}/{}{}/'.format(residues[i].get_full_id()[2], residues[i].id[1], '')  # Residue selection in Pymol syntax
        cmd.alter(residue_string, "b={:2f}".format(feature[i]))  # Substitute the B-score with the LIP score
        if lip_indexes[i] == 1:
            cmd.show_as("spheres", "chain {} and resi {}".format(residues[i].get_full_id()[2], residues[i].id[1]))

    cmd.spectrum("b", palette="rainbow", selection="(all)")  # color by B-factor values


#######################################################
#################  MAIN PROGRAM   #####################
#######################################################
prot_dataset = ProteinDataset()
prot_dataset.parse()

# create the top-level parser
parser = argparse.ArgumentParser(prog='lips_predictor')
subparsers = parser.add_subparsers(help='sub-command help')

# create the parser for the "plots" command
parser_plots = subparsers.add_parser('plots', help='Show boxplots and scatterplots about the features used in this program')
parser_plots.add_argument('-p','--path', nargs='?', default ="../sets/rand_all_w4d4.5.txt" ,help='Specify the training set to plot. Default is "./rand_all_w4d4.5.txt"')
# set default function
parser_plots.set_defaults(func=show_plots)

# create the parser for the "downpdb" command
parser_down_pdb = subparsers.add_parser('downpdb', help='Download pdb files.')
parser_down_pdb.add_argument('-a','--all', default=False ,  action='store_true', help='Download pdb file for every entry in the dataset.')
parser_down_pdb.add_argument('-i', type= str, metavar='id', nargs='+' , help='Specify the ids to download.')
# set default function
parser_down_pdb.set_defaults(func=download_pdb_files)

# create the parser for the "downring" command
parser_down_ring = subparsers.add_parser('downring', help='Download ring files.')
parser_down_ring.add_argument('-a','--all', default=False ,  action='store_true', help='Download ring file for every entry in the dataset.')
parser_down_ring.add_argument('-i', type= str, metavar='id', nargs='+' , help='Specify the ids to download.')
# set default function
parser_down_ring.set_defaults(func=download_ring_files)

# create the parser for the "show3d" command
parser_show3d= subparsers.add_parser('show3d', help='Show lip residues and a feature in pymol.')
parser_show3d.add_argument('pdb_id', type= str, help='Specify the id to show.')
parser_show3d.add_argument('feature_name', choices= prot_dataset.features_names, help='Feature to show.')
parser_show3d.add_argument('-w','--win_length', type= int, nargs='?', default=4, help='Specify the win_length.')
parser_show3d.add_argument('-c','--contact_threshold', type= float, nargs='?', default=6, help='Specify the contact threshold.')
# set default function
parser_show3d.set_defaults(func=feature3D)



# create the parser for the "results3d" command
parser_results3d= subparsers.add_parser('results3d', help='Show lip residues and a feature in pymol.')
parser_results3d.add_argument('pdb_id', type= str, help='Specify the id to show.')
parser_results3d.add_argument('-p','--path', nargs='?', default = '../results.txt' ,help='Specify path for the output file. Default is "../results.txt"')

# set default function
parser_results3d.set_defaults(func=results3D)


# create the parser for the "rand_data" command
parser_rand_dataset= subparsers.add_parser('rand_data', help='Generate a random dataset taking nex examples from every chain.')
parser_rand_dataset.add_argument('pdb_ids', type= str, nargs='+', help='Ids used to generate examples')
parser_rand_dataset.add_argument('-w','--win_length', type= int, nargs='?', default=4, help='Specify the win_length.')
parser_rand_dataset.add_argument('-c','--contact_threshold', type= float, nargs='?', default=6, help='Specify the contact threshold.')
parser_rand_dataset.add_argument('-x','--n_examples', type= int, nargs='?', default=10, help='Number of examples selected for every chain.')
parser_rand_dataset.add_argument('-b','--balance', default = False, action='store_true', help='Balance number of positive and negative examples.')
parser_rand_dataset.add_argument('-p','--path', nargs='?' ,help='Specify path for the output file. Default is "../sets/training.txt"')

# set default function
parser_rand_dataset.set_defaults(func=random_dataset)

# create the parser for the "rand_data" command
parser_run= subparsers.add_parser('run', help='Run the program with parameters specified in configuration file.')
parser_run.add_argument('-p','--path', type= str, nargs='?', default="../config.json", help='Path to configuration file')

# set default function
parser_run.set_defaults(func=run_from_config)

# parse arguments
args = parser.parse_args()
# call default function giving arguments
args.func(args)




#long_prot =["1cee","1dev","1dow","1fqj","1g3j","1hrt","1i7w","1j2j","1jsu","1l8c","1p4q","1q68","1rf8","1sc5","1sqq","1tba","1th1","1xtg","1zoq","2a6q","2auh","2c1t","2o8a"]
#short_prot = []
#for prot in prot_list:
#    if prot not in long_prot:
#        short_prot.append(prot)
#
#prot_list_test = short_prot[6:-1] + long_prot[6:-1]
#prot_list_tr = short_prot[0:6] + long_prot[0:6]


#residues, X_tr, y_tr= prot_dataset.generate_random_examples(prot_list, 4, 4.5, 10)




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