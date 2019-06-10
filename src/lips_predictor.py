#!/usr/bin/env python
from dataset import ProteinDataset
import models
import argparse
from os.path import isfile
from os import getcwd
import json
import re

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


# parse the configuration file
def parse_config_file(path = main_folder_path+"config.json"):
    # if config file exists
    if isfile(path):
        # initialize a collection
        json_obj = {}
        # load the collection using json
        with open(path, 'r') as file:
            json_obj = json.load(file)
            file.close()
        # return the collection
        return json_obj

    # else print error
    print "can't find configuration file."
    return 


# Run the program from configuration file
def run_from_config(args):
    # parse config file
    params = parse_config_file(main_folder_path+args.path)

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
            res, X, y = prot_dataset.balance_neg_pos(res, X, y , params["positive-lb"])
        # output dataset
        prot_dataset.training_set_out(X,y,main_folder_path+params["training-set-path"])


    #### model fitting ####

    # if flag is set to true fit the model

    if params["fit-model"]:
        # parse training set
        training_set = prot_dataset.training_set_in(main_folder_path+params["training-set-path"])

        # initialize the standard classifier
        predictor = models.make_predictor(
            model_type = params["model-type"], 
            config = params[params["model-type"]], 
            training_set = training_set, 
            features = params["features"]
            )

        models.model_out(predictor, main_folder_path+params["trained-model-path"])


    #### prediction ####

    # if flag is set to true predict given ids
    if params["predict"]:

        # load model from file
        predictor = models.model_in(main_folder_path+params["trained-model-path"])
        pdb_ids = params["predict-ids"]
        # if specified in command, sobstiture
        if args.pdb_ids:
            pdb_ids = args.pdb_ids
        # if the first element is 'all'
        if pdb_ids[0] == "all":
            # clear list
            pdb_ids = []
            # append every id that is not in the training set
            for pdb in prot_dataset.get_prot_list() :
                if pdb not in params["training-ids"]:
                    pdb_ids.append(pdb)


        # check if pdb and ring file exists
        prot_dataset.download_pdb(pdb_ids)
        prot_dataset.download_ring(pdb_ids)

        # run predict command with given parameters
        models.predict(
            clf = predictor, 
            pdb_ids = pdb_ids, 
            features = params["features"],
            short_win = params["short-window"], 
            large_win = params["large-window"], 
            contact_threshold = params["contact-threshold"],
            path = main_folder_path+params["result-file"],
            blur = params["probability-blur"],
            blur_w = params["probability-blur-len"]
            )

    return


# plots features of a given training set comparing lips and not lips
def show_plots(args):
    import matplotlib.pyplot as plt
    # geth path
    path = main_folder_path+args.path

    # if given path doesn't exists print an error
    if not isfile(path):
        print "No file found."
        return 

    # else parse the file as a dataframe
    df = prot_dataset.training_set_in(path)
    # for each feature extracted
    for feature in prot_dataset.features_names:
        # show a boxplot of the feature values seprated by lip flag 
        df.boxplot(column=feature, by="y")
        plt.show()
    # show a scatterplot of inter vs intra chain contacts colored by lip flag
    df.plot.scatter(x="Intra", y="Inter", c="y", colormap='viridis')
    plt.show()
    
    # show a scatterplot of large window mean length vs large window angle colored by lip flag
    df.plot.scatter(x="L_Dist/Seq_Len", y="L_Ang", c="y", colormap='viridis')
    plt.show()

    # show a scatterplot of short window length vs short window angle colored by lip flag
    df.plot.scatter(x="S_Dist", y="S_Ang", c="y", colormap='viridis')
    plt.show()
    return 

# download pdb files by id
def download_pdb_files(args):
    # if all flag is set it will download all pdb in dataset
    if args.all:
        prot_dataset.download_pdb(prot_dataset.get_prot_list())
    # else if ids are given it will download them
    elif args.i :
        prot_dataset.download_pdb(args.i)
    else:
        print "No id given."

# download ring files by id
def download_ring_files(args):
    # if all flag is set it will download all pdb in dataset
    if args.all:
        prot_dataset.download_ring(prot_dataset.get_prot_list())
    # else if ids are given it will download them
    elif args.i :
        prot_dataset.download_ring(args.i)
    else:
        print "No id given."


# show a pdb structure colored with lips probability
def results3D(args):
    import pymol
    from pymol import cmd, util

    # get path adn pdb id
    path = main_folder_path+args.path
    pdb_id = args.pdb_id

    # check if result file exists else print error
    if not isfile(path):
        print "File \"{}\" has not been found.".format(path)
        return

    # open file
    file = open(path,'r')
    pymol.finish_launching()  # Open Pymol (not necessary from pyMOL 2.1)
    # fetch pdb id
    cmd.fetch(pdb_id, pdb_id)  # Download the PDB
    # Hide lines
    cmd.hide("lines", pdb_id)
    # Show ribbon
    cmd.show("cartoon", pdb_id)
    # Set all B-factor to 0.0
    cmd.alter("(all)", "b=0.0")

    # multiple pdb ids can be found in results files, set this to true when desired id is found
    disp_flag = False
    # for each line in the file
    for line in file:
        # if reading the desired id
        if disp_flag:
            # if encountering another id start line set flag to false again and exit
            if re.match(">", line):
                disp_flag = False
                break
            # else it is a line of desired pdb id
            else:
                # split the line
                identifier, prob, binary_prob = line.split()
                model, chain, index, insertion_code, name = identifier.split('/')
                # Residue selection in Pymol syntax
                residue_string = '{}/{}{}/'.format(chain, index, '')
                 # Substitute the B-score with the LIP score
                cmd.alter(residue_string, "b={:2f}".format(float(prob)))

        # if the line was the right start of desired id set flag to true
        if re.match(">{}".format(pdb_id), line):
            disp_flag = True
            
    # color all residues by their B-factor
    cmd.spectrum("b", palette="rainbow", selection="(all)")
    # close file
    file.close()
    return 




# show values of a feature in pymol
def feature3D(args):
    import pymol
    from pymol import cmd, util

    # get parameters
    pdb_id = args.pdb_id
    short_win = args.short_win
    large_win = args.large_win
    contact_threshold = args.contact_threshold
    feature_name = args.feature_name

    # compute features for each residue in a labelled chain
    residues, features, lip_indexes = prot_dataset.generate_test(pdb_id, short_win= short_win, large_win = large_win, contact_threshold= contact_threshold)
    # make a dataframe
    df = prot_dataset.as_dataframe(features,lip_indexes)
    # get the array with desired feature
    feature = df[feature_name]

    # Open Pymol (not necessary from pyMOL 2.1)
    pymol.finish_launching()
    # Download the PDB
    cmd.fetch(pdb_id, pdb_id)
    # Hide lines
    cmd.hide("lines", pdb_id)
    # Show ribbon
    cmd.show("cartoon", pdb_id)
    # Set all B-factor to 0.0
    cmd.alter("(all)", "b=0.0")

    # for each residue 
    for i in range(0,len(residues)):
        # Residue selection in Pymol syntax
        residue_string = '{}/{}{}/'.format(residues[i].get_full_id()[2], residues[i].id[1], '')
        # Substitute the B-score with the feature score
        cmd.alter(residue_string, "b={:2f}".format(feature[i]))
        # If the residue is labelled as LIP show it as sphere
        if lip_indexes[i] == 1:
            cmd.show_as("spheres", "chain {} and resi {}".format(residues[i].get_full_id()[2], residues[i].id[1]))

    # color by B-factor values
    cmd.spectrum("b", palette="rainbow", selection="(all)")


#######################################################
#################  MAIN PROGRAM   #####################
#######################################################
# parse the dataset
prot_dataset = ProteinDataset()
prot_dataset.parse()

# create the top-level parser
parser = argparse.ArgumentParser(prog='lips_predictor')
subparsers = parser.add_subparsers(help='sub-command help')

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

# create the parser for the "plots" command
parser_plots = subparsers.add_parser('plots', help='Shows boxplots and scatterplots about the features used in this program, compraing LIPs and not LIPs examples.')
parser_plots.add_argument('-p','--path', nargs='?', default = "sets/plots.txt" ,help='Specify the training set to plot. Default is "../sets/plots.txt"')
# set default function
parser_plots.set_defaults(func=show_plots)

# create the parser for the "show3d" command
parser_show3d= subparsers.add_parser('show3d', help='Shows the pdb structure in pymol, coloring residues by values of a feature and representing LIP residues as spheres.')
parser_show3d.add_argument('pdb_id', type= str, help='The pdb id that pymol has to fetch.')
parser_show3d.add_argument('feature_name', choices= prot_dataset.features_names, help='Name of the feature to show.')
parser_show3d.add_argument('-sw','--short_win', type= int, nargs='?', default=4, help='Length of short window to be used in feature extraction. Default is 4.')
parser_show3d.add_argument('-lw','--large_win', type= int, nargs='?', default=30, help='Length of large window to be used in feature extraction. Default is 30.')
parser_show3d.add_argument('-c','--contact_threshold', type= float, nargs='?', default=5, help='Contact threshold to be used in feature extraction. Default is 5.')
# set default function
parser_show3d.set_defaults(func=feature3D)

# create the parser for the "results3d" command
parser_results3d= subparsers.add_parser('results3d', help='Shows the pdb structure in pymol, coloring residues with probability of being LIPs Blue (not LIP) to Red (is LIP).')
parser_results3d.add_argument('pdb_id', type= str, help='Specify the id to show.')
parser_results3d.add_argument('-p','--path', nargs='?', default = 'results.txt' ,help='Specify path of the results file that contains the predictions of the desired pdb id. Default is "../results.txt"')
# set default function
parser_results3d.set_defaults(func=results3D)


# create the parser for the "rand_data" command
parser_run= subparsers.add_parser('run', help='Run the program with parameters specified in configuration file.')
parser_run.add_argument('-p','--path', type= str, nargs='?', default="config.json", help='Path to configuration file. Default is "../config.json"')
parser_run.add_argument('-i', '--pdb_ids', type= str, metavar='id', nargs='+' , help='Specify the ids to predict, overriding configuration file. If the first id is "all" every id in dataset is computed.')
# set default function
parser_run.set_defaults(func=run_from_config)

# parse arguments
args = parser.parse_args()
# call default function giving arguments
args.func(args)
