#!/usr/bin/env python
from dataset import ProteinDataset
import matplotlib.pyplot as plt
import argparse

def show_plots(args):
    df = prot_dataset.training_set_in(args.path)
    
    for feature in prot_dataset.features_names:
        df.boxplot(column=feature, by="Lip_Flag")
        plt.show()
    
    df.plot.scatter(x="Intra_CC", y="Inter_CC", c="Lip_Flag", colormap='viridis')
    plt.show()
    
    df.plot.scatter(x="Chain_Dist/Seq_Len", y="Chain_Angle", c="Lip_Flag", colormap='viridis')
    plt.show()

    df.plot.scatter(x="Dist", y="Angle", c="Lip_Flag", colormap='viridis')
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
    if args.pdb_ids[0] == "all":
        args.pdb_ids = prot_dataset.get_prot_list()
    res, X, y = prot_dataset.generate_random_examples(args.pdb_ids,win_length= args.win_length,contact_threshold = args.contact_threshold, ex_per_chain=args.n_examples )
    if args.balance:
        res, X, y = prot_dataset.balance_neg_pos(res, X, y )
    prot_dataset.training_set_out(X,y,args.path)

def feature3D(args):
    import pymol
    from pymol import cmd, util

    residues, features, lip_indexes = prot_dataset.generate_test(args.pdb_id, args.win_length, args.contact_threshold)
    df = prot_dataset.as_dataframe(features,lip_indexes)
    feature = df[args.feature_name]

    pymol.finish_launching()  # Open Pymol (not necessary from pyMOL 2.1)
    cmd.fetch(args.pdb_id, args.pdb_id)  # Download the PDB

    cmd.hide("lines", args.pdb_id)  # Hide lines
    cmd.show("cartoon", args.pdb_id)  # Show ribbon

    cmd.alter("(all)", "b=0.0")  # Set all B-factor to 0.0

    pymol.finish_launching()  # Open Pymol (not necessary from pyMOL 2.1)


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
parser_plots.add_argument('-p','--path', nargs='?', default ="./rand_all_w4d4.5.txt" ,help='Specify the training set to plot. Default is "./rand_all_w4d4.5.txt"')
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

# create the parser for the "rand_data" command
parser_rand_dataset= subparsers.add_parser('rand_data', help='Generate a random dataset taking nex examples from every chain.')
parser_rand_dataset.add_argument('pdb_ids', type= str, nargs='+', help='Ids used to generate examples')
parser_rand_dataset.add_argument('-w','--win_length', type= int, nargs='?', default=4, help='Specify the win_length.')
parser_rand_dataset.add_argument('-c','--contact_threshold', type= float, nargs='?', default=6, help='Specify the contact threshold.')
parser_rand_dataset.add_argument('-x','--n_examples', type= int, nargs='?', default=10, help='Number of examples selected for every chain.')
parser_rand_dataset.add_argument('-b','--balance', default = False, action='store_true', help='Balance number of positive and negative examples.')
parser_rand_dataset.add_argument('-p','--path', nargs='?', default ="./training.txt" ,help='Specify path for the output file. Default is "./training.txt"')

# set default function
parser_rand_dataset.set_defaults(func=random_dataset)

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