# Lips Predictor Code Documentation
## Project Structure
In the main folder can be found :
- **src**         folder: containing all the code developed.
- **pdb_files**   folder: where the program looks for existing pdb files.
- **ring_files**  folder: where the program looks for RING files. [About RING](http://protein.bio.unipd.it/ring/).
- **sets**        folder: that is the default location of training sets.
- **img**         folder: containing some images collected during developement.
- **tests**       folder: that is the default location of files used wihile performing tests during developement.

In sets files it is also stored lips_dataset.txt that is the high level dataset used for training.

In **src** folder there are:
* *lips_predictor.py*: the main file that contains the parser and calls functions in other files to execute commands.
* *dataset.py*: the file that contains all usefull functions and classes needed in dataset generation and i/o. To make it more easy to use a class called ProteinDataset holds all the functions used in other files.
* *models.py*: the file that constains all usefull functions needed to train, predict and evaluate models.

## Usage instructions
This program is designed to be used with Python 2 to interact with pymol.

To run the program navigate to src folder and write **python lips_predictor.py -h**

This project uses relative paths but variables are set to find right paths even if running the program from main directory.
When a path is given in input the program appends the main folder location in front of it.
For example to point to *sets* folder the right path is "*sets/a_set.txt*".

The program needs pdb and ring files situated in the correct folder in order to work properly.

There are 5 main commands that can be used:

#### downpdb
Download a given list of pdb structures. Optionally can download the files for each pdb in *lips_dataset.txt*.

#### downring
Download a given list of RING contacts. Optionally can download the files for each pdb in *lips_dataset.txt*.

#### plots
Shows boxplots and scatterplots about the features used in this program, compraing LIPs and not LIPs examples.
It uses a training set file to display the plots. *plots.txt* is given as default but it takes a path as an optional parameter.

#### show3d
Shows the pdb structure in pymol, coloring residues by values of a feature and representing LIP residues as spheres.
Features are computed during execution, so it is possible to set parameters as window sizes and contacts thresholds.
Run with '-h' for more informations on parameters.

#### results3d
Shows the pdb structure in pymol, coloring residues with probability of being LIPs from Blue (not LIP) to Red (is LIP).
Data for the visualization is taken from a result file, default is 'results.txt'

#### run
This is the main command of the program. It takes parameters from a configuration file and performs desired actions.

The file looks like this:
```json
{
	"short-window":4,
	"large-window":30,
	"contact-threshold":5,

	"make-training-set" : false,
	"examples-per-chain" : 100,
	"balance" : false,
	"positive-lb": 40,
	"training-set-path" : "sets/default_training_set.txt",
	"training-ids" : [
		"all"
	],

	"fit-model" : true,
	"model-type": "random-forest",
	"trained-model-path" : "trained-model.sav",	
	"features": [
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
		],

	"predict" : true,
	"predict-ids":[
		"all"
	],
	"probability-blur": true,
	"probability-blur-len": 6,
	"result-file": "results.txt",


	"decision-tree": {
		"max-depth":7,
		"random-state":1
	},

	"random-forest": {
		"n-estimators":10,
		"max-depth":null,
		"random-state":1
	},


	"multi-layer-perceptron": {
		"hidden-layers":[2,2,2,2,2,2],
		"activation":"relu",
		"solver":"adam",
		"alpha":0.0001,
		"learning-rate-init":0.001,
		"max-iter":500,
		"tol":1e-4,
		"random-state":1
	}
}
```
It is divided in 3 differesnt steps that can be turned on/off as pleased:
* Training set generation
* Predictor fitting
* Prediction
*making-training-set*, *fit-model* and *predict* are the switches of theyr respective step of computation.
Every step produces an output file that is used in the next step, so if the file exists the step can be turned off.

```json
	"short-window":4,
	"large-window":30,
	"contact-threshold":5,
```
These are the parameters that manage feature extraction. More infomation in **Features** section below.

Training set is produced selecting a random samples from every chain of specified ids. In config files it is possible to chose the number of examples per chain and if to balance number of positive vs. negative examples. *positive-lb* stands for lower-bound of positives pecent. 
```json
	"examples-per-chain" : 100,
	"balance" : false,
	"positive-lb": 40,
```

The model type can be choosen from the 3 possibilities: *decision-tree*, *random-forest*, *multi-layer-perceptron*. Every model has it's own parameters that can be set in theyr own subsection.
To train a model the desired subset of features can be choosen.

During *prediction step* a set of pdb ids can be given. If 'all' is the first id, all pdb files in *lips_dataset.txt* that are not also in the traing set are predicted.

As final parameters there are boosting options:
```json
	"probability-blur": true,
	"probability-blur-len": 6,
```
If set to true applyes a blur filtering over probability predictions by chain. Blurring uses a *probability-blur-len* x 2 + 1 window.


## Features

Boxplots about discriminance of the features can be seen using command *python lips_predictor.py plots*. 

##### IrIa_CC
Stands for Inter over Intra Chain Contacts. Contacts are computed using RING software and it is possible to configure a threshold. 
This feature is one of the most important to discriminate between LIP and not LIP residues, cause an higher ratio is expected.
It is possible to configure a threshold that is usefull to discard too long contacts.

##### Intra
Stands for Intra-Chain Contacts. Can be usefull to discriminate between LIP and not LIP residues, but seems less important than *IrIa_CC*.
This feature is extracted from shorter window and meaned by number of residues in the window.

##### Inter
Stands for Inter-Chain Contacts. Can be usefull to discriminate between LIP and not LIP residues, but seems less important than *IrIa_CC*.
This feature is extracted from shorter window and meaned by number of residues in the window.

##### L_CC
Stands for Long-Chain Contacts. It is computed considering contacts threshold: if it is longer than half the threshold distance, a contact is considered long. Seems that there is no significative difference between LIP and not LIP residues.
This feature is extracted from shorter window and meaned by number of residues in the window.

##### S_CC
Stands for Short-Chain Contacts. It is computed considering contacts threshold: if it is shorter than half the threshold distance, a contact is considered short. Seems that there is no significative difference between LIP and not LIP residues.
This feature is extracted from shorter window and meaned by number of residues in the window.

##### f_I, f_II, f_III, f_IV, f_V
These five features are extracted from:
> Solving the protein sequence metric problem
> - William R. Atchley, Jieping Zhao, Andrew D. Fernandes, and Tanja Drüke
Every feature is a cluster of chmical-phisical propertyes with high covariance. They don't seems to be so usefult for this task.
Each feature is extracted from shorter window and meaned by number of residues in the window.

##### S_Dist
Stands for Short Window Distance. It represents the distance (in Armstrong) between the frist and last residue in the shorter window.
Togheter with *S_Ang* can be very good to recognise secondary structures, but doesn't seems to discriminate LIP.
This feature is extracted from shorter window and meaned by number of residues in the window.

##### S_Ang
Stands for Short Window Angle. It represents the angle (degrees) between frist, central and last residue in the shorter window.
Togheter with *S_Dist* can be very good to recognise secondary structures, but doesn't seems to discriminate LIP.

##### S_Ang/Dist
Stands for Short Window Angle over Mean Distance. Can be very good to recognise secondary structures, but doesn't seems to discriminate LIP.

##### L_Seq_Len
Stands for Large Window Sequence Length and it is the raw count of residues in the largest window. It is perfectly normal that it variates between window length and half window length (at the beginning of the chain), but when lower it covers the entire chain. In this case from the dataset seems there is a very high probability of being a LIP.

##### L_Dist/Seq_Len
Stands for Large Window Mean Distance over Sequence Length. *L_Dist* is the distance (in Armstrong) between the frist and last residue in the larger window. It is one of the most usefull features for LIP prediction task cause it represent a measure of linearity of the structure.

##### L_Ang
Stands for Large Window Angle. It represents the angle (degrees) between frist, central and last residue in the larger window. It can be considered as another measure of linearity of the structure but it doesn't perform as good as *L_Dist/Seq_Len*.

##### L_Ang*Dist
Stands for Large Window Angle multiplied Mean Distance. It is a union of Large Window linearity measures. It seems to improve a little the score compared to *L_Dist/Seq_Len*.