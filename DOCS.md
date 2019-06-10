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

Features are:

- IrIa_CC
- Intra
- Inter
- L_CC
- S_CC
- f_I
- f_II
- f_III
- f_IV
- f_V
- S_Dist
- S_Ang
- S_Ang/Dist
- L_Seq_Len
- L_Dist/Seq_Len
- L_Ang
- L_Ang/Dist