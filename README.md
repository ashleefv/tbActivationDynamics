# tbActivationDynamics
tbActivationDynamics: open-source Python code for a mathematical model of tuberculosis granuloma activation

A computational model that aims to capture the effect of MMP-1 on tuberculosis (TB) disease outcome including initial activation (uncontrolled infection), latent infection (steady state), and the transition from latent to active TB after biological perturbations.

[![DOI](https://zenodo.org/badge/102527971.svg)](https://zenodo.org/badge/latestdoi/102527971) 

## Overview
This model attempts to model the reactivation of a tuberculosis granuloma.
The model is built on top of a model of the initial immune response by Sud, D.; Bigbee, C.; Flynn, J.L.; Kirschner, D.E. Contribution of CD8+ T cells to control of Mycobacterium tuberculosis infection. J. Immunol. 2008, 176, 4296-4314.
This model includes the upregulation of MMP-1 inside of a tuberculosis granuloma, the breakdown of collagen catalyzed by MMP-1,
and the movement of tuberculosis bacteria through the deteriorated granuloma.


## TB Granuloma Activation Mathematical Model
### Code Authors
Steve Ruggiero, Ashlee N. Ford Versypt, 
School of Chemical Engineering,
Oklahoma State University.
Corresponding author: A. N. Ford Versypt, ashleefv@okstate.edu

### Related publication for model details
Steve Ruggiero*, Minu Pilvankar*, Ashlee N. Ford Versypt, Mathematical Modeling of Tuberculosis Granuloma Activation, Processes, 5, 79, 2017. DOI: 10.3390/pr5040079

### Main files

* ModelTestScript.py
   Runs the model using the parameters in ColDet-model-parameters.csv and plots the results. 
   The graphs are saved to jpeg files and are displayed when run using an ipython console.
	
* GeneDeletion.py
   Runs the model with a single gene deleted. Each of the genes deleted in the script
   are used to plot the different outcomes that can occur when a gene is deleted.
	
* k_Cgraphs.py
   Runs the model with varying values of k_C. This script produces graphs showing
   the effect of varying k_C on disease progression.
	
* LocalSensitivityAnalysis.py
   Performs a local sensativity analysis on the model. This script tweaks each parameter's
   value, runs the model in parallel, and records the ratio between the percent change in 
   bacterial leakage and the percent change of the parameter.
	
* LocalSensitivityBarGraph.py
   Takes the pickled output from LocalSensitivityAnalysis.py and produces a bar graph from the results.
	
* TcellDepletion.py
   Runs the model for a brief period, and then gradually reduces parameters associated
   with the production of T cells. The script also produces graphs of the results of the model.

### Dependent & supplemental files

* ColDetModel.py
	This file contains the class containing the developed model.
	
* ColDet-model-parameters.csv
	Contains the model parameters that are read by ImportParamFromFile.py.
	Changes to this file change the parameters used in model runs.
	
* GeneDeletionPaper.py
	A version of the GeneDeletion.py script, but run with 3 specific genes deleted.
	This reproduces the graph included in the paper.

* ImportParamFromFile.py
	Contains a function used to parse parameters from a .csv file.
	
* LocalSensitivity.p
	Pickled python variable containing the results of the local sensitivity analysis.
	Used in the LocalSensitivityBarGraph.py script.
	
* RunModel.py
	Initializes parameters and runs the model.
	
* SensativityAnalysis.csv
	Spreadsheet containg the results of the local sensativity analysis.
	
* SimpleParallelFuncs.py
	This file contains functions written to make use of the Multiprocessing Python
	package easier to use.

   
## Example Usage
Run one of the scripts under main files to produce model output. 
To edit the parameters of the model, open ColDet-model-parameters.csv in Excel or an other spreadsheet or text editor.
