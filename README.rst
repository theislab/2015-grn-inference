:Author:
  Andrea Ocone 
  Institute of Computational Biology, Helmholtz Zentrum Muenchen
:Email: <andrea.ocone@helmholtz-muenchen.de>
:Date: April, 13. 2015

SOURCES
=======

Here is a brief description of the main MATLAB and c++ files contained in the
directory in the sequence that needs to be run for producing the toy model network interference as in the paper.

``ToggleSwitch/generationData.m``
  Generates the toy model of differentiation and a sample snapshot data from it based on a desired sampling strategy 

``ToggleSwitch/DiffusionMaps.m``
  Produces the reduced dimension representation of data and plots of diffusion map for the snapshot data by calling the function diffusion_maps_nn.m 

``ToggleSwitch/pseudotime/ordering.m``
  First, assigns each cell in the snapshot data to one (or more if shared) of the branches captured by diffusion maps. It asks the user to manually select the initial cell and the final cell of each branch on the diffusion map plot.

  Second, Wanderlust is performed in this file to get an ordering of cells along each branch. multiple replicates of Wanderlust ordering is the output of this file.

``ToggleSwitch/pseudotime/replicates/generatePseudoTS.m``
  Organises multiple replicates files (from multiple Wanderlust sessions) to be used in the network inference. Then it runs runGPreg.m . 

``ToggleSwitch/pseudotime/replicates/GPs/runGPreg.m``
  Creates the emulators for input genes from the multiple pseudo-time series from the previous step.  

``ToggleSwitch/inference/geneC(or geneD)/ANDpm(or ORpm)/paramest.cpp``
  Uses the definition of the ODE models in the paramestFuns.h file and does parameter estimation with MCMC. 

``ToggleSwitch/inference/geneC(or geneD)/BayesFactor/ANDpm (or ORpm)/computeBF.cpp``
  Computes Bayes factor

PROBLEMS AND COMMENTS
=====================

A more user friendly version of the code will be available soon.
Please address any problem or comment to: <andrea.ocone@helmholtz-muenchen.de>
