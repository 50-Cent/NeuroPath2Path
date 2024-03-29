
_________________________ NeuroPath2Path ________________________

1. Set the Directory of data.

Open neuroClassifyInterclass.m (similarly neuroClassifyIntraclass.m for itraclass classification). 
Set the newpath as the path where the directory containing the dataset is present.


2. Run jobClassifyInter.slurm (similarly jobClassifyIntra.slurm). 
The default parameter values are given. 
User can choose different parameters. 
To do that, open the .slurm file and change the arguments of neuroClassifyInterclass(). 
Ex: neuroClassifyInterclass(0.2,15,5,3,9) indicates that this instance of experiment tests on the 20% perclass-randomly sampled data. The binwidth of perclass-path distribution is 15. For a test neuron, the closest candidates are searched through step length 5 (search space reduction) and continue till +/- (9*5) paths from the number of paths of the test neuron. Majority vote isperformed on 3 nearest neighbors.   

For details, please read the paper.

3. Once the job is completed, open OUTPUT.log



% File objectives

% classification
readinput         :  Read swc files from directory
computenP2P         :  Compute path features
computeDSTInter & Intra     : Path alignment + distance between paths
munkres + distanceNeuronNew : Distance between a pair of neurons and path 			   correspondence
jobClassifyInter and jobClassifyIntra : job script files to run in cluster   
neuroClassifyInterclass and neuroClassifyIntraclass: main files


% importance weights
jobIW : jobscript to get the importance weights
getFeatureWeights : Obtain the importance weights 
