# SMORE
Spatial Motif Recognition

[data](https://github.com/zsamadi/SMORE/tree/main/data)

A sample input data table with 7 fields, 
1. NID: node ID, which is the animal ID or section ID in case of the cells being from multiples different tissues. This sample dataset has to IDs.
2. nodeType: cell type or cluster ID, this sample dataset has 12 cell types.
3. Centroid_X, Centroid_Y: spatial coordinates of the cells. Add Centroid_Z if third dimension is available as well.
4. gene_i: the expression for genes. gene names are arbitrary as long as they're on the columns after spatial coordinates. This sample dataset has expression for three genes, gene_1, gene_2 and gene_3.
   
[run](https://github.com/zsamadi/SMORE/tree/main/run)

This folder contains main execution function for the algorithm. Entering smore without any option in the MATLAB command window will run the algorithm with default settings. The default input data is "sampleData.csv" in "data" folder that contains a random network with motif (A/B)CDE embedded in it. Output results and related figures will be saved in the "output" folder by default. The  settings for the algorithm along with their description and default values are provided in the following table. 
Options |Description | Default 
--- | --- | --- 
__Input\Output__
input | input file | 'sampleData.csv' 
output | output folder | 'output' 
__Create Graph__
gMode | The method for creating the graph(options: 'delaunay', 'knn', 'epsilon') | 'delaunay'
ND | Spatial dimensionality of the input data, if it's set to 0, it looks for the dimensionality of the data from input file (options:0, 2, 3) | 0
nNeighs |Number of neighbors for the knn method | 5
rEps |Epsilon radius for epsilon method of creating the graph| 50
__Generate Control Data__
shuffleMode | Shuffling method for creating the control data (options: 'shuffle', 'kernel') | 'shuffle'
neighDepth | Number of neighbors for kernel shuffling method | 4
__Sample Input Graph__
samplingFreq | sampling frequency for URPEN method | 1
WMotif | motif length | 4
__SMORE Configuration__
nMotifs |number of output motifs, overrides threshold settings.| 5
fixedTypes |Vector of fixed cell types, if set to 0, no cell type is fixed. Single fixed cell types can be entered as integers. For example,  fixedTypes=1 is equivalent to cell type 1 being fixed. fixedTypes=[1, 3] implies that cell types 1 and 3 are fixed | 0
nTrain | number of control data generated for training | 50
nScore |number of control data generated for scoring | 1
isEnrich | Switch for performing enrichment or not (options: true, false) | true
diffMotif |differential enrichment (options: true, false)  | false 
nEval | number of seeds from initial evaluation | 25
nRefine |  number of seeds from initial refinement  | 4
nRefineIter |maximum number of iterations for enrichment | 20

you can use following command to modify smore options. 


`smore(input='sampleDataEb.csv', output='output', WMotif=4, nMotifs=5, fixedTypes=0, nTrain=50, nScore=1, shuffleMode='shuffle', samplingFreq=1.0, gMode='delaunay', nNeighs=5, rEps=0, neighDepth=4, isEnrich=true, diffMotif=false, nEval=25, nRefine=4, nRefineIter=20)`


[func](https://github.com/zsamadi/SMORE/tree/main/func)

Contains related functions for running smore. This folder should be at the same folder as run folder or its path should be added to the MATLAB path. 

[output](https://github.com/zsamadi/SMORE/tree/main/output)

Default folder for saving the output results. The folder currently contains output results of running the algorithm on the default input of "sampleData.csv" file. 






