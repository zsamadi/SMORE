# SMORE
Spatial Motif Recognition

[data](https://github.com/zsamadi/SMORE/tree/main/data)

A sample input data table with 4 fields, 
1. SID: Section or TIssue ID, which is the animal ID or section ID in case of the cells being from multiples different tissues. This sample dataset has two IDs.
2. cellType: Cell type or cluster ID, this sample dataset has 12 cell types.
3. Centroid_X, Centroid_Y: Spatial coordinates of the cells. Add Centroid_Z if third dimension is available as well.

A sample input gene expression table with 3 fields, 
1. gene_i: The expression for genes. gene names are arbitrary. This sample dataset has expression for three genes, gene_1, gene_2 and gene_3.
   
[run](https://github.com/zsamadi/SMORE/tree/main/run)

This folder contains main execution function for the algorithm. Going into this folder and entering smore without any option in the MATLAB command window will run the algorithm with default settings. The default input data is "sampleDataEb.csv" in "data" folder that contains a random network with motif (A/B)CDE embedded in it. Output results and related figures will be saved in the "output" folder by default. The  settings for the algorithm along with their description and default values are provided in the following table. 
Options |Description | Default 
--- | --- | --- 
__Input\Output__
input | input file | 'sampleData.csv' 
gexFilename | input Gene expression file, can be csv file or h5 file in case of sparse input data | 'sampleDataGEx.csv' 
gAnotFilename| text file containing Gene  names. Only required in case of sparse input gene expression data| '';
output | output folder | 'output' 
__Create Graph__
gMode | The method for creating the graph(options: 'delaunay', 'knn', 'epsilon') | 'delaunay'
ND | Spatial dimensionality of the input data, if it's set to 0, it looks for the dimensionality of the data from input file (options:0, 2, 3) | 0
nNeighs |Number of neighbors for the knn method | 5
rEps |Epsilon radius for epsilon method of creating the graph. If it's set to 0, the radius is computed to be 4 times the average distance from the nearest neighbor. | 0
__Generate Control Data__
shuffleMode | Shuffling method for creating the control data (options: 'shuffle', 'kernel') | 'shuffle'
neighDepth | Number of neighbors for kernel shuffling method | 4
__Sample Input Graph__
samplingFreq | sampling frequency for URPEN method | 1
WMotif | motif length | 4
__SMORE Configuration__
nMotifs |number of output motifs| 5
fixedTypes |Vector of fixed cell types, if set to 0, no cell type is fixed. Single fixed cell types can be entered as integers. For example,  fixedTypes=1 is equivalent to cell type 1 being fixed. fixedTypes=[1, 3] implies that cell types 1 and 3 are fixed | 0
nTrain | number of control data generated for training | 50
nScore |number of control data generated for scoring | 10
isEnrich | Switch for performing enrichment or not (options: true, false) | true
diffMotif |differential enrichment (options: true, false)  | false 
nEval | number of seeds from initial evaluation | 25
nRefine |  number of seeds from initial refinement  | 4
nRefineIter |maximum number of iterations for enrichment | 20
__Gene expression analysis__
doGEA  |switch to perform gene expression analysis (options: true, false)  | true
iSPGEx | sparse gene expression input data (options: true, false)  | false
isGEByTissue | perform gene expression for each SID separately (options: true, false)  | false
gePvalMin  | Minimum pvalue for selected heatmap plots of gene expression  | 0.05
geRandTest  |Rand test for gene expression analysis, useful for finding threshold for significant cases (options: true, false)  | false


You can use the following command to modify all possible smore options: 

```
smore(input='sampleDataEb.csv', gexFilename='sampleDataEb.csv',gAnotFilename='', output='output', WMotif=4, nMotifs=5, fixedTypes=0, nTrain=10, nScore=10, ...
   shuffleMode='shuffle', samplingFreq=1.0, gMode='delaunay', nNeighs=5, rEps=0, neighDepth=4, isEnrich=true, ...
   diffMotif=false, nEval=25, nRefine=4, nRefineIter=20, doGEA=true, iSPGEx=false, isGEByTissue=false,gePvalMin=0.05, geRandTest=false)
```


[func](https://github.com/zsamadi/SMORE/tree/main/func)

Contains related functions for running smore. This folder should be at the same folder as run folder or its path should be added to the MATLAB path. 

[output](https://github.com/zsamadi/SMORE/tree/main/output)

Default folder for saving the output results. The folder currently contains output results of running the algorithm on the default input of "sampleDataEb.csv" file. 

* __GraphWithTypes__ and __GraphWithShuffleTypes__:  The tissue map highlighted with the original cell types and shuffled cell types, respectively. 

* __primShuffleBg__: Backgrounf frequencies for primary and control data. Global shuffling results in the same background frequency. 

* __mLogoi_W_X_Y_Z__: Sequence logo for the motif "i" with consensus seed "WXYZ". 

* __pvaluW4__: Plot of log pvalues of the output motifs. 

* __motifGroupOnGraphiW__: Highlights of the motif "i" of length "W" on the graph tissue. 

* __readme__: The Settings for the current experiment along with the alphabets for the cell types. 

* __outputResults__: Output PWM matrix results along with the seeds that constitute  each motif. 

* __genEx__: The folder that contains gene expression output results. 

* __MTAddr-i-p-c__:  the probability distribution of the delta median with motif dmedian highlighted for the case of (i=motif number, p= position in the motif, c=cell type).

* __volcano__: Volcano plot for log pvalues versus delta median. 

* __geaTable__: table containing gene expression case details, with columns _geneName_: gene name, _motifNumber_: motif number,	_motifPosition_: position in the motif, 	_cellType_: cell type, _pvalue_: log pvalue of the expression difference, _deltaMedian_: delta median,	_motifMedian_: median expression of the cells in the motif, 	_typeMedian_: median expression of the cells of the type, _numCellMotif_: number of cells in the motif, _numCellType_: total number of the cells of the type. 

* __allHeatmap__: Heatmap for all cases. 

* __selectedSortedHeatmap__: Heatmap for selected cases based on the input gePvalMin threshold. 

* __selectedDMedianHeatmap__: Heatmap of delta median values for the selected cases. 

## Questions
Please send questions or possible issues with running the code on your data to Zain Samadi (zainsamadi at ucla.edu). 

## Citing
Detailed description for the algorithm is provided in:


Samadi Z., Askary A.,
"Spatial motifs reveal patterns in cellular architecture of complex tissues", 
BioRxiv, April. 2024. [[full text]](https://www.biorxiv.org/content/10.1101/2024.04.08.588586v1)






