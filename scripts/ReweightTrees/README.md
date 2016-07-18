## Aplying the weights (Andrey)

1. Use *weightTest.py* script first. It will read the input tree form Alexandra and create
large Histograms which map the (file number, event number) to a single weight.  Such
histogram is created for each weight in the input tree.  
2. Run *addWeights.py*. It will insert the weights into the Limits trees. The output is a
new file which contains weights in a separate leaves of a tree.
3. Then, you would want to *hadd* all new files together - this is your new input for the
Limit code.  The limit code itself of course needs to be modified in order to use those
weights.

PS. For script (1), I think a better way to do this would be to use TTree::BuildIndex
method. So, instead of creating histograms one can have a tree indexed by (file,evt)
combination. I have not tried to implement it though, not sure if it will work.
https://root.cern.ch/doc/master/classTTree.html#a3f6b5bb591ff7a5bd0b06eea6c12b998

PPS. The reason for creating TH1 histograam is to not have to loop over the tree every
time for each event - that would take ages. Geting BinContent() of a histogram is much
faster.


## Producing the weights (by Alexandra)
The input files in lxplus are in
inputLM = "/afs/cern.ch/work/z/zghiche/public/ForXanda/NRDir_LM_350/"
inputHM = "/afs/cern.ch/work/z/zghiche/public/ForXanda/NRDir_HM_350/"

We have 13 samples, to each one we have two limit-trees, one for each category (LOWMASS/HIGHMASS). 
We do the new samples in 3 steps

1) read all the 2D histos from root files where they are saved, save them in multidimensional matrices

2) pass over all the events, make a matrix with the same binning of the input. (this is fast)
All the events here are added to the 2D matrix. NO weighting at this point. 

  => We have now  64190.0 in the low mass category and 138303.0 in the high mass one

3) pass over all the events, event-by-event I find in each bin it is, and make the weights dividing the matrix of point (1), by the ones from point (2) and fill new limit-trees for the new samples. 
They will be saved in 3 folders:

lambdaonly 
V3outliers
V3benchmarks

I hope the names inside are self explanatory
In the "V3benchmarks" folder are created as well a file 

V3_LT_output_GluGluToHHTo2B2G_box_validation0_13TeV-madgraph_LOWMASS.root (HIGHMASS)

This is just the sample 0 of v1 (box-only) made with the steps 1,2,3. To we test against the real sample. 
