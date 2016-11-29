##
## Configs for effs plotter
##

stepLegs = {
0:"Total",
1:">= 2#gamma + 2j", #Two photons and two jets requirement
2:"Trigger", #Trigger requirement
3:"#gamma's Kin-Sel", #Photons kinematic selection
4:"#gamma's Selection", #Photons ID selection
5:">= 2 jets", #Two jets again, does nothing
6:"Jet sel 1", #Kinematic preselection on jets
7:"Jets Selection" #Dijet selection
}
steps = [1, 2, 4, 7]

#For Resonant
'''
outname = "res_effs.pdf"
xaxis = [
250, 260, 270, 280, 300, 320, 340, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900
]
xtitle = "M_{X} [GeV]"
ffolder = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/dev-rafael-Nov4/CMSSW_8_0_20/src/flashgg/bbggTools/test/RunJobs/Signal_Nov24/Hadd/"
files = [
"output_GluGluToRadionToHHTo2B2G_M-250_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-260_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-270_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-280_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-300_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-320_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-340_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-350_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-400_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-450_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-500_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-550_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-600_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-650_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-700_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-750_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-800_narrow_13TeV-madgraph.root",
"output_GluGluToRadionToHHTo2B2G_M-900_narrow_13TeV-madgraph.root"
]
'''

#For NonResonant
outname = "nonres_effs.pdf"
xaxis = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
xtitle = "Shape Benchmark Points"
ffolder = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/dev-rafael-Nov4/CMSSW_8_0_20/src/flashgg/bbggTools/test/RunJobs/Signal_Nov24/Hadd/"
files = [
"output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root",
"output_GluGluToHHTo2B2G_node_box_13TeV-madgraph.root",
"output_GluGluToHHTo2B2G_node_2_13TeV-madgraph.root",
"output_GluGluToHHTo2B2G_node_3_13TeV-madgraph.root",
"output_GluGluToHHTo2B2G_node_4_13TeV-madgraph.root",
"output_GluGluToHHTo2B2G_node_5_13TeV-madgraph.root",
"output_GluGluToHHTo2B2G_node_6_13TeV-madgraph.root",
"output_GluGluToHHTo2B2G_node_7_13TeV-madgraph.root",
"output_GluGluToHHTo2B2G_node_8_13TeV-madgraph.root",
"output_GluGluToHHTo2B2G_node_9_13TeV-madgraph.root",
"output_GluGluToHHTo2B2G_node_10_13TeV-madgraph.root",
"output_GluGluToHHTo2B2G_node_11_13TeV-madgraph.root",
"output_GluGluToHHTo2B2G_node_12_13TeV-madgraph.root",
"output_GluGluToHHTo2B2G_node_13_13TeV-madgraph.root",
]

extraSteps = [

"genTotalWeight*((leadingJet_bDis > 0.8 && subleadingJet_bDis < 0.8)||(leadingJet_bDis < 0.8 && subleadingJet_bDis > 0.8))",
"genTotalWeight*(leadingJet_bDis > 0.8 && subleadingJet_bDis > 0.8)"

]

extraStepsN = [8, 9]
extraLegs = { 8:"MPC", 9:"HPC"}
stepLegs.update(extraLegs)

