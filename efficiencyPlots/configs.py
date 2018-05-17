##
## Configs for effs plotter
##

DEBUG=0

outLoc = ""

colors = ["#4A2545", "#389187","#DB995A", "#064789", "#B33951", "#FF1D15","#F0A202"]

stepLegs = {
0:"Total",
1:" 2#gamma + 2j", #Two photons and two jets requirement
2:"Online", #Trigger requirement
3:"#gamma's Kin-Sel", #Photons kinematic selection
4:"Diphoton", #Photons ID selection
5:">= 2 jets", #Two jets again, does nothing
6:"Jet sel 1", #Kinematic preselection on jets
7:"Dijet" #Dijet selection
}
steps = [2, 4, 7]

isRes = 1
isRadion = 0

ffolder = "/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/Mar82018_ForPubli_RafStyle/Signal/Hadd/"
#ffolder = '/tmp/rateixei/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/May2_Mjj70to190_NewCatMVA/EGML_Signal_Mjj70_NewMVA/Hadd/'
#ffolder = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/flashgg_tag-Moriond17-v5/CMSSW_8_0_26_patch1/src/flashgg/bbggTools/test/RunJobs/AllSignal_Feb25_mjj80_pt25/Hadd/"
outname = "res_effs_mjj70_MVAcat_spin2.pdf"

if isRes==1:
  #For Resonant
  header = "#font[61]{pp#rightarrowX#rightarrowHH#rightarrow#gamma#gammab#bar{b} (spin-2)}: selection steps"
  if isRadion: header = header.replace('spin-2', 'spin-0')
  if isRadion: outname = outname.replace('spin2', 'spin0')
  #  250, 260, 270, 280, 300, 350, 400, 450, 500, 600, 700, 800, 900
  xaxis = [
    260, 270, 280, 300, 350, 400, 450, 500, 600, 700, 800, 900
  ]
  xmin = 210
  xmax = 950
  ymax = 1.1
  ymin = 0
  yoff = 0.0
  legymin = 0.65
  xtitle = "m_{X} [GeV]"
#  xtitle = "Resonance Mass [GeV] (spin-2)"
#  if isRadion: xtitle = xtitle.replace('spin-2', 'spin-0')
  drawOpt = "PL"
  #"/afs/cern.ch/work/r/rateixei/work/DiHiggs/dev-rafael-Nov4/CMSSW_8_0_20/src/flashgg/bbggTools/test/RunJobs/Signal_Nov24/Hadd/"
  files_1 = [
#"output_GluGluToBulkGravitonToHHTo2B2G_M-250_narrow_13TeV-madgraph-0.root",
"output_GluGluToBulkGravitonToHHTo2B2G_M-260_narrow_13TeV-madgraph_0.root",
"output_GluGluToBulkGravitonToHHTo2B2G_M-270_narrow_13TeV-madgraph_0.root",
"output_GluGluToBulkGravitonToHHTo2B2G_M-280_narrow_13TeV-madgraph_0.root",
"output_GluGluToBulkGravitonToHHTo2B2G_M-300_narrow_13TeV-madgraph_0.root",
"output_GluGluToBulkGravitonToHHTo2B2G_M-350_narrow_13TeV-madgraph_0.root",
"output_GluGluToBulkGravitonToHHTo2B2G_M-400_narrow_13TeV-madgraph_0.root",
"output_GluGluToBulkGravitonToHHTo2B2G_M-450_narrow_13TeV-madgraph_0.root",
"output_GluGluToBulkGravitonToHHTo2B2G_M-500_narrow_13TeV-madgraph_0.root",
"output_GluGluToBulkGravitonToHHTo2B2G_M-600_narrow_13TeV-madgraph_0.root",
"output_GluGluToBulkGravitonToHHTo2B2G_M-700_narrow_13TeV-madgraph_0.root",
"output_GluGluToBulkGravitonToHHTo2B2G_M-800_narrow_13TeV-madgraph_0.root",
"output_GluGluToBulkGravitonToHHTo2B2G_M-900_narrow_13TeV-madgraph_0.root"
  ]
  files = files_1
  if isRadion:
    files = [ ff.replace('BulkGraviton', 'Radion') for ff in files_1 ]

  extraSteps = [
#"genTotalWeight*( isSignal && (leadingJet_bDis > 0.92 && subleadingJet_bDis > 0.8 && subleadingJet_bDis < 0.92)||(subleadingJet_bDis > 0.92 && leadingJet_bDis > 0.8 && leadingJet_bDis < 0.92))+genTotalWeight*( isSignal && leadingJet_bDis > 0.92 && subleadingJet_bDis > 0.92)",
"(genTotalWeight*(isSignal)*(ResHHTagger_LM > 0.7))",
"(genTotalWeight*(isSignal)*(ResHHTagger_HM > 0.))"
#"genTotalWeight*( isSignal && leadingJet_bDis > 0.92 && subleadingJet_bDis > 0.92)"

  ]

  extraLegs = {
  int(len(stepLegs)+0):"MPC+HPC (LM)",
  int(len(stepLegs)+1):"MPC+HPC (HM)"
  #int(len(stepLegs)+1):"HPC"
  }

else:
  #For NonResonant
#  outname = "nonres_effs_hm.pdf"
  outname = "nonres_effs_mjj70_MVAcat_nonononres.pdf"
  steps = []
  header = "#font[61]{pp#rightarrowHH#rightarrowb#bar{b}#gamma#gamma} Mass Categorization"
  xaxis = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 'SM', '#kappa_{#lambda} = 0']
  xmin = -1
  xmax = 14
  ymin = 0
  ymax = 0.4
  yoff = 0.1
  legymin = 0.75
  xtitle = "Shape Benchmark Points"
  drawOpt = "P"
#  ffolder = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/dev-rafael-Nov4/CMSSW_8_0_20/src/flashgg/bbggTools/test/RunJobs/Signal_Nov24/Hadd/"
  ffolder = "/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/May2_Mjj70to190_NewCatMVA/EGML_Signal_Mjj70_NewMVA/Hadd/"
  files = [
'output_LT_NR_Nodes_All_merged_kl_10p0_kt_1p5_cg_0p0_c2_-1p0_c2g_0p0.root',
'output_LT_NR_Nodes_All_merged_kl_15p0_kt_1p0_cg_0p0_c2_1p0_c2g_0p0.root',
'output_LT_NR_Nodes_All_merged_kl_15p0_kt_1p0_cg_-1p0_c2_0p0_c2g_1p0.root',
'output_LT_NR_Nodes_All_merged_kl_1p0_kt_1p0_cg_0p0_c2_-1p5_c2g_-0p8.root',
'output_LT_NR_Nodes_All_merged_kl_1p0_kt_1p0_cg_-0p6_c2_1p0_c2g_0p6.root',
'output_LT_NR_Nodes_All_merged_kl_1p0_kt_1p0_cg_0p8_c2_0p0_c2g_-1p0.root',
'output_LT_NR_Nodes_All_merged_kl_1p0_kt_1p0_cg_-0p8_c2_0p5_c2g_0p6.root',
'output_LT_NR_Nodes_All_merged_kl_2p4_kt_1p0_cg_0p2_c2_0p0_c2g_-0p2.root',
'output_LT_NR_Nodes_All_merged_kl_2p4_kt_1p0_cg_1p0_c2_0p0_c2g_-1p0.root',
'output_LT_NR_Nodes_All_merged_kl_-3p5_kt_1p5_cg_0p0_c2_-3p0_c2g_0p0.root',
'output_LT_NR_Nodes_All_merged_kl_5p0_kt_1p0_cg_0p2_c2_0p0_c2g_-0p2.root',
'output_LT_NR_Nodes_All_merged_kl_7p5_kt_1p0_cg_0p0_c2_-1p0_c2g_0p0.root',
'output_LT_NR_Nodes_All_merged_kl_1p0_kt_1p0_cg_0p0_c2_0p0_c2g_0p0.root',
'output_LT_NR_Nodes_All_merged_kl_0p0_kt_1p0_cg_0p0_c2_0p0_c2g_0p0.root'
  ]

  extraSteps = [
#"(NR_weight*(isSignal)*(HHTagger_HM > 0.6)*((diHiggsCandidate.M()-dijetCandidate.M()-diphotonCandidate.M()+250)>350))",
#"(NR_weight*(isSignal)*(HHTagger_LM > 0.6)*((diHiggsCandidate.M()-dijetCandidate.M()-diphotonCandidate.M()+250)<350))"
"(NR_weight*(isSignal)*(HHTagger_HM > -10.6)*((diHiggsCandidate.M()-dijetCandidate.M()-diphotonCandidate.M()+250)>350))",
"(NR_weight*(isSignal)*(HHTagger_LM > -10.6)*((diHiggsCandidate.M()-dijetCandidate.M()-diphotonCandidate.M()+250)<350))"
#"genTotalWeight*( isSignal && (diHiggsCandidate.M()-dijetCandidate.M()+125) > 350 )",
#"genTotalWeight*( isSignal && (diHiggsCandidate.M()-dijetCandidate.M()+125) > 350  && (leadingJet_bDis > 0.92 && subleadingJet_bDis > 0.8 && subleadingJet_bDis < 0.92)||(subleadingJet_bDis > 0.92 && leadingJet_bDis > 0.8 && leadingJet_bDis < 0.92))",
#"genTotalWeight*( isSignal && (diHiggsCandidate.M()-dijetCandidate.M()+125) > 350  && leadingJet_bDis > 0.92 && subleadingJet_bDis > 0.92)"
  ]

  extraLegs = {
#  int(len(stepLegs)):"High Mass region",
#  int(len(stepLegs)+1):"MPC",
  int(len(stepLegs)):"MPC+HPC (High Mass)",
  int(len(stepLegs)+1):"MPC+HPC (Low Mass)"
  }


extraStepsN = [(len(stepLegs)+i) for i,x in enumerate(extraSteps)]
stepLegs.update(extraLegs)

