##################################################
##################################################
## Configuration parameters to run MakeStack.py ##
##################################################
##################################################

doBlind = False
doSignalRegion = True
doJetCR = False

isPhoCR = False
addHiggs = True
hideData = False
addbbH = True
dyjets = False

#do pile up reweighting
doPUweight = True

year = ""

doSignal = True

#btagging working poing
# 0.46 - loose
# 0.80 - medium
# 0.935 - tight
BTAG = 0.

#Luminosity to normalize backgrounds
lumi = 35900#pb
MCSF = 1.0
signalFactor = 1
#List of datasets to be used (cross section information defined there)
data_file = open("datasets/datasets80X_ApprovalNoBkg_Restricted.json")

#number of bins in histograms
nbin = 30
dr = "sqrt( (leadingPhoton.Eta() - subleadingPhoton.Eta())*(leadingPhoton.Eta() - subleadingPhoton.Eta()) + (leadingPhoton.Phi() - subleadingPhoton.Phi())*(leadingPhoton.Phi() - subleadingPhoton.Phi()) )"

#plots will be saved in dirName
prefix = ""
dirSuffix = "June7_36ifb_UNBLIND_Mjj70_MVA_NOBKG"
dirPrefix = "/afs/cern.ch/work/m/mgouzevi/private/LIMITS/CLEAN_EPS_GGBB/SAMPLES_PRODUCER/CMSSW_8_0_26_patch1/src/flashgg/bbggTools/stackPlots/"
dirName = dirPrefix + dirSuffix

#Location of root files for each invidivual samples. Name of the root files is defined in datasets/datasets(76).json
loc = '/eos/cms/store/group/phys_higgs/resonant_HH/RunII/FlatTrees/2016/May2_Mjj70to190_NewCatMVA'
higgsLocation = loc + '/EGML_Background_Mjj70_NewMVA/Hadd/'
bkgLocation = loc + '/EGML_Background_Mjj70_NewMVA/Hadd/'
signalLocation = loc + '/EGML_Signal_Mjj70_NewMVA/Hadd/'
dataLocation = loc + '/EGML_Data_Mjj70_NewMVA/Hadd/'


#plots to be made
plots = []
##plots.append(["dicandidate_Mass", "diHiggsCandidate.M()", "m_{#gamma#gamma jj} [GeV]", 34, 150, 1000])
plots.append(["MXprime", "diHiggsCandidate.M() - dijetCandidate.M() - diphotonCandidate.M() + 250.", "#tilde{M}_{X} [GeV]", 30, 250, 1150])
plots.append(["diPho_Mass", "diphotonCandidate.M()", "m_{#gamma#gamma} [GeV]", 80, 100, 140])
plots.append(["diJet_Mass", "dijetCandidate.M()", "m_{jj} [GeV]", 24, 70, 190])
plots.append(["HHTagger", "HHTagger", "Classification MVA", 54, -1.08, 1.08])
plots.append(["CosTheta_gg", "fabs(CosTheta_gg)", "| cos #theta_{#gamma#gamma} |", 10, 0, 1])
plots.append(["CosTheta_bb", "fabs(CosTheta_bb)", "| cos #theta_{jj} |", 10, 0, 1])
plots.append(["CosThetaStar_CS", "fabs(CosThetaStar_CS)", "| cos #theta^{CS}_{HH} |", 10, 0, 1])


#plots.append(["HHTagger", "HHTagger", "Categorization MVA", 100, -1, 1])
#plots.append(["HHTagger_LM", "HHTagger_LM", "Categorization MVA (Low Mass Training)", 50, -1, 1])
#plots.append(["HHTagger_HM", "HHTagger_HM", "Categorization MVA (High Mass Training)", 50, -1, 1])
#plots.append(["MXprime_binned", "diHiggsCandidate.M() - dijetCandidate.M() - diphotonCandidate.M() + 250.", "#tilde{M}_{X} (GeV)", 80, 200, 300])


#plots.append(["PhotonIDMVA2", "(subleadingPhotonIDMVA)", "2 Photon #gammaMVA discriminant", nbin, 0.2, 1])
#plots.append(["PhotonIDMVA", "(leadingPhotonIDMVA+subleadingPhotonIDMVA)", "Sum Photon #gammaMVA discriminant", nbin, 0, 2])
#plots.append(["PhotonIDMVA1", "(leadingPhotonIDMVA)", "1 Photon #gammaMVA discriminant", nbin, 0.2, 1])
#plots.append(["leadingPhoton_pt", "leadingPhoton.pt()", "p_{T}(#gamma_{1}) [GeV]", 50, 30, 150])

#plots.append(["MX", "diHiggsCandidate.M() - dijetCandidate.M() + 125.", "#tilde{M}_{X} (GeV)", 40, 200, 1000])
#plots.append(["binnedMX", "diHiggsCandidate.M() - dijetCandidate.M() + 125.", "#tilde{M}_{X} (GeV)", nbin, 100, 1000])

#plots.append(["CosThetaStar", "CosThetaStar", "Cos(#theta*)", nbin, -1, 1])

#plots.append(["CosTheta_bbgg", "CosTheta_bbgg", "Cos(#theta_{jj#gamma#gamma})", nbin, -1, 1])
#plots.append(["CosTheta_ggbb", "CosTheta_ggbb", "Cos(#theta_{#gamma#gammajj})", nbin, -1, 1])
#plots.append(["Phi", "Phi0", "#Phi", nbin, -3.5, 3.5])
#plots.append(["Phi1", "Phi1", "#Phi_{1}", nbin, -3.5, 3.5])
#plots.append(["DiJetDiPho_DR", "DiJetDiPho_DR", "#DeltaR(#gamma#gamma,jj)", nbin, 0, 4])
#plots.append(["nvtx", "nvtx", "Number of vertices", 50, 0, 50])
#plots.append(["costhetastar_cs", "fabs(CosThetaStar_CS)", "|cos#theta*|_{CS}", nbin, 0, 1])
#plots.append(["costhetastar", "fabs(CosThetaStar)", "|cos#theta*|", nbin, 0, 1])
#plots.append(["j1ratio_dijet", "leadingJet.Pt()/dijetCandidate.M()", "p_{T}(j_{1})/M(jj)", nbin, 0.1, 1.5])
#plots.append(["dijet_deta", "fabs(leadingJet.Eta() - subleadingJet.Eta())", "#Delta#eta between jets", nbin, 0, 5])
#plots.append(["diPho_Mass", "diphotonCandidate.M()", "M(#gamma#gamma) [GeV]", nbin, 100, 180])
#plots.append(["costhetastar", "fabs(CosThetaStar)", "|cos#theta*|", nbin, 0, 1])
#plots.append(["diPho_Mass", "diphotonCandidate.M()", "DiPhoton Candidate Mass (GeV)", nbin, 100, 180])
#plots.append(["diPho_Mass_HM", "diphotonCandidate.M()", "DiPhoton Candidate Mass (GeV)", nbin, 80, 2000])
#plots.append(["diJet_Mass_Limit", "dijetCandidate.M()", "DiJet Candidate Mass (GeV)", nbin, 60, 180])
#plots.append(["diJet_Mass_HM", "dijetCandidate.M()", "DiJet Candidate Mass (GeV)", nbin, 80, 2000])
#plots.append(["leadingJet_pt", "leadingJet.pt()", "p_{T}(j_{1}) [GeV]", nbin, 15, 200] )
#plots.append(["subleadingJet_pt", "subleadingJet.pt()", "p_{T}(j_{2}) (GeV)", nbin, 15, 80] )
#plots.append(["leadingJet_eta", "leadingJet.Eta()", "#eta(j_{1})", nbin, -3, 3] )
#plots.append(["leadingPhoton_eta", "leadingPhoton.Eta()", "#eta(#gamma_{1})", nbin, -3, 3] )
#plots.append(["MX", "diHiggsCandidate.M() - dijetCandidate.M() + 125.", "#tilde{M}_{X} (GeV)", nbin, 100, 1000])
#plots.append(["MKF", "diHiggsCandidate_KF.M()", "M_{KinFit}(bb#gamma#gamma) (GeV)", nbin, 100, 1000])
#plots.append(["dicandidate_Mass_Limit", "diHiggsCandidate.M()", "M(hh) [GeV]", nbin, 225, 350])
#plots.append(["dicandidate_Mass_HM", "diHiggsCandidate.M()", "M(hh) [GeV]", nbin, 250, 5000])
#plots.append(["btagSum", "leadingJet_bDis+subleadingJet_bDis", "Sum of b-tag of jet pair", nbin, 0, 2])
#plots.append(["subleadingPhoton_pt", "subleadingPhoton.pt()", "p_{T}(#gamma_{2}) [GeV]", nbin, 30, 150])
#plots.append(["subleadingPhoton_eta", "subleadingPhoton.eta()", "#eta(#gamma_{2})", nbin, -3, 3])
#plots.append(["dr_photons", dr, "#DeltaR between photons", nbin, 0, 10])
#plots.append(["leadingPho_MVA", "customLeadingPhotonIDMVA", "Leading Photon #gammaMVA discriminant", nbin, 0, 1])
#plots.append(["subleadingPho_MVA", "customSubLeadingPhotonIDMVA", "SubLeading Photon #gammaMVA discriminant", nbin, 0, 1])


#cuts to be used to make plots
Cut = " isSignal && diphotonCandidate.M() > 100 && diphotonCandidate.M() < 180 "
Cut += " && dijetCandidate.M() > 70 && dijetCandidate.M() < 190 "
#Cut += " & (diHiggsCandidate.M()-dijetCandidate.M()-diphotonCandidate.M()+250) < 350 "
#Cut += " & HHTagger_HM < 0.97 && HHTagger_HM > 0.6 "
#Cut += " & HHTagger_LM < 0.985 && HHTagger_LM > 0.6 && leadingJet_bDis > 0.55 && subleadingJet_bDis > 0.55 "
#Cut += " && (((leadingJet_bDis > 0.8 && subleadingJet_bDis > 0.8) && (leadingJet_bDis < 0.92))+((leadingJet_bDis > 0.8 && subleadingJet_bDis > 0.8) && (subleadingJet_bDis < 0.92)))"
#Cut += " && (diHiggsCandidate.M() - dijetCandidate.M() - diphotonCandidate.M() + 250) < 350"
#Cut += " && (diHiggsCandidate.M() - dijetCandidate.M() + 125.) > 280 && (diHiggsCandidate.M() - dijetCandidate.M() + 125.) < 320"
#Cut += " && diHiggsCandidate.M() > 280 && diHiggsCandidate.M() < 320"
#Cut += " && (leadingJet.pt()/dijetCandidate.M()) > 0.3333"
#Cut += " && (subleadingJet.pt()/dijetCandidate.M()) > 0.25"
#Cut += " && leadingPhoton.pt() > 35 && subleadingPhoton.pt() > 35 "
#Cut += " && leadingJet.pt() > 35 && subleadingJet.pt() > 35 "
#Cut += " && leadingJet.Eta() < 2. && leadingJet.Eta() > -2"
#Cut += " && leadingPhotonID[1] == 1 "
#Cut += " && subleadingPhotonID[1] == 1 "
#Cut += " && leadingPhotonISO[1] == 1 "
#Cut += " && subleadingPhotonISO[1] == 1 "
#Cut += " && leadingPhotonEVeto == 0 "
#Cut += " && subleadingPhotonEVeto == 0 "
#Cut += " && leadingPhotonEVeto == 1 "
#Cut += " && subleadingPhotonEVeto == 1 "
#Cut += " && leadingPhotonISO[0] == 1 "
#Cut += " && subleadingPhotonISO[0] == 1 "
#Cut += " && leadingJet.pt() > 45 && subleadingJet.pt() > 45 "
#Cut += " && leadingJet_bDis > 0.9 && subleadingJet_bDis > 0.9"
