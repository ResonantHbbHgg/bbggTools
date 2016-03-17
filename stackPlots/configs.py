##################################################
##################################################
## Configuration parameters to run MakeStack.py ##
##################################################
##################################################

doBlind = True
doSignalRegion = True
doJetCR = False
isPhoCR = False
addHiggs = True
hideData = False
addbbH = False

lumi = 2700#pb

data_file = open("datasets/datasets76X.json")

bkgsToAdd = [['vbf', "VBF H(#gamma#gamma)", 1], ['vh', "VH(#gamma#gamma)", 1], ['ggf', "ggH(#gamma#gamma)", 1], ['tth', "ttH(#gamma#gamma)", 1], ["DYJet", "Z/#gamma*+Jets", 1], ["GJets", "#gamma+Jets", 2], ["DiPhoJet", "#gamma#gamma+Jets", 1], ["QCD", "QCD", 2] ]

nbin = 40
dr = "sqrt( (leadingPhoton.Eta() - subleadingPhoton.Eta())*(leadingPhoton.Eta() - subleadingPhoton.Eta()) + (leadingPhoton.Phi() - subleadingPhoton.Phi())*(leadingPhoton.Phi() - subleadingPhoton.Phi()) )"

prefix = ""
dirSuffix = "1_3_0_SR_PhoMVAID_76X"
dirPrefix = "/afs/cern.ch/user/r/rateixei/www/HHBBGG/"
dirName = dirPrefix + dirSuffix

higgsLocation = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/flg76X/CMSSW_7_6_3/src/flashgg/bbggTools/test/RunJobs/mva_hig/Hadd/"
bkgLocation = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/flg76X/CMSSW_7_6_3/src/flashgg/bbggTools/test/RunJobs/mva_bkg/Hadd/"
signalLocation = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/flg76X/CMSSW_7_6_3/src/flashgg/bbggTools/test/RunJobs/mva_sig/Hadd/"
dataLocation = "/afs/cern.ch/work/r/rateixei/work/DiHiggs/flg76X/CMSSW_7_6_3/src/flashgg/bbggTools/test/RunJobs/mva_data/Hadd/"

plots = []
plots.append(["j1ratio_dijet", "leadingJet.Pt()/dijetCandidate.M()", "p_{T}(j_{1})/M(jj)", nbin, 0.1, 1.5])
plots.append(["dijet_deta", "fabs(leadingJet.Eta() - subleadingJet.Eta())", "#Delta#eta between jets", nbin, 0, 5])
plots.append(["diPho_Mass", "diphotonCandidate.M()", "DiPhoton Candidate Mass (GeV)", nbin, 100, 180])
#plots.append(["costhetastar", "fabs(CosThetaStar)", "|cos#theta*|", nbin, 0, 1])
plots.append(["leadingPhoton_pt", "leadingPhoton.pt()", "Leading Photon Pt (GeV)", nbin, 30, 300])
plots.append(["nvtx", "nvtx", "Number of vertices", 30, 0, 30])
plots.append(["diPho_Mass", "diphotonCandidate.M()", "DiPhoton Candidate Mass (GeV)", nbin, 100, 180])
plots.append(["diPho_Mass_HM", "diphotonCandidate.M()", "DiPhoton Candidate Mass (GeV)", nbin, 80, 2000])
plots.append(["diJet_Mass", "dijetCandidate.M()", "DiJet Candidate Mass (GeV)", nbin, 60, 180])
plots.append(["diJet_Mass_Limit", "dijetCandidate.M()", "DiJet Candidate Mass (GeV)", nbin, 60, 180])
plots.append(["diJet_Mass_HM", "dijetCandidate.M()", "DiJet Candidate Mass (GeV)", nbin, 80, 2000])
plots.append(["leadingJet_pt", "leadingJet.pt()", "Leading Jet Pt (GeV)", nbin, 15, 200] )
plots.append(["subleadingJet_pt", "subleadingJet.pt()", "Subleading Jet Pt (GeV)", nbin, 15, 80] )
plots.append(["leadingJet_eta", "leadingJet.Eta()", "Leading Jet Eta", nbin, -3, 3] )
plots.append(["leadingPhoton_eta", "leadingPhoton.Eta()", "Leading Photon Eta", nbin, -3, 3] )
plots.append(["dicandidate_Mass", "diHiggsCandidate.M()", "DiHiggs Candidate Mass (GeV)", nbin, 100, 1000])
plots.append(["MX", "diHiggsCandidate.M() - dijetCandidate.M() + 125.", "M_{X} (GeV)", nbin, 100, 1000])
plots.append(["MKF", "diHiggsCandidate_KF.M()", "M_{KinFit}(bb#gamma#gamma) (GeV)", nbin, 100, 1000])
plots.append(["dicandidate_Mass_Limit", "diHiggsCandidate.M()", "DiHiggs Candidate Mass (GeV)", nbin, 225, 350])
plots.append(["dicandidate_Mass_HM", "diHiggsCandidate.M()", "DiHiggs Candidate Mass (GeV)", nbin, 250, 5000])
plots.append(["btagSum", "leadingJet_bDis+subleadingJet_bDis", "Sum of b-tag of jet pair", nbin, 0, 2])
plots.append(["subleadingPhoton_pt", "subleadingPhoton.pt()", "Subleading Photon Pt (GeV)", nbin, 30, 150])
plots.append(["subleadingPhoton_eta", "subleadingPhoton.eta()", "Subleading Photon Eta)", nbin, -3, 3])
plots.append(["dr_photons", dr, "#DeltaR between photons", nbin, 0, 10])
plots.append(["leadingPho_MVA", "leadingPhotonIDMVA", "Leading Photon #gammaMVA discriminant", nbin, 0, 1])
plots.append(["subleadingPho_MVA", "subleadingPhotonIDMVA", "SubLeading Photon #gammaMVA discriminant", nbin, 0, 1])

Cut = " diphotonCandidate.M() > 100 && diphotonCandidate.M() < 180"
Cut += " && dijetCandidate.M() > 60 && dijetCandidate.M() < 180"
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
