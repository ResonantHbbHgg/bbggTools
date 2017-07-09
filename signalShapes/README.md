###Make limit-like signal shapes for testing and development without running the full limit setting procedure

```
parser.add_argument('-F', '--file', dest='File', required=True, type=str)
parser.add_argument('-L', '--label', dest='Label', required=True, type=str)
parser.add_argument('-C', '--cut', dest='cut', type=str, default=' (1>0) ')
parser.add_argument('--fixPars', dest='fpars', type=str, default='')
parser.add_argument('--ggH', dest='ggh', action='store_true', default=False)

python MakeSignalShape.py -F ../test/RunJobs/EGML_NewMVATraining_Signal_v10/Hadd/output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root \
-C "isSignal && MX < 350 && HHTagger_LM > 0.98 && leadingJet_bDis > 0.56 && subleadingJet_bDis > 0.85" \
-L "SMHH_LowMass_HPC"

python MakeSignalShape.py -F ../test/RunJobs/EGML_NewMVATraining_Signal_v10/Hadd/output_GluGluToHHTo2B2G_node_SM_13TeV-madgraph.root \
-C "isSignal && MX < 350 && HHTagger_LM < 0.98 && HHTagger_LM > 0.5  && leadingJet_bDis > 0.56 && subleadingJet_bDis > 0.85" \
-L "SMHH_LowMass_MPC"

```
