import ROOT, sys
from bin.tools import train

#### input variables to be used for training ###
variables = [
    'emax_l00',
    'e2max_l00',
    'eocore_l00',
    'edmax_l00',
    'w3st_l00',
    'w21st_l00',

    'emax_l01',
    'e2max_l01',
    'eocore_l01',
    'edmax_l01',
    'w9st_l01',
    'w41st_l01',

    'emax_l02',
    'e2max_l02',
    'eocore_l02',
    'edmax_l02',
    'w3st_l02',
    'w21st_l02',

    'emax_l03',
    'e2max_l03',
    'eocore_l03',
    'edmax_l03',
    'w3st_l03',
    'w21st_l03',

    'emax_l04',
    'e2max_l04',
    'eocore_l04',
    'edmax_l04',
    'w3st_l04',
    'w21st_l04',

    'e_l01',
    'e_l21',
    'e_l31',
    'e_l41',
    'e_l51',
    'e_l61',
    'e_l71',

    'e_l0T',
    'e_l1T',
    'e_l2T',
    'e_l3T',
    'e_l4T',
    'e_l5T',
    'e_l6T',
    'e_l7T',
]

#### define signal and background trees + selection to be applied ###

eospath = '/eos/experiment/fcc/hh/simulation/samples/v03_ecal_half2nd3rdLayerSmallPhiAllLayers/singlePart/'
#eospath = '/afs/cern.ch/user/d/djamin/fcc_work/heppy/FCChhAnalyses/output/tagger/W_top_vs_QCD_tagger/'
# to make these special trees ,execute scripts/do_jet_list.py
bkgTree = '/afs/cern.ch/user/r/rastein/public/pionrejection/v03_ecal_half2nd3rdLayerSmallPhiAllLayers/pi0/10GeV/pi0_half2nd3rdLayerSmallPhiAllLayers_10GeV.root'
sigTree = '/eos/experiment/fcc/hh/simulation/samples/v03_ecal_half2nd3rdLayerSmallPhiAllLayers/singlePart/photon/bFieldOn/eta0/10GeV/ntup/pi0rejection/photon_half2nd3rdLayerSmallPhiAllLayers_10GeV.root'

# take only hadronic decays
SIGcuts = ''
#SIGcuts = 'fullhad_fullhadsemilep_lep_decays<2 && Jet_trk02_tau21>0 && Jet_trk02_tau31>0 && Jet_trk02_tau32>0'
BKGcuts = SIGcuts
cuts    = [SIGcuts,BKGcuts]

#### define MVA method ###

method = 'BDT'
#nTrain = 50000
nTrain = 0
label = 'pi0_vs_photon_half2nd3rdLayerSmallPhiAllLayers_10GeV'


#### train the MVA ####
train(bkgTree, sigTree, variables, method, nTrain, label, cuts)


#### check the output with the TMVAGUI ####

#outFileName = 'TMVA_{}_{}.root'.format(method, label)
#ROOT.TMVA.TMVAGui(outFileName, label)
#raw_input('Press Enter to continue...')
