import ROOT, sys
from bin.tools import train

#### input variables to be used for training ###
variables = [
    'emax_l00',
    'e2max_l00',
    'eocore_l00',
    'edmax_l00',
    'w9st_l00',
    'w41st_l00',

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
    'w9st_l02',
    'w41st_l02',

    'emax_l03',
    'e2max_l03',
    'eocore_l03',
    'edmax_l03',
    'w9st_l03',
    'w41st_l03',

    'e_l01',
    'e_l20',
    'e_l30',

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

eospath = '/afs/cern.ch/user/r/rastein/public/pionrejection/v03_ecal_half2ndLayer/'
#eospath = '/afs/cern.ch/user/d/djamin/fcc_work/heppy/FCChhAnalyses/output/tagger/W_top_vs_QCD_tagger/'
# to make these special trees ,execute scripts/do_jet_list.py
bkgTree = eospath + 'pi0/10GeV/pi0_half2ndLayer_10GeV.root'
sigTree = eospath + 'photon/10GeV/photon_half2ndLayer_10GeV.root'

# take only hadronic decays
SIGcuts = ''
#SIGcuts = 'fullhad_fullhadsemilep_lep_decays<2 && Jet_trk02_tau21>0 && Jet_trk02_tau31>0 && Jet_trk02_tau32>0'
BKGcuts = SIGcuts
cuts    = [SIGcuts,BKGcuts]

#### define MVA method ###

method = 'BDT'
#nTrain = 50000
nTrain = 0
label = 'pi0_vs_photon_half2ndLayer_10GeV'


#### train the MVA ####
train(bkgTree, sigTree, variables, method, nTrain, label, cuts)


#### check the output with the TMVAGUI ####

#outFileName = 'TMVA_{}_{}.root'.format(method, label)
#ROOT.TMVA.TMVAGui(outFileName, label)
#raw_input('Press Enter to continue...')
