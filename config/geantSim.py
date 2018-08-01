import argparse
simparser = argparse.ArgumentParser()
simparser.add_argument("--bFieldOff", action='store_true', help="Switch OFF magnetic field (default: B field ON)")

simparser.add_argument('-s','--seed', type=int, help='Seed for the random number generator', required=True)
simparser.add_argument('-N','--numEvents', type=int, help='Number of simulation events to run', required=True)
simparser.add_argument('--outName', type=str, help='Name of the output file', required=True)
simparser.add_argument('--detectorPath', type=str, help='Path to detectors', default = "/cvmfs/fcc.cern.ch/sw/releases/0.9.1/x86_64-slc6-gcc62-opt/linux-scientificcernslc6-x86_64/gcc-6.2.0/fccsw-0.9.1-c5dqdyv4gt5smfxxwoluqj2pjrdqvjuj")
simparser.add_argument('--detectorSufix', type=str, help='name of ecal: FCChh_ECalBarrel_withCryostat<suffix>', default = "")


genTypeGroup = simparser.add_mutually_exclusive_group(required = True) # Type of events to generate
genTypeGroup.add_argument("--singlePart", action='store_true', help="Single particle events")
genTypeGroup.add_argument("--pythia", action='store_true', help="Events generated with Pythia")
genTypeGroup.add_argument("--useVertexSmearTool", action='store_true', help="Use the Gaudi Vertex Smearing Tool")

from math import pi
singlePartGroup = simparser.add_argument_group('Single particles')
import sys
singlePartGroup.add_argument('-e','--energy', type=int, required='--singlePart' in sys.argv, help='Energy of particle in GeV')
singlePartGroup.add_argument('--etaMin', type=float, default=0., help='Minimal pseudorapidity')
singlePartGroup.add_argument('--etaMax', type=float, default=0., help='Maximal pseudorapidity')
singlePartGroup.add_argument('--phiMin', type=float, default=0, help='Minimal azimuthal angle')
singlePartGroup.add_argument('--phiMax', type=float, default=2*pi, help='Maximal azimuthal angle')
singlePartGroup.add_argument('--particle', type=int, required='--singlePart' in sys.argv, help='Particle type (PDG)')

pythiaGroup = simparser.add_argument_group('Pythia','Common for min bias and LHE')
pythiaGroup.add_argument('-c', '--card', type=str, default='Generation/data/Pythia_minbias_pp_100TeV.cmd', help='Path to Pythia card (default: PythiaCards/default.cmd)')

simargs, _ = simparser.parse_known_args()

print "=================================="
print "==      GENERAL SETTINGS       ==="
print "=================================="
magnetic_field = not simargs.bFieldOff
num_events = simargs.numEvents
seed = simargs.seed
output_name = simargs.outName
path_to_detector = simargs.detectorPath
detector_sufix = simargs.detectorSufix
print "B field: ", magnetic_field
print "number of events = ", num_events
print "seed: ", seed
print "output name: ", output_name
print "detectors are taken from: ", path_to_detector
print "name of the ECal barrel geometry: ", detector_sufix
if simargs.singlePart:
    energy = simargs.energy
    etaMin = simargs.etaMin
    etaMax = simargs.etaMax
    phiMin = simargs.phiMin
    phiMax = simargs.phiMax
    pdg = simargs.particle
    particle_geant_names = {11: 'e-', -11: 'e+', -13: 'mu+', 13: 'mu-', 22: 'gamma', 111: 'pi0', 211: 'pi+', -211: 'pi-', 130: 'kaon0L'}
    print "=================================="
    print "==       SINGLE PARTICLES      ==="
    print "=================================="
    print "particle PDG, name: ", pdg, " ", particle_geant_names[pdg]
    print "energy: ", energy, "GeV"
    if etaMin == etaMax:
        print "eta: ", etaMin
    else:
        print "eta: from ", etaMin, " to ", etaMax
    if phiMin == phiMax:
        print "phi: ", phiMin
    else:
        print "phi: from ", phiMin, " to ", phiMax
elif simargs.pythia:
    card = simargs.card
    print "=================================="
    print "==            PYTHIA           ==="
    print "=================================="
    print "card = ", card
print "=================================="

from Gaudi.Configuration import *
from GaudiKernel import SystemOfUnits as units
##############################################################################################################
#######                                         GEOMETRY                                         #############
##############################################################################################################
detectors_to_use=[path_to_detector+'/Detector/DetFCChhBaseline1/compact/FCChh_DectEmptyMaster.xml',
                  path_to_detector+'/Detector/DetFCChhTrackerTkLayout/compact/Tracker.xml',
                  path_to_detector+'/Detector/DetFCChhECalInclined/compact/FCChh_ECalBarrel_withCryostat'+str(detector_sufix)+'.xml',
                  path_to_detector+'/Detector/DetFCChhHCalTile/compact/FCChh_HCalBarrel_TileCal.xml',
                  path_to_detector+'/Detector/DetFCChhHCalTile/compact/FCChh_HCalExtendedBarrel_TileCal.xml',
                  path_to_detector+'/Detector/DetFCChhCalDiscs/compact/Endcaps_coneCryo.xml',
                  path_to_detector+'/Detector/DetFCChhCalDiscs/compact/Forward_coneCryo.xml',
                  path_to_detector+'/Detector/DetFCChhTailCatcher/compact/FCChh_TailCatcher.xml',
                  path_to_detector+'/Detector/DetFCChhBaseline1/compact/FCChh_Solenoids.xml',
                  path_to_detector+'/Detector/DetFCChhBaseline1/compact/FCChh_Shielding.xml']

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc", detectors = detectors_to_use, OutputLevel = WARNING)

# ECAL readouts
ecalBarrelReadoutName = "ECalBarrelEta"
ecalBarrelReadoutNamePhiEta = "ECalBarrelPhiEta"
ecalEndcapReadoutName = "EMECPhiEta"
ecalFwdReadoutName = "EMFwdPhiEta"
# HCAL readouts
hcalBarrelReadoutName = "HCalBarrelReadout"
hcalExtBarrelReadoutName = "HCalExtBarrelReadout"
hcalEndcapReadoutName = "HECPhiEta"
hcalFwdReadoutName = "HFwdPhiEta"
# Tail Catcher readout
tailCatcherReadoutName = "Muons_Readout"
# layers to be merged in endcaps, field name of the readout
ecalEndcapNumberOfLayersToMerge = [2] + [2] + [4]*37 + [5]
hcalEndcapNumberOfLayersToMerge = [2] + [4]*19 + [5]
identifierName = "layer"
volumeName = "layer"
##############################################################################################################
#######                                        SIMULATION                                        #############
##############################################################################################################
# Setting random seed, will be propagated to Geant
from Configurables import  RndmGenSvc
from GaudiSvc.GaudiSvcConf import HepRndm__Engine_CLHEP__RanluxEngine_
randomEngine = eval('HepRndm__Engine_CLHEP__RanluxEngine_')
randomEngine = randomEngine('RndmGenSvc.Engine')
randomEngine.Seeds = [seed]

# Magnetic field
from Configurables import SimG4ConstantMagneticFieldTool
if magnetic_field:
    field = SimG4ConstantMagneticFieldTool("bField", FieldOn=True, FieldZMax=20*units.m, IntegratorStepper="ClassicalRK4")
else:
    field = SimG4ConstantMagneticFieldTool("bField", FieldOn=False)

from Configurables import SimG4Svc
geantservice = SimG4Svc("SimG4Svc", detector='SimG4DD4hepDetector', physicslist="SimG4FtfpBert", actions="SimG4FullSimActions", magneticField=field)
# range cut
geantservice.g4PostInitCommands += ["/run/setCut 0.1 mm"]

from Configurables import SimG4Alg, SimG4SaveCalHits, SimG4SingleParticleGeneratorTool, SimG4SaveTrackerHits

savetrackertool = SimG4SaveTrackerHits("saveTrackerHits", readoutNames = ["TrackerBarrelReadout", "TrackerEndcapReadout"]) 
savetrackertool.positionedTrackHits.Path = "TrackerPositionedHits"
savetrackertool.trackHits.Path = "TrackerHits"
savetrackertool.digiTrackHits.Path = "TrackerDigiPostPoint"
saveecaltool = SimG4SaveCalHits("saveECalBarrelHits",readoutNames = [ecalBarrelReadoutName])
saveecaltool.positionedCaloHits.Path = "ECalBarrelPositionedHits"
saveecaltool.caloHits.Path = "ECalBarrelHits"
savehcaltool = SimG4SaveCalHits("saveHCalBarrelHits",readoutNames = [hcalBarrelReadoutName])
savehcaltool.positionedCaloHits.Path = "HCalBarrelPositionedHits"
savehcaltool.caloHits.Path = "HCalBarrelHits"
outputHitsTools = ["SimG4SaveCalHits/saveECalBarrelHits", "SimG4SaveCalHits/saveHCalBarrelHits"]
saveexthcaltool = SimG4SaveCalHits("saveHCalExtBarrelHits",readoutNames = [hcalExtBarrelReadoutName])
saveexthcaltool.positionedCaloHits.Path = "HCalExtBarrelPositionedHits"
saveexthcaltool.caloHits.Path = "HCalExtBarrelHits"
savecalendcaptool = SimG4SaveCalHits("saveECalEndcapHits", readoutNames = [ecalEndcapReadoutName])
savecalendcaptool.positionedCaloHits.Path = "ECalEndcapPositionedHits"
savecalendcaptool.caloHits.Path = "ECalEndcapHits"
savecalfwdtool = SimG4SaveCalHits("saveECalFwdHits", readoutNames = [ecalFwdReadoutName])
savecalfwdtool.positionedCaloHits.Path = "ECalFwdPositionedHits"
savecalfwdtool.caloHits.Path = "ECalFwdHits"
savehcalendcaptool = SimG4SaveCalHits("saveHCalEndcapHits", readoutNames = [hcalEndcapReadoutName])
savehcalendcaptool.positionedCaloHits.Path = "HCalEndcapPositionedHits"
savehcalendcaptool.caloHits.Path = "HCalEndcapHits"
savehcalfwdtool = SimG4SaveCalHits("saveHCalFwdHits", readoutNames = [hcalFwdReadoutName])
savehcalfwdtool.positionedCaloHits.Path = "HCalFwdPositionedHits"
savehcalfwdtool.caloHits.Path = "HCalFwdHits"
savetailcatchertool = SimG4SaveCalHits("saveTailCatcherHits", readoutNames = [tailCatcherReadoutName])
savetailcatchertool.positionedCaloHits = "TailCatcherPositionedHits"
savetailcatchertool.caloHits = "TailCatcherHits"
outputHitsTools += ["SimG4SaveCalHits/saveECalEndcapHits","SimG4SaveCalHits/saveECalFwdHits",
                    "SimG4SaveCalHits/saveHCalExtBarrelHits", "SimG4SaveCalHits/saveHCalEndcapHits",
                    "SimG4SaveCalHits/saveHCalFwdHits", "SimG4SaveCalHits/saveTailCatcherHits", "SimG4SaveTrackerHits/saveTrackerHits"]

geantsim = SimG4Alg("SimG4Alg", outputs = outputHitsTools)

if simargs.singlePart:
    from Configurables import SimG4SingleParticleGeneratorTool
    pgun=SimG4SingleParticleGeneratorTool("SimG4SingleParticleGeneratorTool", saveEdm=True,
                                          particleName=particle_geant_names[pdg], energyMin=energy * 1000, energyMax=energy * 1000,
                                          etaMin=etaMin, etaMax=etaMax, phiMin = phiMin, phiMax = phiMax)
    geantsim.eventProvider = pgun
else:
    from Configurables import PythiaInterface, GenAlg, GaussSmearVertex
    smeartool = GaussSmearVertex("GaussSmearVertex")
    if simargs.useVertexSmearTool:
      smeartool.xVertexSigma = 0.5*units.mm
      smeartool.yVertexSigma = 0.5*units.mm
      smeartool.zVertexSigma = 40*units.mm
      smeartool.tVertexSigma = 180*units.picosecond

    pythia8gentool = PythiaInterface("Pythia8",Filename=card)
    pythia8gen = GenAlg("Pythia8", SignalProvider=pythia8gentool, VertexSmearingTool=smeartool)
    pythia8gen.hepmc.Path = "hepmc"
    from Configurables import HepMCToEDMConverter
    hepmc_converter = HepMCToEDMConverter("Converter", hepmcStatusList=[]) # save all the particles from Pythia
    hepmc_converter.hepmc.Path="hepmc"
    hepmc_converter.genparticles.Path="allGenParticles"
    hepmc_converter.genvertices.Path="GenVertices"
    from Configurables import GenParticleFilter
### Filters generated particles
# accept is a list of particle statuses that should be accepted
    genfilter = GenParticleFilter("StableParticles", accept=[1], OutputLevel=DEBUG)
    genfilter.allGenParticles.Path = "allGenParticles"
    genfilter.filteredGenParticles.Path = "GenParticles"
    from Configurables import SimG4PrimariesFromEdmTool
    particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
    particle_converter.genParticles.Path = "GenParticles"
    geantsim.eventProvider = particle_converter


##############################################################################################################
#######                                       DIGITISATION                                       #############
##############################################################################################################

# Calibration constants
from Configurables import CalibrateInLayersTool, CalibrateCaloHitsTool
calibEcalBarrel = CalibrateInLayersTool("CalibrateEcalBarrel",
                                        # sampling fraction obtained using SamplingFractionInLayers from DetStudies package
                                        samplingFraction = [0.299041341789] + [0.1306220735]+ [0.163243999965]  + [0.186360269398] + [0.203778124831] + [0.216211280314] + [0.227140796653] + [0.243315422934],
                                        readoutName = ecalBarrelReadoutName,
                                        layerFieldName = "layer")
calibHcells = CalibrateCaloHitsTool("CalibrateHCal", invSamplingFraction="41.66")
calibEcalEndcap = CalibrateCaloHitsTool("CalibrateECalEndcap", invSamplingFraction="13.89")
calibEcalFwd = CalibrateCaloHitsTool("CalibrateECalFwd", invSamplingFraction="303.03")
calibHcalEndcap = CalibrateCaloHitsTool("CalibrateHCalEndcap", invSamplingFraction="33.62")
calibHcalFwd = CalibrateCaloHitsTool("CalibrateHCalFwd", invSamplingFraction="1207.7")

# Create cells
from Configurables import CreateCaloCells
# -> ECal barrel
# 1. step - merge hits into cells with default Eta segmentation
createEcalBarrelCellsStep1 = CreateCaloCells("EcalBarrelCellsStep1",
                                             doCellCalibration=True,
                                             calibTool=calibEcalBarrel,
                                             addCellNoise=False, filterCellNoise=False)
createEcalBarrelCellsStep1.hits.Path="ECalBarrelHits"
createEcalBarrelCellsStep1.cells.Path="ECalBarrelCellsStep1"
# 2. step - rewrite the cellId using the Phi-Eta segmentation (remove 'module')
# 2.1. retrieve phi positions from centres of cells
from Configurables import CreateVolumeCaloPositions
positionsEcalBarrel = CreateVolumeCaloPositions("positionsEcalBarrel")
positionsEcalBarrel.hits.Path = "ECalBarrelCellsStep1"
positionsEcalBarrel.positionedHits.Path = "ECalBarrelPositions"
# 2.2. assign cells into phi bins
from Configurables import RedoSegmentation
resegmentEcalBarrel = RedoSegmentation("ReSegmentationEcalBarrel",
                                       oldReadoutName = 'ECalBarrelEta',
                                       oldSegmentationIds = ['module'],
                                       newReadoutName = 'ECalBarrelPhiEta',
                                       inhits = "ECalBarrelPositions",
                                       outhits = "ECalBarrelCellsStep2")
# 3. step - merge cells in the same phi bin
createEcalBarrelCells = CreateCaloCells("CreateECalBarrelCells",
                                        doCellCalibration=False, # already calibrated in step 1
                                        addCellNoise=False, filterCellNoise=False,
                                        hits="ECalBarrelCellsStep2",
                                        cells="ECalBarrelCells")

# Create Ecal cells in endcaps                              
# 1. step - merge layer IDs                                 
# 2. step - create cells                                    
from Configurables import MergeLayers
mergeLayersEcalEndcap = MergeLayers("MergeLayersEcalEndcap",
                   # take the bitfield description from the geometry service
                   readout = ecalEndcapReadoutName,
                   # cells in which field should be merged                  
                   identifier = identifierName,
                   volumeName = volumeName,
                   # how many cells to merge                                
                   merge = ecalEndcapNumberOfLayersToMerge,
                   OutputLevel = INFO)
mergeLayersEcalEndcap.inhits.Path = "ECalEndcapHits"
mergeLayersEcalEndcap.outhits.Path = "mergedECalEndcapHits"
createEcalEndcapCells = CreateCaloCells("CreateEcalEndcapCells",
                                        doCellCalibration=True,
                                        calibTool=calibEcalEndcap,
                                        addCellNoise=False, filterCellNoise=False)
createEcalEndcapCells.hits.Path="mergedECalEndcapHits"
createEcalEndcapCells.cells.Path="ECalEndcapCells"
# -> Ecal forward
createEcalFwdCells = CreateCaloCells("CreateEcalFwdCells",
                                     doCellCalibration=True,
                                     calibTool=calibEcalFwd,
                                     addCellNoise=False, filterCellNoise=False)
createEcalFwdCells.hits.Path="ECalFwdHits"
createEcalFwdCells.cells.Path="ECalFwdCells"

# -> Hcal barrel
createHcalCells = CreateCaloCells("CreateHCalBarrelCells",
                                  doCellCalibration=True,
                                  calibTool=calibHcells,
                                  addCellNoise = False, filterCellNoise = False,
                                  hits="HCalBarrelHits",
                                  cells="HCalBarrelCells")

# -> Hcal extended barrel
createExtHcalCells = CreateCaloCells("CreateExtHcalCaloCells",
                                     doCellCalibration=True,
                                     calibTool=calibHcells,
                                     addCellNoise = False, filterCellNoise = False,
                                     hits="HCalExtBarrelHits",
                                     cells="HCalExtBarrelCells")

# Create Hcal cells in endcaps
# 1. step - merge layer IDs
# 2. step - create cells
mergeLayersHcalEndcap = MergeLayers("MergeLayersHcalEndcap",
                   # take the bitfield description from the geometry service
                   readout = hcalEndcapReadoutName,
                   # cells in which field should be merged                  
                   identifier = identifierName,
                   volumeName = volumeName,
                   # how many cells to merge                                
                   merge = hcalEndcapNumberOfLayersToMerge,
                   OutputLevel = INFO)
mergeLayersHcalEndcap.inhits.Path = "HCalEndcapHits"
mergeLayersHcalEndcap.outhits.Path = "mergedHCalEndcapHits"
createHcalEndcapCells = CreateCaloCells("CreateHcalEndcapCells",
                                        doCellCalibration=True,
                                        calibTool=calibHcalEndcap,
                                        addCellNoise=False, filterCellNoise=False)
createHcalEndcapCells.hits.Path="mergedHCalEndcapHits"
createHcalEndcapCells.cells.Path="HCalEndcapCells"
# -> Hcal forward 
createHcalFwdCells = CreateCaloCells("CreateHcalFwdCaloCells",
                                     doCellCalibration=True,
                                     calibTool=calibHcalFwd,
                                     addCellNoise=False, filterCellNoise=False)
createHcalFwdCells.hits.Path="HCalFwdHits"
createHcalFwdCells.cells.Path="HCalFwdCells"
# -> Tail Catcher
createTailCatcherCells = CreateCaloCells("CreateTailCatcherCells",
                                         doCellCalibration=False,
                                         addCellNoise = False, filterCellNoise = False,
                                         OutputLevel = INFO,
                                         hits="TailCatcherHits",
                                         cells="TailCatcherCells")

# PODIO algorithm
from Configurables import ApplicationMgr, FCCDataSvc, PodioOutput
podioevent = FCCDataSvc("EventDataSvc")
out = PodioOutput("out")
out.outputCommands = ["drop *",
                      "keep GenParticles",
                      "keep GenVertices",
                      "keep TrackerHits",
                      "keep TrackerPositionedHits",
                      "keep TrackerDigiPostPoint",
                      "keep ECalBarrelCells",
                      "keep ECalEndcapCells",
                      "keep ECalFwdCells",
                      "keep HCalBarrelCells",
                      "keep HCalExtBarrelCells",
                      "keep HCalEndcapCells",
                      "keep HCalFwdCells",
                      "keep TailCatcherCells"]
out.filename = output_name

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
geantsim.AuditExecute = True
createEcalBarrelCellsStep1.AuditExecute = True
positionsEcalBarrel.AuditExecute = True
resegmentEcalBarrel.AuditExecute = True
createEcalBarrelCells.AuditExecute = True
createHcalCells.AuditExecute = True
mergeLayersEcalEndcap.AuditExecute = True
createEcalEndcapCells.AuditExecute = True
mergeLayersHcalEndcap.AuditExecute = True
createHcalEndcapCells.AuditExecute = True
createHcalFwdCells.AuditExecute = True
createTailCatcherCells.AuditExecute = True
out.AuditExecute = True

list_of_algorithms = [geantsim,
                      createEcalBarrelCellsStep1,
                      positionsEcalBarrel,
                      resegmentEcalBarrel,
                      createEcalBarrelCells,
                      createHcalCells,
                      mergeLayersEcalEndcap,
                      createEcalEndcapCells,
                      createEcalFwdCells,
                      createExtHcalCells,
                      mergeLayersHcalEndcap,
                      createHcalEndcapCells,
                      createHcalFwdCells,
                      createTailCatcherCells,
                      out]

if simargs.pythia:
    list_of_algorithms = [pythia8gen, hepmc_converter, genfilter] + list_of_algorithms

ApplicationMgr(
    TopAlg = list_of_algorithms,
    EvtSel = 'NONE',
    EvtMax   = num_events,
    ExtSvc = [podioevent, geoservice, geantservice],
)
