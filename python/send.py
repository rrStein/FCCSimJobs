path_to_INIT = '/cvmfs/fcc.cern.ch/sw/views/releases/0.9.1/x86_64-slc6-gcc62-opt/setup.sh'
path_to_LHE = '/afs/cern.ch/work/h/helsens/public/FCCsoft/FlatGunLHEventProducer/'
path_to_FCCSW = '/cvmfs/fcc.cern.ch/sw/releases/0.9.1/x86_64-slc6-gcc62-opt/linux-scientificcernslc6-x86_64/gcc-6.2.0/fccsw-0.9.1-c5dqdyv4gt5smfxxwoluqj2pjrdqvjuj'
yamldir='/afs/cern.ch/work/h/helsens/public/FCCDicts/yaml/FCC/simu/'

import glob, os, sys
import commands
import time
import random
import json
import argparse
import imp
import decimal
import re
import yaml

from math import pi
from datetime import datetime

import utils as ut
import makeyaml as my
import users as us

#__________________________________________________________
def warning(msg, ifConfirm = False):
    print "=================================="
    print "==           WARNING           ==="
    print "=================================="
    print msg
    print "=================================="
    if ifConfirm:
        choice = raw_input("Do you want to exit and adjust arguments? [y/n] ")
        if choice.lower() == 'y':
            exit()

#__________________________________________________________
def getInputFiles(path):
    files = []
    All_files = glob.glob("%s/merge.yaml"%path)
    if len(All_files)==0:
        print 'there is no files checked for process %s exit'%path
        sys.exit(3)

    if len(All_files)>1:
        print 'not sure how you can get two merged files %s exit'%path
        sys.exit(3)

    with open(All_files[0], 'r') as stream:
        try:
            tmpf = yaml.load(stream)
            outfiles=tmpf['merge']['outfiles']
            for f in outfiles:
                fname='%s%s'%(tmpf['merge']['outdir'],f[0])
                files.append(fname)
        except yaml.YAMLError as exc:
            print(exc)

    return files

#__________________________________________________________
def takeOnlyNonexistingFiles(files, output):
    All_files = glob.glob("%s/output*.yaml"%output)
    if len(All_files)==0: return files, len(files)
    newlist=[]
    for f in files:
        found=False
        tocheck=f.split('/')[-1]
        tocheck=tocheck.replace('.root','')
        for ff in All_files:
            if tocheck in ff: found=True
        if not found:newlist.append(f)
    return newlist, len(newlist)


#__________________________________________________________
def getJobInfo(argv):
    if '--recPositions' in argv:
        default_options = 'config/recPositions.py'
        job_type = "ntup/positions"
        short_job_type = "recPos"
        return default_options,job_type,short_job_type,False

    elif '--pionRejection' in argv:
        default_options = ''
        job_type = "ntup/pi0rejection"
        short_job_type = "pi0rej"
        return default_options,job_type,short_job_type,False

    elif '--recSlidingWindow' in argv:
        default_options = 'config/recSlidingWindow.py'
        if '--noise' in argv:
            job_type = "reco/slidingWindow/electronicsNoise"
        elif '--addPileupNoise' in argv:
            job_type = "reco/slidingWindow/pileupNoise"
        else:
            job_type = "reco/slidingWindow/noNoise"
        short_job_type = "recWin"
        return default_options,job_type,short_job_type,False

    elif '--estimatePileup' in sys.argv:
        default_options = 'config/preparePileup.py'
        job_type = "ana/pileup_cluster"
        short_job_type = "pileup"
        return default_options,job_type,short_job_type,False

    elif '--mergePileup' in sys.argv:
        default_options = 'config/overlayPileup.py'
        job_type = "simuPU"
        short_job_type = "pileup"
        return default_options,job_type,short_job_type,False

    elif '--recTopoClusters' in argv:
        default_options = 'config/recTopoClusters.py'
        job_type = "ntup/topoClusters"
        if '--noise' in argv:
            job_type = "reco/topoClusters/electronicsNoise"
        elif '--addPileupNoise' in argv:
            job_type = "reco/topoClusters/pileupNoise/"
        else:
            job_type = "reco/topoClusters/noNoise"
        short_job_type = 'recTopo'
        return default_options,job_type,short_job_type,False

    elif '--ntuple' in argv:
        default_options = '....' # TODO how ?
        job_type = "ntup"
        short_job_type = "ntup"
        return default_options,job_type,short_job_type,False

    elif '--trackerPerformance' in argv:
        default_options = 'config/geantSim_trackerPerformance.py'
        if "--tripletTracker" in argv:
          job_type="simu/trkPerf_triplet"
          short_job_type = "sim"
        else:
          job_type="simu/trkPerf_v3_03"
          short_job_type = "sim"
        return default_options,job_type,short_job_type,True

    else:
        default_options = 'config/geantSim.py'
        job_type = "simu"
        short_job_type = "sim"
        return default_options,job_type,short_job_type,True


#__________________________________________________________
if __name__=="__main__":
    Dir = os.getcwd()

    user=os.environ['USER']
    userext=-999999
    for key, value in us.users.iteritems():
        if key==user:
            userext=value
    if userext<0:
        print 'user not known ',user,'   exit (needs to be added to users.py)'
        sys.exit(3)

    parser = argparse.ArgumentParser()
    parser.add_argument('--local', type=str, help='Use local FCCSW installation, need to provide a file with path_to_INIT and or path_to_FCCSW')
    parser.add_argument('--version', type=str, default = "v03", help='Specify the version of FCCSimJobs')

    parser.add_argument("--bFieldOff", action='store_true', help="Switch OFF magnetic field (default: B field ON)")
    parser.add_argument("--pythiaSmearVertex", action='store_true', help="Write vertex smearing parameters to pythia config file")
    parser.add_argument("--no_eoscopy", action='store_true',  help="DON'T copy result files to eos")

    parser.add_argument('-n','--numEvents', type=int, help='Number of simulation events per job', required='--recPositions' not in sys.argv and '--recSlidingWindow' not in sys.argv and '--recTopoClusters' not in sys.argv and '--ntuple' not in sys.argv)
    parser.add_argument('-N','--numJobs', type=int, default = 10, help='Number of jobs to submit')
    parser.add_argument('-o','--output', type=str, help='Path of the output on EOS', default="/eos/experiment/fcc/hh/simulation/samples/")
    parser.add_argument('-l','--log', type=str, help='Path of the logs', default = "BatchOutputs/")


    jobTypeGroup = parser.add_mutually_exclusive_group() # Type of job: simulation or reconstruction
    jobTypeGroup.add_argument("--sim", action='store_true', help="Simulation (default)")
    jobTypeGroup.add_argument("--pionRejection", action='store_true', help="Convert to pion0 rejection variables")
    jobTypeGroup.add_argument("--recPositions", action='store_true', help="Generate positions of cells with deposited energy")
    jobTypeGroup.add_argument("--recSlidingWindow", action='store_true', help="Reconstruction with sliding window")
    jobTypeGroup.add_argument("--recTopoClusters", action='store_true', help="Reconstruction with topo-clusters")
    jobTypeGroup.add_argument("--estimatePileup", action='store_true', help="Estimate pileup from minbias events")
    jobTypeGroup.add_argument("--mergePileup", action='store_true', help="Estimate pileup from minbias events")
    jobTypeGroup.add_argument("--ntuple", action='store_true', help="Conversion to ntuple")
    jobTypeGroup.add_argument("--trackerPerformance", action='store_true', help="Tracker-only performance studies")
    # Add noise on cluster level
    parser.add_argument("--noise", action='store_true', help="Add electronics noise")
    parser.add_argument("--addPileupNoise", action='store_true', help="Add pile-up noise in qudrature to electronics noise")
    parser.add_argument('--pileup', type=int,  required = '--mergePileup' in sys.argv or '--addPileupNoise' in sys.argv, help='Pileup')
    parser.add_argument("--tripletTracker", action="store_true", help="Use triplet tracker layout instead of baseline")

    default_options,job_type,short_job_type,sim = getJobInfo(sys.argv)
    parser.add_argument('--jobOptions', type=str, default = default_options, help='Name of the job options run by FCCSW (default config/geantSim.py')

    genTypeGroup = parser.add_mutually_exclusive_group(required = True) # Type of events to generate
    genTypeGroup.add_argument("--singlePart", action='store_true', help="Single particle events")
    genTypeGroup.add_argument("--physics", action='store_true', help="Physics even run with Pythia")
    genTypeGroup.add_argument("--LHE", action='store_true', help="LHE events + Pythia")

    singlePartGroup = parser.add_argument_group('Single particles')
    singlePartGroup.add_argument('-e','--energy', type=int, required = '--singlePart' in sys.argv, help='Energy of particle in GeV')
    singlePartGroup.add_argument('--etaMin', type=float, default=0., help='Minimal pseudorapidity')
    singlePartGroup.add_argument('--etaMax', type=float, default=0., help='Maximal pseudorapidity')
    singlePartGroup.add_argument('--phiMin', type=float, default=0, help='Minimal azimuthal angle')
    singlePartGroup.add_argument('--phiMax', type=float, default=2*pi, help='Maximal azimuthal angle')
    singlePartGroup.add_argument('--particle', type=int,  required = '--singlePart' in sys.argv, help='Particle type (PDG)')

    physicsGroup = parser.add_argument_group('Physics','Physics process')
    lhe_physics = ["ljets","cjets","bjets","top","Wqq","Zqq","Hbb"]
    SM_physics = ["MinBias","Haa","Zee", "H4e"]
    physicsGroup.add_argument('--process', type=str, required = '--physics' in sys.argv, help='Process type', choices = lhe_physics + SM_physics)
    physicsGroup.add_argument('--pt', type=int,
                              required = ('ljets' in sys.argv or 'top' in sys.argv or 'Wqq' in sys.argv
                                          or 'Zqq' in sys.argv or 'Hbb' in sys.argv),
                              help='Transverse momentum of simulated jets (valid only for Hbb, Wqq, Zqq, t, ljet)')

    recoSlidingWinGroup = parser.add_argument_group('RecoSlidingWindow','Sliding window reconstruction')
    recoSlidingWinGroup.add_argument('--winEta', type=int, default=7, help='Size of the final cluster in eta')
    recoSlidingWinGroup.add_argument('--winPhi', type=int, default=19, help='Size of the final cluster in phi')
    recoSlidingWinGroup.add_argument('--enThreshold', type=float, default=3., help='Energy threshold of the seeding clusters [GeV]')

    recoTopoClusterGroup = parser.add_argument_group('RecoTopoCluster','Topological cluster reconstruction')
    recoTopoClusterGroup.add_argument('--sigma1', type=int, default=4, help='Energy threshold [in number of sigmas] for seeding')
    recoTopoClusterGroup.add_argument('--sigma2', type=int, default=2, help='Energy threshold [in number of sigmas] for neighbours')
    recoTopoClusterGroup.add_argument('--sigma3', type=int, default=0, help='Energy threshold [in number of sigmas] for last neighbours')

    args, _ = parser.parse_known_args()

    batchGroup = parser.add_mutually_exclusive_group(required = True) # Where to submit jobs
    batchGroup.add_argument("--lsf", action='store_true', help="Submit with LSF")
    batchGroup.add_argument("--condor", action='store_true', help="Submit with condor")
    batchGroup.add_argument("--no_submit", action='store_true', help="Just generate scripts without submitting")

    lsfGroup = parser.add_argument_group('Pythia','Common for min bias and LHE')
    lsfGroup.add_argument('-q', '--queue', type=str, default='8nh', help='lxbatch queue (default: 8nh)')

    args, _ = parser.parse_known_args()

    print "=================================="
    print "==      GENERAL SETTINGS       ==="
    print "=================================="

    if '--local' in sys.argv:
        print "FCCSW is taken from : ", args.local ,"\n"
        path=imp.load_source('path', args.local)
        try :
            path_to_INIT = path.path_to_INIT
        except AttributeError, e:
            pass
        try :
            path_to_FCCSW = path.path_to_FCCSW
        except AttributeError, e:
            pass

        print 'path_to_INIT : ',path_to_INIT
        print 'path_to_FCCSW: ',path_to_FCCSW

    version = args.version

    print 'FCCSim version: ',version
    magnetic_field = not args.bFieldOff
    b_field_str = "bFieldOn" if not args.bFieldOff else "bFieldOff"
    num_events = args.numEvents if sim else -1 # if reconstruction is done use -1 to run over all events in file
    num_jobs = args.numJobs
    job_options = args.jobOptions
    output_path = args.output

    ## Use Pileup valid only in certain cases!
    if args.addPileupNoise: # pileup is used to specify level of the pileup noise
        job_type = job_type.replace("pileupNoise", "pileupNoise/PU"+str(args.pileup))
        short_job_type += "PU"+str(args.pileup)
        num_events = args.numEvents
    elif args.mergePileup: # merge cells from simulated events
        job_type = "simuPU"+str(args.pileup)
        short_job_type += "PU"+str(args.pileup)
        num_events = args.numEvents
    elif (args.recPositions or args.recTopoClusters or args.recSlidingWindow) and args.pileup: # if reconstruction is run on pileup (mixed) events
        job_type = job_type.replace("ntup", "ntupPU"+str(args.pileup)).replace("reco", "recoPU"+str(args.pileup))
        short_job_type += "PU"+str(args.pileup)
    elif args.pileup:
        warning("'--pileup "+str(args.pileup)+"' is specified but no usecase was found. Please remove it or update job sending script.")
    if args.recSlidingWindow:
        job_type += '/eta' + str(args.winEta) + 'phi' + str(args.winPhi) + 'en' + str(args.enThreshold) + '/'
        short_job_type += '_eta' + str(args.winEta) + 'phi' + str(args.winPhi) + 'en' + str(args.enThreshold)
    elif args.recTopoClusters:
        job_type += '/' + str(args.sigma1) + str(args.sigma2) + str(args.sigma3) + '/'
        short_job_type += '_' + str(args.sigma1) + str(args.sigma2) + str(args.sigma3)

    print "B field: ", magnetic_field
    print "number of events = ", num_events
    print "output path: ", output_path
    print "FCCSW job options: ", job_options
    if args.singlePart:
        energy = args.energy
        etaMin = args.etaMin
        etaMax = args.etaMax
        phiMin = args.phiMin
        phiMax = args.phiMax
        pdg = args.particle
        particle_human_names = {11: 'electron', -11: 'positron', -13: 'mup', 13: 'mum', 22: 'photon', 111: 'pi0', 211: 'pip', -211: 'pim', 130: "K0L"}
        print "=================================="
        print "==       SINGLE PARTICLES      ==="
        print "=================================="
        print "particle PDG, name: ", pdg, " ", particle_human_names[pdg]
        print "energy: ", energy, "GeV"
        if etaMin == etaMax:
            print "eta: ", etaMin
            eta_str = "eta" + str(decimal.Decimal(str(etaMin)).normalize())
        else:
            print "eta: from ", etaMin, " to ", etaMax
            if abs(etaMin) == etaMax:
                eta_str = "etaTo" + str(decimal.Decimal(str(etaMax)).normalize())
            else:
                eta_str = "etaFrom" + str(decimal.Decimal(str(etaMin)).normalize()).replace('-','M') + "To" + str(decimal.Decimal(str(etaMax)).normalize()).replace('-','M')
        if phiMin == phiMax:
            print "phi: ", phiMin
            eta_str += "_phi" + str(decimal.Decimal(str(phiMin)).normalize())
        else:
            print "phi: from ", phiMin, " to ", phiMax
            if not(phiMin == 0 and phiMax == 2*pi):
                eta_str += "_phiFrom" + str(decimal.Decimal(str(phiMin)).normalize()).replace('-','M') + "To" + str(decimal.Decimal(str(phiMax)).normalize()).replace('-','M')
        job_dir = os.path.join("singlePart", particle_human_names[pdg], b_field_str, eta_str, str(energy) + "GeV")
    elif args.physics:
        print "=================================="
        print "==           PHYSICS           ==="
        print "=================================="
        process = args.process
        print "process: ", process
        job_dir = "physics/" + process + "/" + b_field_str + "/"
        if process in lhe_physics:
            LHE = True
            pt = args.pt
            print "LHE: ", LHE
            print "jet pt: ", pt
            job_dir += "etaTo0.5/" + str(pt) + "GeV/"
        else:
            LHE = False
            job_dir += "etaFull/"
        print "=================================="
    thebatch=''
    if args.condor:
        print "submitting to condor ..."
        thebatch='condor'
    elif args.lsf:
        print "submitting to lxbatch ..."
        queue = args.queue
        print "queue: ", queue
        thebatch='lsf'

    rundir = os.getcwd()
    nbjobsSub=0

    if args.mergePileup and not args.local == "inits/pileup.py":
        warning("Please note that '--mergePileup' is not supported for FCCSW v0.9.1. Make sure that you use suitable software version (recommended: '--local inits/pileup.py')", True)
    if args.estimatePileup and not args.local == "inits/pileup.py":
        warning("Please note that '--preparePileup' is not supported for FCCSW v0.9.1. Make sure that you use suitable software version (recommended: '--local inits/pileup.py')", True)
    if args.recPositions and not args.local == "inits/reco.py":
        warning("Please note that '--recPositions' is not supported for FCCSW v0.9.1. Make sure that you use suitable software version (recommended: '--local inits/reco.py')", True)
    if args.recTopoClusters and not args.local == "inits/reco.py":
        warning("Please note that '--recTopoClusters' is not supported for FCCSW v0.9.1. Make sure that you use suitable software version (recommended: '--local inits/reco.py')", True)

    # first make sure the output path for root files exists
    outdir = os.path.join( output_path, version, job_dir, job_type)
    print "Output will be stored in ... ", outdir
    if not sim:
        # if signal events are used (simu/) or mixed pileup events (simuPU.../)
        if not args.mergePileup and not args.addPileupNoise and args.pileup and not (args.pileup == 0):
            inputID = os.path.join(yamldir, version, job_dir, 'simuPU'+str(args.pileup))
        else:
            inputID = os.path.join(yamldir, version, job_dir, 'simu')
        outputID = os.path.join(yamldir, version, job_dir, job_type)

        input_files = getInputFiles(inputID)
        input_files, instatus = takeOnlyNonexistingFiles(input_files, outputID)

        if instatus == 0:
            warning("Directory contains no files")
            exit()
        if instatus < num_jobs:
            num_jobs = instatus
            warning("Directory contains only " + str(instatus) + " files, using all for the reconstruction")
    # for the purpose of mixing pileup events all inputs must be passed to the config
    # merging pileup events will be done randomly from a given event pool
    if args.mergePileup:
        all_inputs = ""
        for f in input_files:
            all_inputs += " " + f # event pool = all inputs
        seed=ut.getuid() # to generate new output name
        print 'seed  ',seed
        outfile = 'output_%i.root'%(seed)
        print "Name of the output file: ", outfile

    os.system("mkdir -p %s"%outdir)
    # first make sure the output path exists
    current_dir = os.getcwd()
    logdir = args.log if os.path.isabs(args.log) else current_dir + '/' +  args.log
    job_options = job_options if os.path.isabs(job_options) else current_dir + '/' +  job_options
    print  "Saving the logs in ... ",logdir

    # if Pythia used, figure out the Pythia card
    if sim and args.physics:
        if LHE:
            card = current_dir + "/PythiaCards/pythia_pp_LHE.cmd"
        else:
            card = current_dir + "/PythiaCards/pythia_pp_" + process + ".cmd"
        print card
        print os.path.isfile(card)

    seed=''
    uid=os.path.join(version, job_dir, job_type)
    for i in xrange(num_jobs):
        if sim:
            seed=ut.getuid()
            yamldir_process = '%s/%s'%(yamldir,uid)
            if not ut.dir_exist(yamldir_process):
                os.system("mkdir -p %s"%yamldir_process)
            myyaml = my.makeyaml(yamldir_process, seed)
            if not myyaml:
                print 'job %s already exists'%str(seed)
                continue

            print 'seed  ',seed
            outfile = 'output_%i.root'%(seed)
            print "Name of the output file: ", outfile
        else:
            infile = os.path.basename(input_files[i])
            outfile = infile
            print "Name of the input file: ", infile

            yamldir_process = '%s/%s'%(yamldir,uid)
            if not ut.dir_exist(yamldir_process):
                os.system("mkdir -p %s"%yamldir_process)

            infile_split = re.split(r'[_.]',infile)
            seed = infile_split[1]
            myyaml = my.makeyaml(yamldir_process, seed)
            if not myyaml:
                print 'job %s already exists'%str(seed)
                continue


        if args.physics:
            frunname = 'job_%s_%s_%s.sh'%(short_job_type,process,str(seed))
        else:
            frunname = 'job_%s_%s_%dGeV_eta%s_%s.sh'%(short_job_type,particle_human_names[pdg],energy,str(decimal.Decimal(str(etaMax)).normalize()),str(seed))

        os.system("mkdir -p %s"%logdir)
        frun = None
        try:
            frun = open(logdir+'/'+frunname, 'w')
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
            time.sleep(10)
            frun = open(logdir+'/'+frunname, 'w')
        print frun


        commands.getstatusoutput('chmod 777 %s/%s'%(logdir,frunname))
        frun.write('#!/bin/bash\n')
        frun.write('unset LD_LIBRARY_PATH\n')
        frun.write('unset PYTHONHOME\n')
        frun.write('unset PYTHONPATH\n')
        frun.write('export JOBDIR=$PWD\n')
        frun.write('source %s\n' % (path_to_INIT))
        # workaround for pythia 8 version mismatch
        frun.write("export PYTHIA8_XML=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/pythia8/230-b1563/x86_64-slc6-gcc62-opt/share/Pythia8/xmldoc/\n")
        frun.write("export PYTHIA8DATA=$PYTHIA8_XML\n")

        # set options to run FCCSW
        common_fccsw_command = '%s/run fccrun.py %s --outName $JOBDIR/%s --numEvents %i'%(path_to_FCCSW,job_options, outfile ,num_events)
        if not magnetic_field:
            common_fccsw_command += ' --bFieldOff'
        if sim:
            common_fccsw_command += ' --seed %i'%(seed)
        if args.noise:
            common_fccsw_command += ' --addElectronicsNoise'
        if args.addPileupNoise:
            common_fccsw_command += ' --addPileupNoise --mu ' + str(args.pileup)
        if args.physics:
            common_fccsw_command += ' --physics'
        if '--local' in sys.argv:
            common_fccsw_command += ' --detectorPath ' + path_to_FCCSW
        if args.mergePileup:
            common_fccsw_command += ' --pileup ' + str(args.pileup)
        if args.recPositions and args.pileup:
            common_fccsw_command += ' --prefixCollections merged '
        if args.recSlidingWindow:
            common_fccsw_command += ' --winEta ' + str(args.winEta) + ' --winPhi ' + str(args.winPhi) + ' --enThreshold ' + str(args.enThreshold) + ' '
        if args.recTopoClusters:
            common_fccsw_command += ' --sigma1 ' + str(args.sigma1) + ' --sigma2 ' + str(args.sigma2) + ' --sigma3 ' + str(args.sigma3) + ' '
        print '-------------------------------------'
        print common_fccsw_command
        print '-------------------------------------'
        if args.tripletTracker:
          common_fccsw_command += '--tripletTracker'
        if sim:
            if args.physics:
                frun.write('cp %s $JOBDIR/card.cmd\n'%(card))
                if LHE:
                    #CLEMENT TO BE CHANGED the eta values etaMin,etaMax -> -1.6,1.6
                    frun.write('cd %s\n' %(path_to_LHE))
                    if process=='ljets':
                        frun.write('python flatGunLHEventProducer.py   --pdg 1 2 3   --guntype pt   --nevts %i   --ecm 100000.   --pmin %f   --pmax %f   --etamin %f   --etamax %f  --seed %i   --output $JOBDIR/events.lhe\n'%(num_events,pt,pt,-0.5,0.5,seed))
                    elif process=='cjets':
                        frun.write('python flatGunLHEventProducer.py   --pdg 4   --guntype pt   --nevts %i   --ecm 100000.   --pmin %f   --pmax %f   --etamin %f   --etamax %f  --seed %i   --output $JOBDIR/events.lhe\n'%(num_events,pt,pt,-1.6,1.6,seed))
                    elif process=='bjets':
                        frun.write('python flatGunLHEventProducer.py   --pdg 5   --guntype pt   --nevts %i   --ecm 100000.   --pmin %f   --pmax %f   --etamin %f   --etamax %f  --seed %i   --output $JOBDIR/events.lhe\n'%(num_events,pt,pt,-1.6,1.6,seed))
                    elif process=='top':
                        frun.write('python flatGunLHEventProducer.py   --pdg 6   --guntype pt   --nevts %i   --ecm 100000.   --pmin %f   --pmax %f   --etamin %f   --etamax %f  --seed %i   --output $JOBDIR/events.lhe\n'%(num_events,pt,pt,-1.6,1.6,seed))
                    elif process=='Zqq':
                        frun.write('python flatGunLHEventProducer.py   --pdg 23   --guntype pt   --nevts %i   --ecm 100000.   --pmin %f   --pmax %f   --etamin %f   --etamax %f  --seed %i   --output $JOBDIR/events.lhe\n'%(num_events,pt,pt,-0.5,0.5,seed))
                    elif process=='Wqq':
                        frun.write('python flatGunLHEventProducer.py   --pdg 24   --guntype pt   --nevts %i   --ecm 100000.   --pmin %f   --pmax %f   --etamin %f   --etamax %f  --seed %i   --output $JOBDIR/events.lhe\n'%(num_events,pt,pt,-0.5,0.5,seed))
                    elif process=='Hbb':
                        frun.write('python flatGunLHEventProducer.py   --pdg 25   --guntype pt   --nevts %i   --ecm 100000.   --pmin %f   --pmax %f   --etamin %f   --etamax %f  --seed %i   --output $JOBDIR/events.lhe\n'%(num_events,pt,pt,-1.6,1.6,seed))
                    else:
                        print 'process does not exists, exit'
                        sys.exit(3)
                    frun.write('gunzip $JOBDIR/events.lhe.gz\n')
                    frun.write('echo "Beams:LHEF = $JOBDIR/events.lhe" >> $JOBDIR/card.cmd\n')
                # run FCCSW using Pythia as generator
                frun.write('cd $JOBDIR\n')
                if args.pythiaSmearVertex:
                  frun.write('echo "Beams:allowVertexSpread = on" >> $JOBDIR/card.cmd\n')
                  frun.write('echo "Beams:sigmaVertexX = 0.5" >> $JOBDIR/card.cmd\n')
                  frun.write('echo "Beams:sigmaVertexY = 0.5" >> $JOBDIR/card.cmd\n')
                  frun.write('echo "Beams:sigmaVertexZ = 40.0" >> $JOBDIR/card.cmd\n')
                  frun.write('echo "Beams:sigmaTime = 0.180" >> $JOBDIR/card.cmd\n')
                frun.write('%s  --pythia --card $JOBDIR/card.cmd \n'%(common_fccsw_command))
            else:
                # run single particles
                frun.write('cd $JOBDIR\n')
                frun.write('%s  --singlePart --particle %i -e %i --etaMin %f --etaMax %f --phiMin %f --phiMax %f\n'%(common_fccsw_command, pdg, energy, etaMin, etaMax, phiMin, phiMax))
        else:
            frun.write('cd $JOBDIR\n')
            if args.mergePileup:
                frun.write('%s --inName %s\n'%(common_fccsw_command, all_inputs))
            elif not(args.pionRejection):
                frun.write('%s --inName %s\n'%(common_fccsw_command, input_files[i]))
        if args.recPositions:
            frun.write('python %s/python/Convert.py edm.root $JOBDIR/%s\n'%(current_dir,outfile))
            frun.write('rm edm.root \n')
        elif args.pionRejection:
            frun.write('python %s/python/Pionrejection_variables.py %s $JOBDIR/%s\n'%(current_dir,input_files[i], outfile))
        elif args.recTopoClusters or args.recSlidingWindow:
            frun.write('python %s/python/Convert.py $JOBDIR/%s $JOBDIR/clusters.root\n'%(current_dir,outfile))
            ntup_path = outdir.replace('/reco', '/ntup')
            if not ut.dir_exist(ntup_path):
                os.system("mkdir -p %s"%(ntup_path))
            frun.write('python /afs/cern.ch/work/h/helsens/public/FCCutils/eoscopy.py $JOBDIR/clusters.root %s/%s\n'%(ntup_path, outfile))
            if args.recSlidingWindow:
                analysis_path = outdir.replace('/reco', '/ana')
                if not ut.dir_exist(analysis_path):
                    os.system("mkdir -p %s"%(analysis_path))
                frun.write('python /afs/cern.ch/work/h/helsens/public/FCCutils/eoscopy.py $JOBDIR/hist_%s %s\n'%(outfile, analysis_path))

        if not args.no_eoscopy:
          frun.write('python /afs/cern.ch/work/h/helsens/public/FCCutils/eoscopy.py $JOBDIR/%s %s\n'%(outfile,outdir))
          frun.write('rm $JOBDIR/%s \n'%(outfile))
        frun.close()

        if args.lsf:
            #cmdBatch="bsub -M 4000000 -R \"pool=40000\" -q %s -cwd%s %s" %(queue, logdir,logdir+'/'+frunname)
            cmdBatch="bsub  -R 'rusage[mem=16000:pool=8000]' -q %s -cwd%s %s" %(queue, logdir,logdir+'/'+frunname)
            batchid=-1
            job,batchid=ut.SubmitToLsf(cmdBatch,10)
        elif args.no_submit:
            job = 0
            print "scripts generated in ", os.path.join(logdir, frunname),
        else:
            os.system("mkdir -p %s/out"%logdir)
            os.system("mkdir -p %s/log"%logdir)
            os.system("mkdir -p %s/err"%logdir)
            # create also .sub file here
            fsubname = frunname.replace('.sh','.sub')
            fsub = None
            try:
                fsub = open(logdir+'/'+fsubname, 'w')
            except IOError as e:
                print "I/O error({0}): {1}".format(e.errno, e.strerror)
                time.sleep(10)
                fsub = open(logdir+'/'+fsubname, 'w')
            print fsub
            commands.getstatusoutput('chmod 777 %s/%s'%(logdir,fsubname))
            fsub.write('executable            = %s/%s\n' %(logdir,frunname))
            fsub.write('arguments             = $(ClusterID) $(ProcId)\n')
            fsub.write('output                = %s/out/job.%s.$(ClusterId).$(ProcId).out\n'%(logdir,str(seed)))
            fsub.write('log                   = %s/log/job.%s.$(ClusterId).log\n'%(logdir,str(seed)))
            fsub.write('error                 = %s/err/job.%s.$(ClusterId).$(ProcId).err\n'%(logdir,str(seed)))
            if args.physics:
                fsub.write('RequestCpus = 8\n')
            else:
                fsub.write('RequestCpus = 4\n')
            fsub.write('+JobFlavour = "nextweek"\n')
            fsub.write('+AccountingGroup = "group_u_FCC.local_gen"\n')
            fsub.write('queue 1\n')
            fsub.close()
            cmdBatch="condor_submit %s/%s"%(logdir.replace(current_dir+"/",''),fsubname)

            print cmdBatch
            batchid=-1
            job,batchid=ut.SubmitToCondor(cmdBatch,10)

        nbjobsSub+=job

    print 'succesfully sent %i  jobs'%nbjobsSub
