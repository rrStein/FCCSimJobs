import os
import sys
import math as m
import numpy as n
import ROOT as r
import os.path
import os
from scipy.signal import argrelextrema, savgol_filter, find_peaks_cwt, argrelmax

from ROOT import gSystem
result=gSystem.Load("libDDCorePlugins")
from ROOT import dd4hep
if result < 0:
    print "No lib loadable!"

system_decoder = dd4hep.DDSegmentation.BitFieldCoder("system:4")
ecalBarrel_decoder = dd4hep.DDSegmentation.BitFieldCoder("system:4,cryo:1,type:3,subtype:3,layer:8,eta:9,phi:10")
hcalBarrel_decoder = dd4hep.DDSegmentation.BitFieldCoder("system:4,module:8,row:9,layer:5")
hcalExtBarrel_decoder = dd4hep.DDSegmentation.BitFieldCoder("system:4,module:8,row:9,layer:5")
ecalEndcap_decoder = dd4hep.DDSegmentation.BitFieldCoder("system:4,subsystem:1,type:3,subtype:3,layer:8,eta:10,phi:10")
hcalEndcap_decoder = dd4hep.DDSegmentation.BitFieldCoder("system:4,subsystem:1,type:3,subtype:3,layer:8,eta:10,phi:10")
ecalFwd_decoder = dd4hep.DDSegmentation.BitFieldCoder("system:4,subsystem:1,type:3,subtype:3,layer:8,eta:11,phi:10")
hcalFwd_decoder = dd4hep.DDSegmentation.BitFieldCoder("system:4,subsystem:1,type:3,subtype:3,layer:8,eta:11,phi:10")
trackerBarrel_decoder = dd4hep.DDSegmentation.BitFieldCoder("system:4,layer:5,module:18,x:-15,z:-15")
trackerEndcap_decoder = dd4hep.DDSegmentation.BitFieldCoder("system:4,posneg:1,disc:5,component:17,x:-15,z:-15")

lastECalBarrelLayer = int(7)
lastECalEndcapLayer = int(39)
lastECalFwdLayer = int(41)

def systemID(cellid):
    return system_decoder.get(cellid, "system")

def benchmarkCorr(ecal, ecal_last, ehad, ehad_first):
    a=0.978
    b=0.479
    c=-0.0000054
    ebench = ecal*a + ehad + b * math.sqrt(math.fabs(a*ecal_last*ehad_first)) + c*(ecal*a)**2
    return ebench

def signy(y):
    if y>0: return 1
    elif y<0: return -1
    return 0

# Function for calculating the shower width on j cells in cell units. imax describes the cells
# cell with maximum energy deposit.

def Shower_width(Energies, Cells, j):

    top = 0
    bot = 0
    k = 0
    imaxphi = int(Cells[n.argmax(Energies)][0])
    imaxeta = int(Cells[n.argmax(Energies)][1])

    # Loop over j total selected cells.
    while k < j:
        # Loop over all detected cells.
        for i in Cells:
            # In case no cells were found the value will be 0.
            try:
                # Condition to choose only cells within j cell units of imax and sum the recorded
                # energies multiplied by the difference in cell and imax squared.
                if int(i[0]) > (imaxphi - j) and int(i[0]) < (imaxphi + j):
                    top += sum(Energies[n.where(Cells == i)[0]])*(((int(i[0])-imaxphi)**2))
                    bot += sum(Energies[n.where(Cells == i)[0]])

                    k += 1

                else:
                    top += 0
                    bot += 0

            except(IndexError):
                top += 0
                bot += 0


    # Divide the sum of energies times the difference in cell and imax squared by the sum of Energies
    # and take the square root for the shower width.
	Wnst = m.sqrt((top/bot))
    return Wnst

# Function for finding the difference of the energy in the cell with the second maximum and the
# energy in the cell with the minimal value between the first and second maximum.
def edmaxy(Emax,E2ndmax,Cells,Energies):

# Define the cells with the 1st and 2nd maximum in terms of eta and phi.
    try:

        cell2max = Cells[n.where(Energies == E2ndmax)[0][0]]
        cell2maxphi = int(cell2max[0])
        cell2maxeta = int(cell2max[1])
        cellmax = Cells[n.where(Energies == Emax)[0][0]]
        cellmaxphi = int(cellmax[0])
        cellmaxeta = int(cellmax[1])

    except(IndexError):

        print "indexerror"
        return -1

    cellminE = []
    i = 0
    if m.fabs(cellmaxphi-cell2maxphi) <= 1:
        return 0.

    else:
        # Loop over all the cells recorded.
        while i < len(Cells):

            # Condition to choose only cells that are on the same phi and eta values as the 1st and
            # 2nd maximum or inbetween these phi and eta values. Record all the energies in these cells.
            if (int(Cells[i][0]) > cell2maxphi and int(Cells[i][0]) < cellmaxphi) or (int(Cells[i][0]) > cellmaxphi and int(Cells[i][0]) < cell2maxphi):
                if (int(Cells[i][1]) >= cell2maxeta and int(Cells[i][1]) <= cellmaxeta) or (int(Cells[i][1]) >= cellmaxeta and int(Cells[i][1]) <= cell2maxeta):
                    cellminE.append(Energies[i])

            i += 1

        # If there were no cells matching the condition then return 0 for edmax.
        if cellminE == []:
            return 0.

        # Otherwise take the minimum value of all the recorded energies in the cells matching the condition.
        else:
            ncellminE = n.array(cellminE)
            minE = n.amin(ncellminE)
            mincellphi = Cells[n.where(Energies == minE)[0][0]][0]
            mincellE = 0.
            E2cellE = 0.
            for c in Cells:
                if c[0] == cell2maxphi:
                    E2cellE += Energies[n.where(Cells == c)[0][0]]
                if c[0] == mincellphi:
                    mincellE += Energies[n.where(Cells == c)[0][0]]
                else:
                    continue

        # Return the difference of the 2nd maximum and the minimum in the valley in between 1st and 2nd maximum.
        edmax = m.sqrt((E2cellE - mincellE)**2)
        return edmax

# Function to find the energy deposited outside of the shower core.
def eocorey(Emax,Cells,Energies,j):

    # Define the cell with maximum energy deposit in terms of phi and eta.
    cellmaxphi = int(Cells[n.argmax(Energies)][0])
    cellmaxeta = int(Cells[n.argmax(Energies)][1])
    E3 = 0.
    E1 = 0.
    i = 0

    # Loop over all recorded cells.
    while i < len(Cells):

        # Choose only cells that are within +/- j cells around (and including) the cell with maximum energy deposit
        # and sum the energies deposited in these cells.
        if int(Cells[i][0])  <= cellmaxphi + j and int(Cells[i][0]) >= cellmaxphi - j:
            E3 += Energies[i]

            # Choose only cells that are within +/- 1 cells around (and including) the cell with maximum energy deposit
            # and sum the energies deposited in these cells.
            if int(Cells[i][0]) <= cellmaxphi + 1 and int(Cells[i][0]) >= cellmaxphi - 1:
                E1 += Energies[i]
        i += 1

    eocore = (E3 - E1) / E1
    return eocore

# Define the variables that are to be found/converted.
ev_num = n.zeros(1, dtype=int)
ev_nRechits = n.zeros(1, dtype=int)
e2max = n.zeros(1, dtype=float)
emax = n.zeros(1, dtype=float)
edmax = n.zeros(1, dtype=float)
eocore = n.zeros(1, dtype=float)
w3st = n.zeros(1, dtype=float)
w21st = n.zeros(1, dtype=float)

if len(sys.argv)!=3:
    print 'usage python Convert.py infile outfile'
infile_name = sys.argv[1]
outfile_name = sys.argv[2]

current_dir = os.getcwd()

if os.path.isfile(outfile_name) == False:
    infile=r.TFile.Open(infile_name)
    intree=infile.Get('events')

    maxEvent = intree.GetEntries()
    print 'Number of events : ',maxEvent

    outfile=r.TFile(outfile_name,"recreate")
    outtree=r.TTree('events','Events')

    # Branches for the discriminating variables of the ecal detector.

    outtree.Branch("e2max", e2max, "e2max/D")
    outtree.Branch("emax", emax, "emax/D")
    outtree.Branch("edmax", edmax, "edmax/D")
    outtree.Branch("eocore", eocore, "eocore/D")
    outtree.Branch("w3st", w3st, "w3st/D")
    outtree.Branch("w21st", w21st, "w21st/D")

# If the output root file already exists, update it.
else:
    infile=r.TFile.Open(infile_name)
    intree=infile.Get('events')
    outfile=r.TFile.Open(outfile_name,"update")
    outtree=outfile.Get('events')

    maxEvent = intree.GetEntries()
    print 'Number of events : ',maxEvent

    # Branches for the discriminating variables of the ecal detector.

    outtree.SetBranchAddress("e2max", e2max)
    outtree.SetBranchAddress("emax", emax)
    outtree.SetBranchAddress("edmax", edmax)
    outtree.SetBranchAddress("eocore", eocore)
    outtree.SetBranchAddress("w3st", w3st)
    outtree.SetBranchAddress("w21st", w21st)

# Loop over all events in the input file.
numEvent = 0
for event in intree:
    ev_num[0] = numEvent
    numHits = 0
    E = .0
    calE = 0
    cal1E = 0
    Ecal_Phi = []
    Ecal_cell = []
    Ecal1_E = []
    Ecal_Eta = []
    cal1Emax = -1.
    cal1E2max = -1.
    cal1Etamax = -1.
    cal1Eta2max = -1.
    cal1Phimax = -1.
    cal1Phi2max = -1.

    # Loop over everything recorded in the ecal barrel cells for each event.
    for c in event.ECalBarrelCells:
        Layer = ecalBarrel_decoder.get(c.core.cellId, "layer")
        # Choose only the second layer of the barrel.
        if Layer == 1 or Layer == 0:

            # Record all energies and cells with their phi and eta values in the layer.
            Eta = ecalBarrel_decoder.get(c.core.cellId, "eta")
            cal1E = c.core.energy
            #eta = int(Eta) - 169*0.01
            Phi = ecalBarrel_decoder.get(c.core.cellId, "phi")
            Ecal_Eta.append(Eta)
            Ecal_cell.append([Phi,Eta])
            Ecal1_E.append(cal1E)
            Ecal_Phi.append(Phi)

            # If the energy found is larger than the 1st or 2nd maximum defined in previous iteration
            # redefine them accordingly
            if cal1E >= cal1Emax or cal1E > cal1E2max:

                # If it is the first iteration just set the values as the 1st maximum
                if len(Ecal1_E) < 2:

                    cal1Emax = cal1E
                    cal1Etamax = Eta
                    cal1Phimax = Phi

                # If the recorded energy is greater than the previous 1st maximum...
                elif cal1E > cal1Emax:
                    # ...and not in the same phi slice, then set the 2nd maximum as the previous 1st maximum
                    # and update the 1st maximum.
                    if Phi != Ecal_Phi[Ecal1_E.index(cal1Emax)]:
                        #if m.fabs(Eta - Ecal_Eta[Ecal1_E.index(cal1E)]) <= 1:
                        cal1E2max = cal1Emax
                        cal1Eta2max = cal1Etamax
                        cal1Phi2max = cal1Phimax
                        cal1Emax = cal1E
                        cal1Etamax = Eta
                        cal1Phimax = Phi

                    # If it is in the same phi slice, then set it as a new 1st maximum without updating the 2nd
                    # maximum.
                    else:

                        cal1Emax = cal1E
                        cal1Etamax = Eta
                        Cal1Phimax = Phi

                # If it is not larger than the 1st maximum then it must be larger than the 2nd maximum and
                # in this case update the 2nd maximum only.
                else:

                    cal1E2max = cal1E
                    cal1Eta2max = Eta
                    cal1Phi2max = Phi


        numHits += 1

    # Convert the lists into numpy arrays.
    Cellids = n.array(Ecal_cell)
    Energies = n.array(Ecal1_E)
    Etas = n.array(Ecal_Eta)
    Phis = n.array(Ecal_Phi)
    Es_eta0 = []
    Es_eta1 = []
    Es_eta_1 = []

    for e in Energies:
        #if Phis[n.where(Energies == e)[0][0]] != Phis[n.where(Energies == cal1Emax)[0][0]]:
        if Etas[n.where(Energies == e)[0][0]] == cal1Etamax:
            Es_eta0.append(e)
        elif (Etas[n.where(Energies == e)[0][0]] - cal1Etamax) == 1:
            Es_eta1.append(e)
        elif (Etas[n.where(Energies == e)[0][0]] - cal1Etamax) == -1:
            Es_eta_1.append(e)
            #    else:
        #        continue
        else:
            continue

    NewE_eta0 = n.array(Es_eta0)
    NewE_eta1 = n.array(Es_eta1)
    NewE_eta_1 = n.array(Es_eta_1)
    #print n.sort(NewEnergies)

    '''if len(Energies) > 17:
        smoothenergies = savgol_filter(Energies,17,4)
    #    smoothenergies = savgol_filter(Energies,17,4)
    elif len(Energies)%2 == 0:
        smoothenergies = savgol_filter(Energies,len(Energies)-1,4)
    else:
        smoothenergies = savgol_filter(Energies,len(Energies),4)

    print len(smoothenergies), len(Energies)'''
    # If there were more than 2 Energies recorded in the barrel, then find the extremum points of
    # these energies and sort the list of energies at these extremum indices.
    try:
        MaxInd0 = argrelextrema(NewE_eta0, n.greater)
        MaxInd1 = argrelextrema(NewE_eta1, n.greater)
        MaxInd_1 = argrelextrema(NewE_eta_1, n.greater)
        #MaxInd2 = find_peaks_cwt(NewEnergies,n.arange(1,5))

        Maxes0 = n.sort(NewE_eta0[MaxInd0])
        Maxes1 = n.sort(NewE_eta1[MaxInd1])
        Maxes_1 = n.sort(NewE_eta_1[MaxInd_1])
        #Maxes2 = n.sort(NewEnergies[MaxInd2])
        #print Maxes
        #print Maxes1
        #print cal1Emax, cal1E2max, Maxes[len(Maxes)-1]
        #print n.amax(Maxes)
    except(IndexError):
        print "no maxima"
    # If there are any extremums, then check if the largest extremum is bigger than the 1st maximum
    # recorded earlier. If it is, update the 2nd and 1st maxima.
    try:
        #print n.float(Maxes0[len(Maxes0)-1]), cal1Emax
        #print n.float(Maxes0[len(Maxes0)-1]), cal1E2max, "\n"
        '''if n.float(Maxes[len(Maxes)-1]) > cal1Emax:
            cal1E2max = cal1Emax
            cal1Phi2max = cal1Phimax
            cal1Eta2max = cal1Etamax
            cal1Emax = Maxes[len(Maxes)-1]
            cal1Phimax = Phis[n.where(Energies == cal1Emax)[0][0]]
            cal1Etamax = Etas[n.where(Energies == cal1Emax)[0][0]]'''

        # If the maximal extremum energy is less than the 1st maximum but greater than the 2nd
        # maximum, set it as the 2nd maximum.
        if n.float(Maxes0[len(Maxes0)-1]) < cal1Emax:
            if n.float(Maxes0[len(Maxes0)-1]) > cal1E2max:
                #print "found"
                cal1E2max = Maxes0[len(Maxes0)-1]
                cal1Phi2max = Phis[n.where(Energies == cal1E2max)[0][0]]
                cal1Eta2max = Etas[n.where(Energies == cal1E2max)[0][0]]

        # In all other cases, set the 2nd maximum as the second largest extremum value found
        # to make sure it is a local maxima not a part of the 1st maximum in another cell.
        if n.float(Maxes0[len(Maxes0)-1]) == cal1Emax:
            cal1E2max = Maxes0[(len(Maxes0)-2)]
            cal1Phi2max = Phis[n.where(Energies == cal1E2max)[0][0]]
            cal1Eta2max = Etas[n.where(Energies == cal1E2max)[0][0]]
        #print cal1E2max
        #print "\n"
    except(IndexError):
        cal1E2max = 0.
        #print "2nd max is 0"
    #except(ValueError):
    #    continue
    #print cal1Emax, cal1E2max
   # if cal1Emax == 0. or cal1E2max == 0:
       # #print cal1Emax, cal1E2max, Energies
    ##print "All recorded energies: ", Energies
    ##print Cellids, "\n", n.unique(Cellids)
    #print "No. of energies: ", len(Energies), "No. of cells: ", len(Phis), "No. of unique cells: ", len(n.unique(Phis))
    #print "1st maximum: ", cal1Emax, "and cell: ", cal1Phimax
    #print "2nd maximum: ", cal1E2max, "and cell: ",  cal1Phi2max
    '''#print "The energies using numpy.where: 1st max: ", Energies[n.where(Energies == cal1Emax)[0][0]], "2nd max: ",Energies[n.where(Energies == cal1E2max)[0][0]]
    #print "The cells using numpy.where: 1st max cell: ", Cellids[n.where(Energies == cal1Emax)[0][0]], "2nd max cell: ", Cellids[n.where(Energies == cal1E2max)[0][0]]
    ev_nRechits[0] = numHits
    #print "\n"
    #print "test n.where: ", Cellids[n.where(Energies == cal1Emax)[0]], Energies[n.where(Energies == cal1Emax)[0]]'''
    #print "\n"

    # If there were more than 2 energies recorded in the barrel, then set calculate the variables and fill
    # them in the branches.
    if len(Energies) > 2:
        w3st[0] = Shower_width(Energies, Cellids, 3)
        w21st[0] = Shower_width(Energies, Cellids, 21)
        eocore[0] = eocorey(cal1Emax, Cellids, Energies, 3)
        e2max[0] = cal1E2max
        emax[0] = cal1Emax
        edmax[0] = edmaxy(cal1Emax, cal1E2max, Cellids, Energies)
    else:
        continue

    outtree.Fill()

    numEvent += 1

outtree.Write()
outfile.Write()
outfile.Close()
