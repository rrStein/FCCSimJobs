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

    try:
        top = 0
        bot = 0
        k = 0
        imaxphi = int(Cells[n.argmax(Energies)][0])
        imaxeta = int(Cells[n.argmax(Energies)][1])

        # Loop over j total selected cells.
        while k < j:
            # Loop over all detected cells.
            for cell in Cells:
                # In case no cells were found the value will be 0.
                try:
                    # Condition to choose only cells within j cell units of imax and sum the recorded
                    # energies multiplied by the difference in cell and imax squared.
                    if int(cell[0]) >= (imaxphi - j) and int(cell[0]) <= (imaxphi + j):
                        top += sum(Energies[n.where(Cells == cell)[0]])*(((int(cell[0])-imaxphi)**2))
                        bot += sum(Energies[n.where(Cells == cell)[0]])

                        k += 1

                    else:
                        top += 0
                        bot += 0

                except(IndexError):
                    top += 0
                    bot += 0


        # Divide the sum of energies times the difference in cell and imax squared by the sum of Energies
        # and take the square root for the shower width.
    	Wnst = top/bot
        return Wnst
    except(IndexError, ValueError):
        return 0.

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

        cellminE = []
        i = 0
        if m.fabs(cellmaxphi-cell2maxphi) <= 1:
            return 0.

        else:
            # Loop over all the cells recorded.
            while i < len(Cells):

                # Condition to choose only cells that are on the same phi and eta values as the 1st and
                # 2nd maximum or inbetween these phi and eta values. Record all the energies in these cells.
                if (int(Cells[i][0]) >= cell2maxphi and int(Cells[i][0]) <= cellmaxphi) or (int(Cells[i][0]) >= cellmaxphi and int(Cells[i][0]) <= cell2maxphi):
                    if (int(Cells[i][1]) >= cell2maxeta and int(Cells[i][1]) <= cellmaxeta) or (int(Cells[i][1]) >= cellmaxeta and int(Cells[i][1]) <= cell2maxeta):
                        if Energies[i] != Emax:
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
            edmax = E2ndmax - minE
            edmax_2 = m.sqrt((E2cellE - mincellE)**2)
            return edmax
    except(IndexError):
        #print "indexerror"
        return 0.

# Function to find the energy deposited outside of the shower core.
def eocorey(Emax,Cells,Energies,j):

    # Define the cell with maximum energy deposit in terms of phi and eta.
    try:
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
    except(IndexError,ValueError):
        return 0.


# Define the variables that are to be found/converted.
ev_num = n.zeros(1, dtype=int)
ev_nRechits = n.zeros(1, dtype=int)

# Variables for first layer.
e2max_l00 = n.zeros(1, dtype=float)
emax_l00 = n.zeros(1, dtype=float)
edmax_l00 = n.zeros(1, dtype=float)
eocore_l00 = n.zeros(1, dtype=float)
w3st_l00 = n.zeros(1, dtype=float)
w21st_l00 = n.zeros(1, dtype=float)

# Variables for second layer.
e2max_l01 = n.zeros(1, dtype=float)
emax_l01 = n.zeros(1, dtype=float)
edmax_l01 = n.zeros(1, dtype=float)
eocore_l01 = n.zeros(1, dtype=float)
w3st_l01 = n.zeros(1, dtype=float)
w21st_l01 = n.zeros(1, dtype=float)

# Variables for layer energies against first layer energy and total energy
e_l10 = n.zeros(1, dtype=float)
e_l20 = n.zeros(1, dtype=float)
e_l30 = n.zeros(1, dtype=float)
e_l40 = n.zeros(1, dtype=float)
e_l50 = n.zeros(1, dtype=float)
e_l60 = n.zeros(1, dtype=float)
e_l70 = n.zeros(1, dtype=float)

e_l0T = n.zeros(1, dtype=float)
e_l1T = n.zeros(1, dtype=float)
e_l2T = n.zeros(1, dtype=float)
e_l3T = n.zeros(1, dtype=float)
e_l4T = n.zeros(1, dtype=float)
e_l5T = n.zeros(1, dtype=float)
e_l6T = n.zeros(1, dtype=float)
e_l7T = n.zeros(1, dtype=float)


if len(sys.argv)!=3:
    print 'usage python Convert.py infile outfile'
infile_name = sys.argv[1]
outfile_name = sys.argv[2]

current_dir = os.getcwd()

#if os.path.isfile(outfile_name) == False:
infile=r.TFile.Open(infile_name)
intree=infile.Get('events')

maxEvent = intree.GetEntries()
print 'Number of events : ',maxEvent

outfile=r.TFile(outfile_name,"recreate")
outtree=r.TTree('events','Events')

# Branches for the discriminating variables of the ecal detector.
outtree.Branch("e2max_l00", e2max_l00, "e2max_l00/D")
outtree.Branch("emax_l00", emax_l00, "emax_l00/D")
outtree.Branch("edmax_l00", edmax_l00, "edmax_l00/D")
outtree.Branch("eocore_l00", eocore_l00, "eocore_l00/D")
outtree.Branch("w3st_l00", w3st_l00, "w3st_l00/D")
outtree.Branch("w21st_l00", w21st_l00, "w21st_l00/D")

outtree.Branch("e2max_l01", e2max_l01, "e2max_l01/D")
outtree.Branch("emax_l01", emax_l01, "emax_l01/D")
outtree.Branch("edmax_l01", edmax_l01, "edmax_l01/D")
outtree.Branch("eocore_l01", eocore_l01, "eocore_l01/D")
outtree.Branch("w3st_l01", w3st_l01, "w3st_l01/D")
outtree.Branch("w21st_l01", w21st_l01, "w21st_l01/D")

outtree.Branch("e_l10", e_l10, "e_l10/D")
outtree.Branch("e_l20", e_l20, "e_l20/D")
outtree.Branch("e_l30", e_l30, "e_l30/D")
outtree.Branch("e_l40", e_l40, "e_l40/D")
outtree.Branch("e_l50", e_l50, "e_l50/D")
outtree.Branch("e_l60", e_l60, "e_l60/D")
outtree.Branch("e_l70", e_l70, "e_l70/D")

outtree.Branch("e_l0T", e_l0T, "e_l0T/D")
outtree.Branch("e_l1T", e_l1T, "e_l1T/D")
outtree.Branch("e_l2T", e_l2T, "e_l2T/D")
outtree.Branch("e_l3T", e_l3T, "e_l3T/D")
outtree.Branch("e_l4T", e_l4T, "e_l4T/D")
outtree.Branch("e_l5T", e_l5T, "e_l5T/D")
outtree.Branch("e_l6T", e_l6T, "e_l6T/D")
outtree.Branch("e_l7T", e_l7T, "e_l7T/D")

#firstcount = 0
# Loop over all events in the input file.
numEvent = 0
for event in intree:
    ev_num[0] = numEvent
    Layers = []
    numHits = 0
    Etot = 0.
    Elayer0 = 0.
    Elayer1 = 0.
    Elayer2 = 0.
    Elayer3 = 0.
    Elayer4 = 0.
    Elayer5 = 0.
    Elayer6 = 0.
    Elayer7 = 0.

    calE = 0
    cal1E = 0

    Ecal0_Phi = []
    Ecal0_cell = []
    Ecal0_E = []
    Ecal0_Eta = []
    cal0Emax = -1.
    cal0E2max = -1.
    cal0Etamax = -1.
    cal0Eta2max = -1.
    cal0Phimax = -1.
    cal0Phi2max = -1.

    Ecal1_Phi = []
    Ecal1_cell = []
    Ecal1_E = []
    Ecal1_Eta = []
    cal1Emax = -1.
    cal1E2max = -1.
    cal1Etamax = -1.
    cal1Eta2max = -1.
    cal1Phimax = -1.
    cal1Phi2max = -1.

    # Loop over everything recorded in the ecal barrel cells for each event.
    for c in event.ECalBarrelCells:
        Layer = ecalBarrel_decoder.get(c.core.cellId, "layer")
        Layers.append(Layer)
        calE = c.core.energy
        Etot += calE
        Eta = ecalBarrel_decoder.get(c.core.cellId, "eta")

        Phi = ecalBarrel_decoder.get(c.core.cellId, "phi")

        # Choose only the first layer of the barrel.
        if Layer == 0:
            Ecal0_Eta.append(Eta)
            Ecal0_cell.append([Phi,Eta])
            Ecal0_E.append(calE)
            Ecal0_Phi.append(Phi)
            Elayer0 += calE
            # If the energy found is larger than the 1st or 2nd maximum defined in previous iteration
            # redefine them accordingly
            if calE >= cal0Emax or calE > cal0E2max:

                # If it is the first iteration just set the values as the 1st maximum
                #if len(Ecal1_E) < 2:
                if calE > cal0Emax:
                    cal0E2max = cal0Emax
                    cal0Eta2max = cal0Etamax
                    cal0Phi2max = cal0Phimax
                    cal0Emax = calE
                    cal0Etamax = Eta
                    cal0Phimax = Phi

                elif calE > cal0E2max:
                    cal0E2max = calE
                    cal0Eta2max = Eta
                    cal0Phi2max = Phi

        # Choose only the second layer of the barrel.
        elif Layer == 1:

            # Record all energies and cells with their phi and eta values in the layer.
            Ecal1_Eta.append(Eta)
            Ecal1_cell.append([Phi,Eta])
            Ecal1_E.append(calE)
            Ecal1_Phi.append(Phi)
            Elayer1 += calE
            # If the energy found is larger than the 1st or 2nd maximum defined in previous iteration
            # redefine them accordingly
            if calE >= cal1Emax or calE > cal1E2max:

                # If it is the first iteration just set the values as the 1st maximum
                #if len(Ecal1_E) < 2:
                if calE > cal1Emax:
                    cal1E2max = cal1Emax
                    cal1Eta2max = cal1Etamax
                    cal1Phi2max = cal1Phimax
                    cal1Emax = calE
                    cal1Etamax = Eta
                    cal1Phimax = Phi
                elif calE > cal1E2max:
                    cal1E2max = calE
                    cal1Eta2max = Eta
                    cal1Phi2max = Phi

        elif Layer == 2:
            Elayer2 += calE

        elif Layer == 3:
            Elayer3 += calE

        elif Layer == 4:
            Elayer4 += calE

        elif Layer == 5:
            Elayer5 += calE

        elif Layer == 6:
            Elayer6 += calE

        elif Layer == 7:
            Elayer7 += calE

        else:
            continue
        numHits += 1

    # Convert the lists into numpy arrays.
    Cellids0 = n.array(Ecal0_cell)
    Energies0 = n.array(Ecal0_E)
    Etas0 = n.array(Ecal0_Eta)
    Phis0 = n.array(Ecal0_Phi)

    Cellids1 = n.array(Ecal1_cell)
    Energies1 = n.array(Ecal1_E)
    Etas1 = n.array(Ecal1_Eta)
    Phis1 = n.array(Ecal1_Phi)

    Es10 = []
    Es00 = []
    Es11 = []
    Es01 = []
    Es1_1 = []
    Es0_1 = []

    for e in Energies0:
        if Etas0[n.where(Energies0 == e)[0][0]] == cal0Etamax:
            Es00.append(e)

        if Etas0[n.where(Energies0 == e)[0][0]] == cal0Etamax+1:
            Es01.append(e)

        if Etas0[n.where(Energies0 == e)[0][0]] == cal0Etamax-1:
            Es0_1.append(e)

        else:
            continue

    for e in Energies1:
        if Etas1[n.where(Energies1 == e)[0][0]] == cal1Etamax:
            Es10.append(e)

        if Etas1[n.where(Energies1 == e)[0][0]] == cal1Etamax+1:
            Es11.append(e)

        if Etas1[n.where(Energies1 == e)[0][0]] == cal1Etamax-1:
            Es1_1.append(e)

        else:
            continue

    E00 = n.array(Es00)
    E10 = n.array(Es10)
    E01 = n.array(Es01)
    E11 = n.array(Es11)
    E0_1 = n.array(Es0_1)
    E1_1 = n.array(Es1_1)
    #print n.sort(NewEnergies)

    # If there were more than 2 Energies recorded in the barrel, then find the extremum points of
    # these energies and sort the list of energies at these extremum indices.
    try:
        Max00 = n.sort(E00[argrelmax(E00)])
        Max01 = n.sort(E01[argrelmax(E01)])
        Max0_1 = n.sort(E0_1[argrelmax(E0_1)])
        try:
            if Phis0[n.where(Energies0 == Max00[-1])[0][0]] == cal0Phimax:
                if Phis0[n.where(Energies0 == Max00[-2])[0][0]] != cal0Phimax:
                    fmax0 = Max00[-2]

            else:
                fmax0 = Max00[-1]
        except(IndexError):
            fmax0 = 0

        try:
            if Phis0[n.where(Energies0 == Max01[-1])[0][0]] == cal0Phimax:
                if Phis0[n.where(Energies0 == Max01[-2])[0][0]] != cal0Phimax:
                    fmax1 = Max01[-2]

            else:
                fmax1 = Max01[-1]
        except(IndexError):
            fmax1 = 0

        try:
            if Phis0[n.where(Energies0 == Max0_1[-1])[0][0]] == cal0Phimax:
                if Phis0[n.where(Energies0 == Max0_1[-2])[0][0]] != cal0Phimax:
                    fmax_1 = Max0_1[-2]

            else:
                fmax_1 = Max0_1[-1]
        except(IndexError):
            fmax_1 = 0

        if (cal0E2max > fmax0) and (cal0E2max > fmax1) and (cal0E2max > fmax_1) and (m.fabs(Phis0[n.where(Energies0 == cal0E2max)[0][0]] - cal0Phimax) > 1):
            if m.fabs(cal0Phimax - Phis0[n.where(Energies0 == cal0E2max)[0][0]]) < 2:
                #print "might need finer granularity", "\n"
                cal0E2max = cal0E2max

        elif (fmax_1 > fmax0) or (fmax1 > fmax0):
            if fmax_1 < fmax1:
                cal0E2max = fmax1
                cal0Phi2max = Phis0[n.where(Energies0 == fmax1)[0][0]]
                cal0Eta2max = Etas0[n.where(Energies0 == fmax1)[0][0]]

            else:
                cal0E2max = fmax_1
                cal0Phi2max = Phis0[n.where(Energies0 == fmax_1)[0][0]]
                cal0Eta2max = Etas0[n.where(Energies0 == fmax_1)[0][0]]

        else:
            cal0E2max = fmax0
            cal0Phi2max = Phis0[n.where(Energies0 == fmax0)[0][0]]
            cal0Eta2max = cal0Etamax

    except(IndexError):
        cal0E2max = 0.0
        #print "No second maxima in first layer"
        #firstcount += 1

    try:

        Max10 = n.sort(E10[argrelmax(E10)])
        Max11 = n.sort(E11[argrelmax(E11)])
        Max1_1 = n.sort(E1_1[argrelmax(E1_1)])

        try:
            if Phis1[n.where(Energies1 == Max10[-1])[0][0]] == cal1Phimax:
                if Phis1[n.where(Energies1 == Max10[-2])[0][0]] != cal1Phimax:
                    smax0 = Max10[-2]

            else:
                smax0 = Max10[-1]
        except(IndexError):
            smax0 = 0

        try:
            if Phis1[n.where(Energies1 == Max11[-1])[0][0]] == cal1Phimax:
                if Phis1[n.where(Energies1 == Max11[-2])[0][0]] != cal1Phimax:
                    smax1 = Max11[-2]

            else:
                smax1 = Max11[-1]
        except(IndexError):
            smax1 = 0

        try:
            if Phis1[n.where(Energies1 == Max1_1[-1])[0][0]] == cal1Phimax:
                if Phis1[n.where(Energies1 == Max1_1[-2])[0][0]] != cal1Phimax:
                    smax_1 = Max1_1[-2]

            else:
                smax_1 = Max1_1[-1]
        except(IndexError):
            smax_1 = 0

        if (cal1E2max > smax0) and (cal1E2max > smax1) and (cal1E2max > smax_1) and (m.fabs(Phis1[n.where(Energies1 == cal1E2max)[0][0]] - cal1Phimax) > 1):
            if m.fabs(cal1Phimax - Phis1[n.where(Energies1 == cal1E2max)[0][0]]) < 2:
                cal1E2max = cal1E2max
                #print "might need finer granularity", "\n"

        elif (smax_1 > smax0) or (smax1 > smax0):
            if smax_1 < smax1:
                cal1E2max = smax1
                cal1Phi2max = Phis1[n.where(Energies1 == smax1)[0][0]]
                cal1Eta2max = Etas1[n.where(Energies1 == smax1)[0][0]]

            else:
                cal1E2max = smax_1
                cal1Phi2max = Phis1[n.where(Energies1 == smax_1)[0][0]]
                cal1Eta2max = Etas1[n.where(Energies1 == smax_1)[0][0]]

        else:
            cal1E2max = smax0
            cal1Phi2max = Phis1[n.where(Energies1 == smax0)[0][0]]
            cal1Eta2max = cal1Etamax

    except(IndexError):
        cal1E2max = -1
        #print "No 2nd maximum found in second layer"

    if cal0E2max > 0.0:
        w3st_l00[0] = Shower_width(Energies0, Cellids0, 1)
        w21st_l00[0] = Shower_width(Energies0, Cellids0, 21)
        eocore_l00[0] = eocorey(cal0Emax, Cellids0, Energies0, 3)
        e2max_l00[0] = cal0E2max
        emax_l00[0] = cal0Emax
        edmax_l00[0] = edmaxy(cal0Emax, cal0E2max, Cellids0, Energies0)
    else:

        w3st_l00[0] = Shower_width(Energies0, Cellids0, 1)
        w21st_l00[0] = Shower_width(Energies0, Cellids0, 21)
        eocore_l00[0] = eocorey(cal0Emax, Cellids0, Energies0, 3)
        e2max_l00[0] = 0
        emax_l00[0] = cal0Emax
        edmax_l00[0] = 0

    if cal1E2max > 0.0:
        w3st_l01[0] = Shower_width(Energies1, Cellids1, 1)
        w21st_l01[0] = Shower_width(Energies1, Cellids1, 21)
        eocore_l01[0] = eocorey(cal1Emax, Cellids1, Energies1, 3)
        e2max_l01[0] = cal1E2max
        emax_l01[0] = cal1Emax
        edmax_l01[0] = edmaxy(cal1Emax, cal1E2max, Cellids1, Energies1)
    else:
        w3st_l01[0] = Shower_width(Energies1, Cellids1, 1)
        w21st_l01[0] = Shower_width(Energies1, Cellids1, 21)
        eocore_l01[0] = eocorey(cal1Emax, Cellids1, Energies1, 3)
        e2max_l01[0] = 0
        emax_l01[0] = cal1Emax
        edmax_l01[0] = 0

    if Elayer0 > 0.0:

        e_l10[0] = Elayer1/Elayer0
        e_l20[0] = Elayer2/Elayer0
        e_l30[0] = Elayer3/Elayer0
        e_l40[0] = Elayer4/Elayer0
        e_l50[0] = Elayer5/Elayer0
        e_l60[0] = Elayer6/Elayer0
        e_l70[0] = Elayer7/Elayer0
    else:
        e_l10[0] = -1
        e_l20[0] = -1
        e_l30[0] = -1
        e_l40[0] = -1
        e_l50[0] = -1
        e_l60[0] = -1
        e_l70[0] = -1

    e_l0T[0] = Elayer0/Etot
    e_l1T[0] = Elayer1/Etot
    e_l2T[0] = Elayer2/Etot
    e_l3T[0] = Elayer3/Etot
    e_l4T[0] = Elayer4/Etot
    e_l5T[0] = Elayer5/Etot
    e_l6T[0] = Elayer6/Etot
    e_l7T[0] = Elayer7/Etot

    outtree.Fill()

    numEvent += 1

outtree.Write()
outfile.Write()
outfile.Close()
