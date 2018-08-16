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
ecalBarrel_decoder = dd4hep.DDSegmentation.BitFieldCoder("system:4,cryo:1,type:3,subtype:3,layer:8,eta:11,phi:11")
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
largestPhiValue = int(1408)

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

def Find2ndMax(Energies,Cells,Distance):
    try:
        Emax = n.amax(Energies)
        E2max = 0
        Etamax = Cells[n.argmax(Energies)][1]
        Phimax = Cells[n.argmax(Energies)][0]
    except(ValueError):
        return 0.

    if len(Energies) % 2 == 0:
        window = len(Energies) - 1
    else:
        window = len(Energies)
    try:
        Filtered_Energies = savgol_filter(Energies, window, 5, mode='nearest')

    except(ValueError):
        try:
            Filtered_Energies = savgol_filter(Energies, window, 3, mode='nearest')

        except(ValueError):
            Filtered_Energies = n.array([])
    try:
        Extrema = n.sort(Energies[argrelextrema(Filtered_Energies,n.greater)])

        if Extrema == n.array([]):
            E2max = 0.0

        else:
            newEta = Cells[n.where(Energies == Extrema[-1])[0][0]][1]
            newPhi = Cells[n.where(Energies == Extrema[-1])[0][0]][0]
            de = int(m.fabs(newEta-Etamax))
            dp = int(m.fabs(newPhi-Phimax))

            if dp > (largestPhiValue/2):
                dp = int(m.fabs(dp - (largestPhiValue+1)))

            dist = m.sqrt(de**2 + dp**2)

            if Extrema[-1] != Emax and dist > Distance:
                E2max = Extrema[-1]
                Phi2max = newPhi
                Eta2max = newEta

            else:
                E2max = Extrema[-2]
                Phi2max = Cells[n.where(Energies == Extrema[-2])[0][0]][0]
                Eta2max = Cells[n.where(Energies == Extrema[-2])[0][0]][1]

    except(IndexError):
        E2max = 0.0

    for e in Energies:
        newEta = Cells[n.where(Energies == e)[0][0]][1]
        newPhi = Cells[n.where(Energies == e)[0][0]][0]
        de = int(m.fabs(newEta-Etamax))
        dp = int(m.fabs(newPhi-Phimax))
        if dp > (largestPhiValue/2):
            dp = int(m.fabs(dp - (largestPhiValue+1)))

        dist = m.sqrt(de**2 + dp**2)

        if e > E2max and e > Emax*0.6 and dist > Distance:
            E2max = e
            Phi2max = newEta
            Eta2max = newPhi


    return E2max

# Function for calculating the shower width on j cells in cell units. imax describes the cells
# cell with maximum energy deposit.

def Shower_width(Energies, Cells, j):
    #topphi = 0
    #botphi = 0
    topeta = 0
    boteta = 0
    #wphi = int(j/1.926)
    #limphi = int((wphi-1)/2)
    limeta = int((j-1)/2)

    try:
        imaxphi = int(Cells[n.argmax(Energies)][0])
        imaxeta = int(Cells[n.argmax(Energies)][1])

        # Loop over all detected cells.
        for cell in Cells:
            phi = int(cell[0])
            eta = int(cell[1])
            dp = int(m.fabs(phi - imaxphi))
            de = int(m.fabs(eta - imaxeta))

            try:
                if (dp <= 1 or dp >= (largestPhiValue-1)) and de <= limeta:
                    topeta += sum(Energies[n.where(Cells == cell)[0]])*(((eta-imaxeta)**2))
                    boteta += sum(Energies[n.where(Cells == cell)[0]])

                else:
                    topeta += 0
                    boteta += 0

            except(IndexError):
                topeta += 0
                boteta += 0

        Weta = (topeta/boteta)
        return Weta
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
        dmaxphi = int(m.fabs(cellmaxphi - cell2maxphi))
        dmaxeta = int(m.fabs(cellmaxeta - cell2maxeta))

        if (dmaxphi <= 1 or dmaxphi == largestPhiValue) and dmaxeta <= 1:
            return 0.

        else:
            # Loop over all the cells recorded.
            while i < len(Cells):
                phi = int(Cells[i][0])
                eta = int(Cells[i][1])
                dp = int(m.fabs(int(Cells[i][0]) - cellmaxphi))
                de = int(m.fabs(int(Cells[i][1]) - cellmaxeta))

                # Condition to choose only cells that are on the same phi and eta values as the 1st and
                # 2nd maximum or inbetween these phi and eta values. Record all the energies in these cells.
                if (phi >= cell2maxphi and phi <= cellmaxphi) or (phi >= cellmaxphi and phi <= cell2maxphi):
                    if (eta >= cell2maxeta and eta <= cellmaxeta) or (eta >= cellmaxeta and eta <= cell2maxeta):
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

            edmax = E2ndmax - minE
            return edmax

    except(IndexError):
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
            de = int(m.fabs(int(Cells[i][1]) - cellmaxeta))
            dp = int(m.fabs(int(Cells[i][0]) - cellmaxphi))
            if dp >= (largestPhiValue+1) - j:
                dp = int(m.fabs(dp-(largestPhiValue+1)))

            # Choose only cells that are within +/- j cells around (and including) the cell with maximum energy deposit
            # and sum the energies deposited in these cells.
            if m.sqrt(de**2 + dp**2) <= j:
                E3 += Energies[i]

                # Choose only cells that are within +/- 1 cells around (and including) the cell with maximum energy deposit
                # and sum the energies deposited in these cells.
                if m.sqrt(de**2 + dp**2) <= 1:
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
w9st_l01 = n.zeros(1, dtype=float)
w41st_l01 = n.zeros(1, dtype=float)

# Variables for third layer.
e2max_l02 = n.zeros(1, dtype=float)
emax_l02 = n.zeros(1, dtype=float)
edmax_l02 = n.zeros(1, dtype=float)
eocore_l02 = n.zeros(1, dtype=float)
w3st_l02 = n.zeros(1, dtype=float)
w21st_l02 = n.zeros(1, dtype=float)

# Variables for fourth layer.
e2max_l03 = n.zeros(1, dtype=float)
emax_l03 = n.zeros(1, dtype=float)
edmax_l03 = n.zeros(1, dtype=float)
eocore_l03 = n.zeros(1, dtype=float)
w3st_l03 = n.zeros(1, dtype=float)
w21st_l03 = n.zeros(1, dtype=float)

# Variables for fifth layer.
e2max_l04 = n.zeros(1, dtype=float)
emax_l04 = n.zeros(1, dtype=float)
edmax_l04 = n.zeros(1, dtype=float)
eocore_l04 = n.zeros(1, dtype=float)
w3st_l04 = n.zeros(1, dtype=float)
w21st_l04 = n.zeros(1, dtype=float)

# Variables for layer energies against first layer energy and total energy
e_l01 = n.zeros(1, dtype=float)
e_l21 = n.zeros(1, dtype=float)
e_l31 = n.zeros(1, dtype=float)
e_l41 = n.zeros(1, dtype=float)
e_l51 = n.zeros(1, dtype=float)
e_l61 = n.zeros(1, dtype=float)
e_l71 = n.zeros(1, dtype=float)

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
outtree.Branch("w9st_l01", w9st_l01, "w9st_l01/D")
outtree.Branch("w41st_l01", w41st_l01, "w41st_l01/D")

outtree.Branch("e2max_l02", e2max_l02, "e2max_l02/D")
outtree.Branch("emax_l02", emax_l02, "emax_l02/D")
outtree.Branch("edmax_l02", edmax_l02, "edmax_l02/D")
outtree.Branch("eocore_l02", eocore_l02, "eocore_l02/D")
outtree.Branch("w3st_l02", w3st_l02, "w3st_l02/D")
outtree.Branch("w21st_l02", w21st_l02, "w21st_l02/D")

outtree.Branch("e2max_l03", e2max_l03, "e2max_l03/D")
outtree.Branch("emax_l03", emax_l03, "emax_l03/D")
outtree.Branch("edmax_l03", edmax_l03, "edmax_l03/D")
outtree.Branch("eocore_l03", eocore_l03, "eocore_l03/D")
outtree.Branch("w3st_l03", w3st_l03, "w3st_l03/D")
outtree.Branch("w21st_l03", w21st_l03, "w21st_l03/D")

outtree.Branch("e2max_l04", e2max_l04, "e2max_l04/D")
outtree.Branch("emax_l04", emax_l04, "emax_l04/D")
outtree.Branch("edmax_l04", edmax_l04, "edmax_l04/D")
outtree.Branch("eocore_l04", eocore_l04, "eocore_l04/D")
outtree.Branch("w3st_l04", w3st_l04, "w3st_l04/D")
outtree.Branch("w21st_l04", w21st_l04, "w21st_l04/D")

outtree.Branch("e_l01", e_l01, "e_l01/D")
outtree.Branch("e_l21", e_l21, "e_l21/D")
outtree.Branch("e_l31", e_l31, "e_l31/D")
outtree.Branch("e_l41", e_l41, "e_l41/D")
outtree.Branch("e_l51", e_l51, "e_l51/D")
outtree.Branch("e_l61", e_l61, "e_l61/D")
outtree.Branch("e_l71", e_l71, "e_l71/D")

outtree.Branch("e_l0T", e_l0T, "e_l0T/D")
outtree.Branch("e_l1T", e_l1T, "e_l1T/D")
outtree.Branch("e_l2T", e_l2T, "e_l2T/D")
outtree.Branch("e_l3T", e_l3T, "e_l3T/D")
outtree.Branch("e_l4T", e_l4T, "e_l4T/D")
outtree.Branch("e_l5T", e_l5T, "e_l5T/D")
outtree.Branch("e_l6T", e_l6T, "e_l6T/D")
outtree.Branch("e_l7T", e_l7T, "e_l7T/D")

# Loop over all events in the input file.
numEvent = 0

for event in intree:
    ev_num[0] = numEvent
    Layers = []
    numHits = 0
    calE = 0
    cal1E = 0

    Etot = 0.
    Elayer0 = 0.
    Elayer1 = 0.
    Elayer2 = 0.
    Elayer3 = 0.
    Elayer4 = 0.
    Elayer5 = 0.
    Elayer6 = 0.
    Elayer7 = 0.

    Ecal0_cell = []
    Ecal0_E = []
    cal0Emax = 0.

    Ecal1_cell = []
    Ecal1_E = []
    cal1Emax = 0.

    Ecal2_cell = []
    Ecal2_E = []
    cal2Emax = 0.

    Ecal3_cell = []
    Ecal3_E = []
    cal3Emax = 0.

    Ecal4_cell = []
    Ecal4_E = []
    cal4Emax = 0.

    # Loop over everything recorded in the ecal barrel cells for each event.
    for c in event.ECalBarrelCells:
        Layer = ecalBarrel_decoder.get(c.core.cellId, "layer")
        Layers.append(Layer)
        calE = c.core.energy
        Etot += calE
        Eta = ecalBarrel_decoder.get(c.core.cellId, "eta")
        Phi = ecalBarrel_decoder.get(c.core.cellId, "phi")
        if Phi > 1408:
            print "Error, Phi got out of hand (greater than 704)"
        if Eta > 1356:
            print "Error, Eta larger than 1356"

        # Choose only the first layer of the barrel.
        if Layer == 0:
            Ecal0_cell.append([Phi,Eta])
            Ecal0_E.append(calE)
            Elayer0 += calE

            # If the energy found is larger than the 1st maximum defined in previous iteration
            # redefine it.
            if calE > cal0Emax:
                cal0Emax = calE

        # Choose only the second layer of the barrel.
        elif Layer == 1:
            Ecal1_cell.append([Phi,Eta])
            Ecal1_E.append(calE)
            Elayer1 += calE

            if calE > cal1Emax:
                cal1Emax = calE

        elif Layer == 2:
            Elayer2 += calE
            Ecal2_cell.append([Phi,Eta])
            Ecal2_E.append(calE)

            if calE > cal2Emax:
                cal2Emax = calE

        elif Layer == 3:
            Elayer3 += calE
            Ecal3_cell.append([Phi,Eta])
            Ecal3_E.append(calE)

            if calE > cal3Emax:
                cal3Emax = calE

        elif Layer == 4:
            Elayer4 += calE
            Ecal4_cell.append([Phi,Eta])
            Ecal4_E.append(calE)

            if calE > cal4Emax:
                cal4Emax = calE

        elif Layer == 5:
            Elayer5 += calE

        elif Layer == 6:
            Elayer6 += calE

        elif Layer == 7:
            Elayer7 += calE
        elif Layer > 7:
            print "Error, too large layers found."

        else:
            continue
        numHits += 1

    # convert the lists into numpy arrays.
    Cellids0 = n.array(Ecal0_cell)
    Energies0 = n.array(Ecal0_E)

    Cellids1 = n.array(Ecal1_cell)
    Energies1 = n.array(Ecal1_E)

    Cellids2 = n.array(Ecal2_cell)
    Energies2 = n.array(Ecal2_E)

    Cellids3 = n.array(Ecal3_cell)
    Energies3 = n.array(Ecal3_E)

    Cellids4 = n.array(Ecal4_cell)
    Energies4 = n.array(Ecal4_E)
    
    emax_l00[0] = cal0Emax
    e2max_l00[0] = Find2ndMax(Energies0,Cellids0,6)
    w3st_l00[0] = Shower_width(Energies0, Cellids0, 3)
    w21st_l00[0] = Shower_width(Energies0, Cellids0, 21)
    eocore_l00[0] = eocorey(cal0Emax, Cellids0, Energies0, 3)
    edmax_l00[0] = edmaxy(cal0Emax, e2max_l00[0], Cellids0, Energies0)

    emax_l01[0] = cal1Emax
    e2max_l01[0] = Find2ndMax(Energies1,Cellids1,12)
    w9st_l01[0] = Shower_width(Energies1, Cellids1, 9)
    w41st_l01[0] = Shower_width(Energies1, Cellids1, 41)
    eocore_l01[0] = eocorey(cal1Emax, Cellids1, Energies1, 5)
    edmax_l01[0] = edmaxy(cal1Emax, e2max_l01[0], Cellids1, Energies1)

    emax_l02[0] = cal2Emax
    e2max_l02[0] = Find2ndMax(Energies2,Cellids2,6)
    w3st_l02[0] = Shower_width(Energies2, Cellids2, 3)
    w21st_l02[0] = Shower_width(Energies2, Cellids2, 21)
    eocore_l02[0] = eocorey(cal2Emax, Cellids2, Energies2, 3)
    edmax_l02[0] = edmaxy(cal2Emax, e2max_l02[0], Cellids2, Energies2)

    emax_l03[0] = cal3Emax
    e2max_l03[0] = Find2ndMax(Energies3,Cellids3,6)
    w3st_l03[0] = Shower_width(Energies3, Cellids3, 3)
    w21st_l03[0] = Shower_width(Energies3, Cellids3, 21)
    eocore_l03[0] = eocorey(cal3Emax, Cellids3, Energies3, 3)
    edmax_l03[0] = edmaxy(cal3Emax, e2max_l03[0], Cellids3, Energies3)

    emax_l04[0] = cal4Emax
    e2max_l04[0] = Find2ndMax(Energies4,Cellids4,6)
    w3st_l04[0] = Shower_width(Energies4, Cellids4, 3)
    w21st_l04[0] = Shower_width(Energies4, Cellids4, 21)
    eocore_l04[0] = eocorey(cal4Emax, Cellids4, Energies4, 3)
    edmax_l04[0] = edmaxy(cal4Emax, e2max_l04[0], Cellids4, Energies4)

    if Elayer1 > 0.1:
        e_l01[0] = Elayer0/Elayer1
        e_l21[0] = Elayer2/Elayer1
        e_l31[0] = Elayer3/Elayer1
        e_l41[0] = Elayer4/Elayer1
        e_l51[0] = Elayer5/Elayer1
        e_l61[0] = Elayer6/Elayer1
        e_l71[0] = Elayer7/Elayer1


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
