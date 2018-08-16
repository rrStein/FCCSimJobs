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
ecalBarrel_decoder = dd4hep.DDSegmentation.BitFieldCoder("system:4,cryo:1,type:3,subtype:3,layer:8,eta:11,phi:10")
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
largestPhiValue = int(704)

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
    topphi = 0
    botphi = 0
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

    Ecal0_Phi = []
    Ecal0_cell = []
    Ecal0_E = []
    Ecal0_Eta = []
    cal0Emax = 0.
    cal0E2max = 0.
    cal0Etamax = 0.
    cal0Eta2max = 0.
    cal0Phimax = 0.
    cal0Phi2max = 0.

    Ecal1_Phi = []
    Ecal1_cell = []
    Ecal1_E = []
    Ecal1_Eta = []
    cal1Emax = 0.
    cal1E2max = 0.
    cal1Etamax = 0.
    cal1Eta2max = 0.
    cal1Phimax = 0.
    cal1Phi2max = 0.

    Ecal2_Phi = []
    Ecal2_cell = []
    Ecal2_E = []
    Ecal2_Eta = []
    cal2Emax = 0.
    cal2E2max = 0.
    cal2Etamax = 0.
    cal2Eta2max = 0.
    cal2Phimax = 0.
    cal2Phi2max = 0.

    Ecal3_Phi = []
    Ecal3_cell = []
    Ecal3_E = []
    Ecal3_Eta = []
    cal3Emax = 0.
    cal3E2max = 0.
    cal3Etamax = 0.
    cal3Eta2max = 0.
    cal3Phimax = 0.
    cal3Phi2max = 0.

    Ecal4_Phi = []
    Ecal4_cell = []
    Ecal4_E = []
    Ecal4_Eta = []
    cal4Emax = 0.
    cal4E2max = 0.
    cal4Etamax = 0.
    cal4Eta2max = 0.
    cal4Phimax = 0.
    cal4Phi2max = 0.

    # Loop over everything recorded in the ecal barrel cells for each event.
    for c in event.ECalBarrelCells:
        Layer = ecalBarrel_decoder.get(c.core.cellId, "layer")
        Layers.append(Layer)
        calE = c.core.energy
        Etot += calE
        Eta = ecalBarrel_decoder.get(c.core.cellId, "eta")
        Phi = ecalBarrel_decoder.get(c.core.cellId, "phi")
        if Phi > 704:
            print "Error, Phi got out of hand (greater than 704)"
        if Eta > 1356:
            print "Error, Eta larger than 1356"

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
            # Record all energies and cells with their phi and eta values in the layer.
            Ecal2_Eta.append(Eta)
            Ecal2_cell.append([Phi,Eta])
            Ecal2_E.append(calE)
            Ecal2_Phi.append(Phi)
            # If the energy found is larger than the 1st or 2nd maximum defined in previous iteration
            # redefine them accordingly
            if calE >= cal2Emax or calE > cal2E2max:
                if calE > cal2Emax:
                    cal2E2max = cal2Emax
                    cal2Eta2max = cal2Etamax
                    cal2Phi2max = cal2Phimax
                    cal2Emax = calE
                    cal2Etamax = Eta
                    cal2Phimax = Phi

                elif calE > cal2E2max:
                    cal2E2max = calE
                    cal2Eta2max = Eta
                    cal2Phi2max = Phi

        elif Layer == 3:
            Elayer3 += calE
            # Record all energies and cells with their phi and eta values in the layer.
            Ecal3_Eta.append(Eta)
            Ecal3_cell.append([Phi,Eta])
            Ecal3_E.append(calE)
            Ecal3_Phi.append(Phi)

            # If the energy found is larger than the 1st or 2nd maximum defined in previous iteration
            # redefine them accordingly
            if calE >= cal3Emax or calE > cal3E2max:
                if calE > cal3Emax:
                    cal3E2max = cal3Emax
                    cal3Eta2max = cal3Etamax
                    cal3Phi2max = cal3Phimax
                    cal3Emax = calE
                    cal3Etamax = Eta
                    cal3Phimax = Phi

                elif calE > cal3E2max:
                    cal3E2max = calE
                    cal3Eta2max = Eta
                    cal3Phi2max = Phi

        elif Layer == 4:
            Elayer4 += calE
            # Record all energies and cells with their phi and eta values in the layer.
            Ecal4_Eta.append(Eta)
            Ecal4_cell.append([Phi,Eta])
            Ecal4_E.append(calE)
            Ecal4_Phi.append(Phi)

            # If the energy found is larger than the 1st or 2nd maximum defined in previous iteration
            # redefine them accordingly
            if calE >= cal4Emax or calE > cal4E2max:
                if calE > cal4Emax:
                    cal4E2max = cal4Emax
                    cal4Eta2max = cal4Etamax
                    cal4Phi2max = cal4Phimax
                    cal4Emax = calE
                    cal4Etamax = Eta
                    cal4Phimax = Phi

                elif calE > cal4E2max:
                    cal4E2max = calE
                    cal4Eta2max = Eta
                    cal4Phi2max = Phi

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
    Etas0 = n.array(Ecal0_Eta)
    Phis0 = n.array(Ecal0_Phi)

    Cellids1 = n.array(Ecal1_cell)
    Energies1 = n.array(Ecal1_E)
    Etas1 = n.array(Ecal1_Eta)
    Phis1 = n.array(Ecal1_Phi)

    Cellids2 = n.array(Ecal2_cell)
    Energies2 = n.array(Ecal2_E)
    Etas2 = n.array(Ecal2_Eta)
    Phis2 = n.array(Ecal2_Phi)

    Cellids3 = n.array(Ecal3_cell)
    Energies3 = n.array(Ecal3_E)
    Etas3 = n.array(Ecal3_Eta)
    Phis3 = n.array(Ecal3_Phi)

    Cellids4 = n.array(Ecal4_cell)
    Energies4 = n.array(Ecal4_E)
    Etas4 = n.array(Ecal4_Eta)
    Phis4 = n.array(Ecal4_Phi)

    '''print "Phi: ", "\n"
    print Phis0, "\n"
    print "Eta: ", "\n"
    print Etas0
    print cal0Etamax, "\n"'''

    Es00 = []
    Es10 = []
    Es20 = []
    Es30 = []
    Es40 = []

    if len(Energies0) % 2 == 0:
        window0 = len(Energies0) - 1

    else:
        window0 = len(Energies0)

    if len(Energies1) % 2 == 0:
        window1 = len(Energies1) - 1

    else:
        window1 = len(Energies1)

    if len(Energies2) % 2 == 0:
        window2 = len(Energies2) - 1

    else:
        window2 = len(Energies2)

    if len(Energies3) % 2 == 0:
        window3 = len(Energies3) - 1

    else:
        window3 = len(Energies3)

    if len(Energies4) % 2 == 0:
        window4 = len(Energies4) - 1

    else:
        window4 = len(Energies4)

    try:
        E00 = savgol_filter(Energies0, window0, 5, mode='nearest')

    except(ValueError):
        try:
            E00 = savgol_filter(Energies0, window0, 3, mode='nearest')

        except(ValueError):
            E00 = n.array([])

    try:
        E10 = savgol_filter(Energies1, window1, 5, mode='nearest')

    except(ValueError):
        try:
            E10 = savgol_filter(Energies1, window1, 3, mode='nearest')

        except(ValueError):
            E10 = n.array([])

    try:
        E20 = savgol_filter(Energies2, window2, 5, mode='nearest')

    except(ValueError):
        try:
            E20 = savgol_filter(Energies2, window2, 3, mode='nearest')

        except(ValueError):
            E20 = n.array([])

    try:
        E30 = savgol_filter(Energies3, window3, 5, mode='nearest')

    except(ValueError):
        try:
            E30 = savgol_filter(Energies3, window3, 3, mode='nearest')

        except(ValueError):
            E30 = n.array([])

    try:
        E40 = savgol_filter(Energies4, window4, 5, mode='nearest')

    except(ValueError):
        try:
            E40 = savgol_filter(Energies4, window4, 3, mode='nearest')

        except(ValueError):
            E40 = n.array([])

    try:
        es00 = n.sort(Energies0[argrelextrema(E00,n.greater)])

        if es00 == n.array([]):
            cal0E2max = 0.0

        else:
            newEta0 = Cellids0[n.where(Energies0 == es00[-1])[0][0]][1]
            newPhi0 = Cellids0[n.where(Energies0 == es00[-1])[0][0]][0]
            de0 = int(m.fabs(newEta0-cal0Etamax))
            dp0 = int(m.fabs(newPhi0-cal0Phimax))
            if dp0 > (largestPhiValue/2):
                dp0 = int(m.fabs(dp0 - (largestPhiValue+1)))

            dist0 = m.sqrt(de0**2 + dp0**2)

            if es00[-1] != cal0Emax and dist0 > 2:
                cal0E2max = es00[-1]
                cal0Phi2max = Phis0[n.where(Energies0 == es00[-1])[0][0]]
                cal0Eta2max = Etas0[n.where(Energies0 == es00[-1])[0][0]]

            else:
                cal0E2max = es00[-2]
                cal0Phi2max = Phis0[n.where(Energies0 == es00[-2])[0][0]]
                cal0Eta2max = Etas0[n.where(Energies0 == es00[-2])[0][0]]

    except(IndexError):
        cal0E2max = 0.0

    try:
        es10 =  n.sort(Energies1[argrelextrema(E10,n.greater)])
        #print es10
        if es10 == n.array([]):
            cal1E2max = 0.0

        else:
            newEta1 = Cellids1[n.where(Energies1 == es10[-1])[0][0]][1]
            newPhi1 = Cellids1[n.where(Energies1 == es10[-1])[0][0]][0]
            de1 = int(m.fabs(newEta1-cal1Etamax))
            dp1 = int(m.fabs(newPhi1-cal1Phimax))
            if dp1 > (largestPhiValue/2):
                dp1 = int(m.fabs(dp1 - (largestPhiValue+1)))

            dist1 = m.sqrt(de1**2 + dp1**2)

            if es10[-1] != cal1Emax and dist1 > 2:
                cal1E2max = es10[-1]
                cal1Phi2max = Phis1[n.where(Energies1 == es10[-1])[0][0]]
                cal1Eta2max = Etas1[n.where(Energies1 == es10[-1])[0][0]]

            else:
                cal1E2max = es10[-2]
                cal1Phi2max = Phis1[n.where(Energies1 == es10[-2])[0][0]]
                cal1Eta2max = Etas1[n.where(Energies1 == es10[-2])[0][0]]

    except(IndexError):
        cal1E2max = 0.0

    try:
        es20 = n.sort(Energies2[argrelextrema(E20,n.greater)])

        if es20 == n.array([]):
            cal2E2max = 0.0

        else:
            newEta2 = Cellids2[n.where(Energies2 == es20[-1])[0][0]][1]
            newPhi2 = Cellids2[n.where(Energies2 == es20[-1])[0][0]][0]
            de2 = int(m.fabs(newEta2-cal2Etamax))
            dp2 = int(m.fabs(newPhi2-cal2Phimax))
            if dp2 > (largestPhiValue/2):
                dp2 = int(m.fabs(dp2 - (largestPhiValue+1)))

            dist2 = m.sqrt(de2**2 + dp2**2)

            if es20[-1] != cal2Emax and dist2 > 2:
                cal2E2max = es20[-1]
                cal2Phi2max = Phis2[n.where(Energies2 == es20[-1])[0][0]]
                cal2Eta2max = Etas2[n.where(Energies2 == es20[-1])[0][0]]

            else:
                cal2E2max = es20[-2]
                cal2Phi2max = Phis2[n.where(Energies2 == es20[-2])[0][0]]
                cal2Eta2max = Etas2[n.where(Energies2 == es20[-2])[0][0]]

    except(IndexError):
        cal2E2max = 0.0

    try:
        es30 =  n.sort(Energies3[argrelextrema(E30,n.greater)])
        #print es30
        if es30 == n.array([]):
            cal3E2max = 0.0

        else:
            newEta3 = Cellids3[n.where(Energies3 == es30[-1])[0][0]][1]
            newPhi3 = Cellids3[n.where(Energies3 == es30[-1])[0][0]][0]
            de3 = int(m.fabs(newEta3-cal3Etamax))
            dp3 = int(m.fabs(newPhi3-cal3Phimax))
            if dp3 > (largestPhiValue/2):
                dp3 = int(m.fabs(dp3 - (largestPhiValue+1)))

            dist3 = m.sqrt(de3**2 + dp3**2)

            if es30[-1] != cal3Emax and dist3 > 2:
                cal3E2max = es30[-1]
                cal3Phi2max = Phis3[n.where(Energies3 == es30[-1])[0][0]]
                cal3Eta2max = Etas3[n.where(Energies3 == es30[-1])[0][0]]

            else:
                cal3E2max = es30[-2]
                cal3Phi2max = Phis3[n.where(Energies3 == es30[-2])[0][0]]
                cal3Eta2max = Etas3[n.where(Energies3 == es30[-2])[0][0]]

    except(IndexError):
        cal3E2max = 0.0

    try:
        es40 =  n.sort(Energies4[argrelextrema(E40,n.greater)])
        #print es40
        if es40 == n.array([]):
            cal4E2max = 0.0

        else:
            newEta4 = Cellids3[n.where(Energies4 == es40[-1])[0][0]][1]
            newPhi4 = Cellids3[n.where(Energies4 == es40[-1])[0][0]][0]
            de4 = int(m.fabs(newEta4-cal4Etamax))
            dp4 = int(m.fabs(newPhi4-cal4Phimax))
            if dp4 > (largestPhiValue/2):
                dp4 = int(m.fabs(dp4 - (largestPhiValue+1)))

            dist4 = m.sqrt(de4**2 + dp4**2)

            if es40[-1] != cal4Emax and dist4 > 2:
                cal4E2max = es40[-1]
                cal4Phi2max = Phis4[n.where(Energies4 == es40[-1])[0][0]]
                cal4Eta2max = Etas4[n.where(Energies4 == es40[-1])[0][0]]

            else:
                cal4E2max = es40[-2]
                cal4Phi2max = Phis4[n.where(Energies4 == es40[-2])[0][0]]
                cal4Eta2max = Etas4[n.where(Energies4 == es40[-2])[0][0]]

    except(IndexError):
        cal4E2max = 0.0

    for e in Energies0:
        newEta0 = Cellids0[n.where(Energies0 == e)[0][0]][1]
        newPhi0 = Cellids0[n.where(Energies0 == e)[0][0]][0]
        de = int(m.fabs(newEta0-cal0Etamax))
        dp = int(m.fabs(newPhi0-cal0Phimax))
        if dp > (largestPhiValue/2):
            dp = int(m.fabs(dp - (largestPhiValue+1)))

        dist = m.sqrt(de**2 + dp**2)

        if e > cal0E2max and e > cal0Emax*0.5 and dist > 2:
            cal0E2max = e
            cal0Phi2max = newEta0
            cal0Eta2max = newPhi0

    for e in Energies1:
        newEta1 = Cellids1[n.where(Energies1 == e)[0][0]][1]
        newPhi1 = Cellids1[n.where(Energies1 == e)[0][0]][0]
        de = int(m.fabs(newEta1-cal1Etamax))
        dp = int(m.fabs(newPhi1-cal1Phimax))
        if dp > (largestPhiValue/2):
            dp = int(m.fabs(dp - (largestPhiValue+1)))

        dist = m.sqrt(de**2 + dp**2)

        if e > cal1E2max and e > cal1Emax*0.5 and dist > 4:
            cal1E2max = e
            cal1Phi2max = newEta1
            cal1Eta2max = newPhi1

    for e in Energies2:
        newEta2 = Cellids2[n.where(Energies2 == e)[0][0]][1]
        newPhi2 = Cellids2[n.where(Energies2 == e)[0][0]][0]
        de = int(m.fabs(newEta2-cal2Etamax))
        dp = int(m.fabs(newPhi2-cal2Phimax))
        if dp > (largestPhiValue/2):
            dp = int(m.fabs(dp - (largestPhiValue+1)))

        dist = m.sqrt(de**2 + dp**2)

        if e > cal2E2max and e > cal2Emax*0.5 and dist > 2:
            cal2E2max = e
            cal2Phi2max = newEta2
            cal2Eta2max = newPhi2

    for e in Energies3:
        newEta3 = Cellids3[n.where(Energies3 == e)[0][0]][1]
        newPhi3 = Cellids3[n.where(Energies3 == e)[0][0]][0]
        de = int(m.fabs(newEta3-cal3Etamax))
        dp = int(m.fabs(newPhi3-cal3Phimax))
        if dp > (largestPhiValue/2):
            dp = int(m.fabs(dp - (largestPhiValue+1)))

        dist = m.sqrt(de**2 + dp**2)

        if e > cal3E2max and e > cal3Emax*0.5 and dist > 2:
            cal3E2max = e
            cal3Phi2max = newEta3
            cal3Eta2max = newPhi3

    for e in Energies4:
        newEta4 = Cellids4[n.where(Energies4 == e)[0][0]][1]
        newPhi4 = Cellids4[n.where(Energies4 == e)[0][0]][0]
        de = int(m.fabs(newEta4-cal4Etamax))
        dp = int(m.fabs(newPhi4-cal4Phimax))
        if dp > (largestPhiValue/2):
            dp = int(m.fabs(dp - (largestPhiValue+1)))

        dist = m.sqrt(de**2 + dp**2)

        if e > cal4E2max and e > cal4Emax*0.5 and dist > 2:
            cal4E2max = e
            cal4Phi2max = newEta4
            cal4Eta2max = newPhi4

    #sorted = n.sort(Energies1)
    #print cal1E2max
    #print cal1Emax
    #print "\n", sorted,Cellids1[n.where(Energies1 == sorted[-1])[0][0]],Cellids1[n.where(Energies1 == sorted[-2])[0][0]],Cellids1[n.where(Energies1 == cal1E2max)[0][0]], "\n"

    if cal0E2max > 0:
        w3st_l00[0] = Shower_width(Energies0, Cellids0, 3)
        w21st_l00[0] = Shower_width(Energies0, Cellids0, 21)
        eocore_l00[0] = eocorey(cal0Emax, Cellids0, Energies0, 3)
        e2max_l00[0] = cal0E2max
        emax_l00[0] = cal0Emax
        edmax_l00[0] = edmaxy(cal0Emax, cal0E2max, Cellids0, Energies0)

    else:
        w3st_l00[0] = Shower_width(Energies0, Cellids0, 3)
        w21st_l00[0] = Shower_width(Energies0, Cellids0, 21)
        eocore_l00[0] = eocorey(cal0Emax, Cellids0, Energies0, 3)
        e2max_l00[0] = 0
        emax_l00[0] = cal0Emax
        edmax_l00[0] = 0

    if cal1E2max > 0:
        w9st_l01[0] = Shower_width(Energies1, Cellids1, 9)
        w41st_l01[0] = Shower_width(Energies1, Cellids1, 41)
        eocore_l01[0] = eocorey(cal1Emax, Cellids1, Energies1, 3)
        e2max_l01[0] = cal1E2max
        emax_l01[0] = cal1Emax
        edmax_l01[0] = edmaxy(cal1Emax, cal1E2max, Cellids1, Energies1)

    else:
        w9st_l01[0] = Shower_width(Energies1, Cellids1, 9)
        w41st_l01[0] = Shower_width(Energies1, Cellids1, 41)
        eocore_l01[0] = eocorey(cal1Emax, Cellids1, Energies1, 3)
        e2max_l01[0] = 0
        emax_l01[0] = cal1Emax
        edmax_l01[0] = 0

    if cal2E2max > 0:
        w3st_l02[0] = Shower_width(Energies2, Cellids2, 3)
        w21st_l02[0] = Shower_width(Energies2, Cellids2, 21)
        eocore_l02[0] = eocorey(cal2Emax, Cellids2, Energies2, 3)
        e2max_l02[0] = cal2E2max
        emax_l02[0] = cal2Emax
        edmax_l02[0] = edmaxy(cal2Emax, cal2E2max, Cellids2, Energies2)

    else:
        w3st_l02[0] = Shower_width(Energies2, Cellids2, 3)
        w21st_l02[0] = Shower_width(Energies2, Cellids2, 21)
        eocore_l02[0] = eocorey(cal2Emax, Cellids2, Energies2, 3)
        e2max_l02[0] = 0
        emax_l02[0] = cal2Emax
        edmax_l02[0] = 0

    if cal3E2max > 0:
        w3st_l03[0] = Shower_width(Energies3, Cellids3, 3)
        w21st_l03[0] = Shower_width(Energies3, Cellids3, 21)
        eocore_l03[0] = eocorey(cal3Emax, Cellids3, Energies3, 3)
        e2max_l03[0] = cal3E2max
        emax_l03[0] = cal3Emax
        edmax_l03[0] = edmaxy(cal3Emax, cal3E2max, Cellids3, Energies3)

    else:
        w3st_l03[0] = Shower_width(Energies3, Cellids3, 3)
        w21st_l03[0] = Shower_width(Energies3, Cellids3, 21)
        eocore_l03[0] = eocorey(cal3Emax, Cellids3, Energies3, 3)
        e2max_l03[0] = 0
        emax_l03[0] = cal3Emax
        edmax_l03[0] = 0

    if cal4E2max > 0:
        w3st_l04[0] = Shower_width(Energies4, Cellids4, 3)
        w21st_l04[0] = Shower_width(Energies4, Cellids4, 21)
        eocore_l04[0] = eocorey(cal4Emax, Cellids4, Energies4, 3)
        e2max_l04[0] = cal4E2max
        emax_l04[0] = cal4Emax
        edmax_l04[0] = edmaxy(cal4Emax, cal4E2max, Cellids4, Energies4)

    else:
        w3st_l04[0] = Shower_width(Energies4, Cellids4, 3)
        w21st_l04[0] = Shower_width(Energies4, Cellids4, 21)
        eocore_l04[0] = eocorey(cal4Emax, Cellids4, Energies4, 3)
        e2max_l04[0] = 0
        emax_l04[0] = cal4Emax
        edmax_l04[0] = 0

    if Elayer1 > 0.1:
        e_l01[0] = Elayer0/Elayer1
        e_l20[0] = Elayer2/Elayer1
        e_l30[0] = Elayer3/Elayer1
        e_l40[0] = Elayer4/Elayer1
        e_l50[0] = Elayer5/Elayer1
        e_l60[0] = Elayer6/Elayer1
        e_l70[0] = Elayer7/Elayer1

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
