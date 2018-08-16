import numpy as n
import ROOT as r
import math as m
from array import array

Energies = n.array([''])
Etas = n.array([-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75])
TrainPath = "/afs/cern.ch/user/r/rastein/TMVAFacility/TMVA_BDT_"
TrainType = "pi0_vs_photon_half2nd3rdLayerSmallPhiAllLayers_"
TrainMethod = "BDT_"
unit = ""
Efficiencies = n.array([0.9])

N = len(Etas)
def ExtractBGRejection(S, Eta):
    pi0Rej = []
    yh = []
    yl = []
    x_ = []
    #for e in Etas:
    #for S in Efficiencies:
    ene = ""
    if Eta == -1.0 or Eta == 0.:
        Eta = int(Eta)


    eta = "eta"+str(Eta)
    if Eta == 0:
        eta = "All15var"
    file = TrainPath+TrainType+ene+unit+eta+".root"
    #print file
    infile = r.TFile.Open(file)
    dir = TrainType+ene+unit+eta+"/Method_"+TrainMethod+TrainType+ene+unit+eta+"/"+TrainMethod+TrainType+ene+unit+eta
    hist = infile.Get(dir+"/MVA_"+TrainMethod+TrainType+ene+unit+eta+"_rejBvsS")
    hist.Draw()
    x = hist.GetXaxis()
    ID = x.FindBin(S)
    diff = n.float(x.GetBinWidth(ID))
    rejection = hist.GetBinContent(ID)
    y1 = hist.GetBinContent(x.FindBin(S+diff))
    y_1 = hist.GetBinContent(x.FindBin(S-diff))
    R = 1/(1-rejection)
    Ry1 = 1/(1-(y1))
    Ry_1 = 1/(1-(y_1))
    pi0Rej.append(R)
    yh.append(Ry_1-R)
    yl.append(R-Ry1)
    x_.append(0.25)
    return R,0.25,Ry_1-R,R-Ry1,diff
            #print "The pi0 rejection for energy: ",e," is: ", round(R,2)," +/- ",round((Ry_1 - Ry1)/2,2)," at ",S*100,"%"," signal efficiency"
            #print Ry1,Ry_1
            #n.asarray(BGRej)
            #n.asarray(Errx)
            #n.asarray(Erry)
x = Etas
N = len(Etas)
y = []
y_err_1 = []
y_err1 = []
x_err = []
sigerr = []
for e in Etas:
    print e
    y.append(ExtractBGRejection(0.9,e)[0])
    x_err.append(ExtractBGRejection(0.9,e)[1])
    y_err1.append(ExtractBGRejection(0.9,e)[2])
    y_err_1.append(ExtractBGRejection(0.9,e)[3])
    #y2 = ExtractBGRejection(0.8,Training_Energies)[0]
    #ex2 = ExtractBGRejection(0.8,Training_Energies)[1]
    #ey2l = ExtractBGRejection(0.8,Training_Energies)[2]
    #ey2h = ExtractBGRejection(0.8,Training_Energies)[3]
    sigerr.append(int(ExtractBGRejection(0.9,e)[4]*100))
    #sigerr2 = int(ExtractBGRejection(0.8,Training_Energies)[4]*100)

canvas2 = r.TCanvas( 'canvas', 'Pion rejection at different {#eta}', 200, 10, 700, 500 )
canvas2.SetTicks(1,1)
canvas2.SetLeftMargin(0.14)
canvas2.SetRightMargin(0.08)
r.gStyle.SetOptStat(0)
canvas2.SetGrid()
canvas2.SetFillColor(0)
canvas2.GetFrame().SetBorderSize(12)

gr1 = r.TGraphAsymmErrors(N,array('d',x),array('d',y),array('d',x_err),array('d',x_err),array('d',y_err_1),array('d',y_err1))
gr1.SetLineWidth(2)

gr1.SetLineColorAlpha(r.kBlack,0.5)
gr1.SetMarkerColor(r.kBlue)

gr1.SetMarkerStyle(2)

gr1.SetMarkerSize(2)
#fr1 = pad1.DrawFrame(5.0,0.0,55.0,5.5)
#fr1.SetXTitle("p_T (GeV)")
#fr1.SetYTitle("Pion Rejection Factor")
#pad1.GetFrame().SetFillColor(10)
#pad1.GetFrame().SetBorderSize(12)
gr1.Draw("AP")


gr1.GetYaxis().SetTitleOffset(1.4)
gr1.GetXaxis().SetTitleOffset(1.40)
gr1.SetTitle("")

xax = r.TPaveText()
gr1.GetXaxis().SetTitle("Pseudorapidity, #scale[1.1]{#eta}")
gr1.GetYaxis().SetTitle("#scale[1.1]{#pi_{0}} Rejection Factor, R")

#maxh = gr1.GetMaximum()
#minh = gr2.GetMinimum()

gr1.SetMaximum(5.0)
gr1.SetMinimum(2.0)

line = r.TLine(-2,3,2,3)
line.SetLineColor(r.kRed)
line.SetLineWidth(2)
line.Draw()
