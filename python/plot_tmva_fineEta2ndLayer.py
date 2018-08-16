import numpy as n
import ROOT as r
import math as m
from array import array

Training_Energies = n.array([10,20,30,40,50])
TrainPath = "/afs/cern.ch/user/r/rastein/TMVAFacility/TMVA_BDT_"
TrainType = "pi0_vs_photon_fineEta2ndLayer_"
TrainMethod = "BDT_"
unit = "GeV"


def ExtractBGRejection(S, Energies):
    Pi0Rej = []
    Errx = []
    Erry1 = []
    Erry2 = []
    for e in Energies:
        ene = str(e)
        file = TrainPath+TrainType+ene+unit+".root"
        print file
        infile = r.TFile.Open(file)
        dir = TrainType+ene+unit+"/Method_"+TrainMethod+TrainType+ene+unit+"/"+TrainMethod+TrainType+ene+unit
        hist = infile.Get(dir+"/MVA_"+TrainMethod+TrainType+ene+unit+"_rejBvsS")
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
        #print Ry1,Ry_1
        Erry1.append(n.float((R-Ry1)))
        Erry2.append(n.float((Ry_1-R)))
        Errx.append(n.float(5))
        Pi0Rej.append(R)
    #n.asarray(BGRej)
    #n.asarray(Errx)
    #n.asarray(Erry)
    return Pi0Rej,Errx,Erry1,Erry2,diff

N = len(Training_Energies)
x = Training_Energies
y1 = ExtractBGRejection(0.8,Training_Energies)[0]
ex1 = ExtractBGRejection(0.8,Training_Energies)[1]
ey1l = ExtractBGRejection(0.8,Training_Energies)[2]
ey1h = ExtractBGRejection(0.8,Training_Energies)[3]
y2 = ExtractBGRejection(0.9,Training_Energies)[0]
ex2 = ExtractBGRejection(0.9,Training_Energies)[1]
ey2l = ExtractBGRejection(0.9,Training_Energies)[2]
ey2h = ExtractBGRejection(0.9,Training_Energies)[3]
sigerr1 = int(ExtractBGRejection(0.8,Training_Energies)[4]*100)
sigerr2 = int(ExtractBGRejection(0.9,Training_Energies)[4]*100)

canvas2 = r.TCanvas( 'canvas', 'Pion rejection at different energies', 200, 10, 700, 500 )
canvas2.SetTicks(1,1)
canvas2.SetLeftMargin(0.14)
canvas2.SetRightMargin(0.08)
r.gStyle.SetOptStat(0)
canvas2.SetGrid()
canvas2.SetFillColor(0)
canvas2.GetFrame().SetBorderSize(12)

gr1 = r.TGraphAsymmErrors(N,array('d',x),array('d',y2),array('d',ex1),array('d',ex1),array('d',ey2l),array('d',ey2h))
gr1.SetLineWidth(2)

gr1.SetLineColorAlpha(r.kBlack,0.5)
gr1.SetMarkerColor(r.kBlue)

gr1.SetMarkerStyle(21)

gr1.SetMarkerSize(1)
#fr1 = pad1.DrawFrame(5.0,0.0,55.0,5.5)
#fr1.SetXTitle("p_T (GeV)")
#fr1.SetYTitle("Pion Rejection Factor")
#pad1.GetFrame().SetFillColor(10)
#pad1.GetFrame().SetBorderSize(12)
gr1.Draw("ALP")


gr1.GetYaxis().SetTitleOffset(1.4)
gr1.GetXaxis().SetTitleOffset(1.40)
gr1.SetTitle("fineEta2ndLayer")

xax = r.TPaveText()
gr1.GetXaxis().SetTitle("p_T (GeV)")
gr1.GetYaxis().SetTitle("Pion Rejection Factor")

#maxh = gr1.GetMaximum()
#minh = gr2.GetMinimum()

gr1.SetMaximum(9.0)
gr1.SetMinimum(2.0)




'''canvas2.cd()
pad2 = r.TPad("overlay","",0,0,1,1)
pad2.SetFillStyle(4000)
pad2.SetFillColor(0)
pad2.SetFrameFillStyle(4000)
pad2.Draw()
pad2.cd()'''

gr2 = r.TGraphAsymmErrors(N,array('d',x),array('d',y1),array('d',ex2),array('d',ex2),array('d',ey1l),array('d',ey1h))
gr2.SetLineWidth(2)
gr2.SetLineColorAlpha(r.kBlack,0.5)
gr2.SetMarkerColor(r.kOrange+1)

gr2.SetMarkerStyle(20)
gr2.SetMarkerSize(1)

#fr2 = pad2.DrawFrame(5.0,0.0,55.0,5.5)
gr2.Draw("LP")

#axis = r.TGaxis(5.0,0.0,55.0,6.0,"+L")
#axis.Draw()
leg = r.TLegend(0.52,0.72,0.92,0.90)
leg.SetFillColorAlpha(r.kYellow-10,0.3)
leg.SetFillStyle(3010)
leg.SetLineColor(1)
leg.SetShadowColor(31)
leg.SetTextSize(0.035)
leg.SetTextFont(42)
leg.AddEntry(gr1,"90% +/- "+str(sigerr2)+"% Signal Efficiency","lp")
leg.AddEntry(gr2,"80% +/- "+str(sigerr1)+"% Signal Efficiency","lp")
leg.Draw("same")

canvas2.RedrawAxis()

canvas2.Modified()
canvas2.Update()
