import numpy as n
import ROOT as r
import math as m
from array import array

Types = ["basic_","fineEta2ndLayer_","half2ndLayer_","half2nd3rdLayer_","half2nd3rdLayerSmallPhiAllLayers_"]

Training_Energies = n.array([10,20,30,40,50])
TrainPath = "/afs/cern.ch/user/r/rastein/TMVAFacility/TMVA_BDT_"
TrainType = "pi0_vs_photon_"
TrainMethod = "BDT_"
unit = "GeV"


def ExtractBGRejection(S, Energies,type):
    Pi0Rej = []
    Errx = []
    Erry1 = []
    Erry2 = []
    for e in Energies:
        ene = str(e)
        file = TrainPath+TrainType+type+ene+unit+".root"
    #    print file
        infile = r.TFile.Open(file)
        dir = TrainType+type+ene+unit+"/Method_"+TrainMethod+TrainType+type+ene+unit+"/"+TrainMethod+TrainType+type+ene+unit
        hist = infile.Get(dir+"/MVA_"+TrainMethod+TrainType+type+ene+unit+"_rejBvsS")
        hist.Draw()
        x = hist.GetXaxis()
        ID = x.FindBin(S)
        diff = n.float(x.GetBinWidth(ID))
        rejection = hist.GetBinContent(ID)
        y1 = hist.GetBinContent(x.FindBin(S+diff))
        y_1 = hist.GetBinContent(x.FindBin(S-diff))
        if rejection >= 1:
            rejection = 0.99
        if y1 >= 1:
            y1 = 0.99
        if y_1 >= 1:
            y_1 = 0.99
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


N = len(Training_Energies)
y = []
x = Training_Energies
ex = []
eyl = []
eyh = []
sigerr = []
for type in Types:
    print type
    y.append(ExtractBGRejection(0.9,Training_Energies,type)[0])
    ex.append(ExtractBGRejection(0.9,Training_Energies,type)[1])
    eyl.append(ExtractBGRejection(0.9,Training_Energies,type)[2])
    eyh.append(ExtractBGRejection(0.9,Training_Energies,type)[3])
    #y2 = ExtractBGRejection(0.8,Training_Energies)[0]
    #ex2 = ExtractBGRejection(0.8,Training_Energies)[1]
    #ey2l = ExtractBGRejection(0.8,Training_Energies)[2]
    #ey2h = ExtractBGRejection(0.8,Training_Energies)[3]
    sigerr.append(int(ExtractBGRejection(0.9,Training_Energies,type)[4]*100))
    #sigerr2 = int(ExtractBGRejection(0.8,Training_Energies)[4]*100)

canvas2 = r.TCanvas( 'canvas', 'Pion rejection at different energies', 200, 10, 700, 500 )
canvas2.SetTicks(1,1)
canvas2.SetLeftMargin(0.14)
canvas2.SetRightMargin(0.08)
canvas2.SetGrid()
r.gStyle.SetOptStat(0)
#canvas2.SetGrid()
canvas2.SetFillColor(0)
canvas2.GetFrame().SetBorderSize(12)

gr1 = r.TGraphAsymmErrors(N,array('d',x),array('d',y[0]),array('d',ex[0]),array('d',ex[0]),array('d',eyl[0]),array('d',eyh[0]))
gr1.SetLineWidth(2)

gr1.SetLineColorAlpha(r.kBlue+1,0.5)
gr1.SetMarkerColor(r.kBlue)

gr1.SetMarkerStyle(21)

gr1.SetMarkerSize(1)
#fr1 = pad1.DrawFrame(5.0,0.0,55.0,5.5)
#fr1.SetXTitle("p_T (GeV)")
#fr1.SetYTitle("Pion Rejection Factor")
#pad1.GetFrame().SetFillColor(10)
#pad1.GetFrame().SetBorderSize(12)
gr1.Draw("ACP")


gr1.GetYaxis().SetTitleOffset(1.4)
gr1.GetXaxis().SetTitleOffset(1.40)
gr1.SetTitle("")

xax = r.TPaveText()
gr1.GetXaxis().SetTitle("Transverse momentum, p_{T} (GeV)")
gr1.GetYaxis().SetTitle("Pion Rejection Factor, R")

#maxh = gr1.GetMaximum()
#minh = gr2.GetMinimum()

gr1.SetMaximum(12.0)
gr1.SetMinimum(1.0)


line = r.TLine(0,3,60,3)
line.SetLineColor(r.kRed)
line.SetLineWidth(2)
line.Draw()

'''canvas2.cd()
pad2 = r.TPad("overlay","",0,0,1,1)
pad2.SetFillStyle(4000)
pad2.SetFillColor(0)
pad2.SetFrameFillStyle(4000)
pad2.Draw()
pad2.cd()'''

gr2 = r.TGraphAsymmErrors(N,array('d',x),array('d',y[1]),array('d',ex[1]),array('d',ex[1]),array('d',eyl[1]),array('d',eyh[1]))
gr2.SetLineWidth(2)
gr2.SetLineColorAlpha(r.kOrange+1,0.5)
gr2.SetMarkerColor(r.kOrange)

gr2.SetMarkerStyle(20)
gr2.SetMarkerSize(1)

#fr2 = pad2.DrawFrame(5.0,0.0,55.0,5.5)
gr2.Draw("CP")

gr3 = r.TGraphAsymmErrors(N,array('d',x),array('d',y[2]),array('d',ex[2]),array('d',ex[2]),array('d',eyl[2]),array('d',eyh[2]))
gr3.SetLineWidth(2)
gr3.SetLineColorAlpha(r.kMagenta+1,0.5)
gr3.SetMarkerColor(r.kMagenta)

gr3.SetMarkerStyle(22)
gr3.SetMarkerSize(1)

#fr2 = pad2.DrawFrame(5.0,0.0,55.0,5.5)
gr3.Draw("CP")

gr4 = r.TGraphAsymmErrors(N,array('d',x),array('d',y[3]),array('d',ex[3]),array('d',ex[3]),array('d',eyl[3]),array('d',eyh[3]))
gr4.SetLineWidth(2)
gr4.SetLineColorAlpha(r.kGreen+1,0.5)
gr4.SetMarkerColor(r.kGreen)

gr4.SetMarkerStyle(33)
gr4.SetMarkerSize(1)

#fr2 = pad2.DrawFrame(5.0,0.0,55.0,5.5)
gr4.Draw("CP")
#axis = r.TGaxis(5.0,0.0,55.0,6.0,"+L")
#axis.Draw()

gr5 = r.TGraphAsymmErrors(N,array('d',x),array('d',y[4]),array('d',ex[4]),array('d',ex[4]),array('d',eyl[4]),array('d',eyh[4]))
gr5.SetLineWidth(2)
gr5.SetLineColorAlpha(r.kGray+2,0.5)
gr5.SetMarkerColor(r.kBlack)

gr5.SetMarkerStyle(34)
gr5.SetMarkerSize(1)

#fr2 = pad2.DrawFrame(5.0,0.0,55.0,5.5)
gr5.Draw("CP")

leg = r.TLegend(0.67,0.68,0.9,0.86)

leg.SetFillColor(0)
leg.SetFillStyle(1001)
leg.SetLineColor(0)
leg.SetShadowColor(10)
leg.SetTextSize(0.035)
leg.SetTextFont(42)
leg.AddEntry(gr1,"Sample 1", "lp")#Delta#eta = 0.01","lp")
leg.AddEntry(gr2,"Sample 2", "lp")#Delta#eta = 0.0025 in 2^{nd} layer","lp")
leg.AddEntry(gr3,"Sample 3", "lp")#+ 2^{nd} layer halved","lp")
leg.AddEntry(gr4,"Sample 4", "lp")#+ 3^{rd} layer halved","lp")
leg.AddEntry(gr5,"Sample 5", "lp")#+ #Delta#phi = 0.0045","lp")
leg.AddEntry(line,"Threshold R_{#pi} = 3","l")
leg.Draw("same")

text = "90% #pm "+str(sigerr[0])+"% Signal Efficiency"
Text = r.TLatex()
Text.SetNDC(r.kTRUE)
Text.SetTextAlign(0);
Text.SetTextSize(0.035)

Text.DrawLatex(0.22, 0.82, text)
canvas2.RedrawAxis()

canvas2.Modified()
canvas2.Update()
canvas2.Print("Superimposed_Geometries.pdf")
