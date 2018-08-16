import numpy as n
import ROOT as r
import math as m
from array import array

Vars = ["5var","10var","15var","20var","25var","30var","35var","40var","45var"]
vars = [5,10,15,20,25,30,35,40,45]
#Variable_names = ["5var","7var","10var","13var","15var","18var","21var","25var","30var","35var","41var"]
Times = [79.83,90.86,123.29,147.28,170.63,189.43,198.55,219.27,250.43]
#Variables = dict(zip(Vars,Variable_names))


Training_Energies = n.array([50])
TrainPath = "/afs/cern.ch/user/r/rastein/TMVAFacility/TMVA_BDT_"
TrainType = "pi0_vs_photon_half2nd3rdLayerSmallPhiAllLayers_"
TrainMethod = "BDT_"
unit = "GeV_"


def ExtractBGRejection(S, Energies,type):
    Pi0Rej = []
    Errx = []
    Erry1 = []
    Erry2 = []
    for e in Energies:
        ene = str(e)
        file = TrainPath+TrainType+ene+unit+type+".root"
    #    print file
        infile = r.TFile.Open(file)
        dir = TrainType+ene+unit+type+"/Method_"+TrainMethod+TrainType+ene+unit+type+"/"+TrainMethod+TrainType+ene+unit+type
        hist = infile.Get(dir+"/MVA_"+TrainMethod+TrainType+ene+unit+type+"_rejBvsS")
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
        Errx.append(n.float(0))
        Pi0Rej.append(R)
    #n.asarray(BGRej)
    #n.asarray(Errx)
    #n.asarray(Erry)
    return Pi0Rej,Errx,Erry1,Erry2,diff

N = len(vars)
y = []
xs = vars
y2 = Times
ex = []
eyl = []
eyh = []
sigerr = []

for var in Vars:
    print var
    y.append(ExtractBGRejection(0.9,Training_Energies,var)[0][0])
    ex.append(ExtractBGRejection(0.9,Training_Energies,var)[1][0])
    eyl.append(ExtractBGRejection(0.9,Training_Energies,var)[2][0])
    eyh.append(ExtractBGRejection(0.9,Training_Energies,var)[3][0])
    #y2 = ExtractBGRejection(0.8,Training_Energies)[0]
    #ex2 = ExtractBGRejection(0.8,Training_Energies)[1]
    #ey2l = ExtractBGRejection(0.8,Training_Energies)[2]
    #ey2h = ExtractBGRejection(0.8,Training_Energies)[3]
    sigerr.append(int(ExtractBGRejection(0.9,Training_Energies,var)[4]*100))
    #sigerr2 = int(ExtractBGRejection(0.8,Training_Energies)[4]*100)

canvas2 = r.TCanvas( 'canvas', '#scale[1.1]{#pi_{0}} rejection for different no. of variables', 200, 10, 700, 500 )
canvas2.SetTicks(1,1)
canvas2.SetLeftMargin(0.14)
canvas2.SetRightMargin(0.08)
r.gStyle.SetOptStat(0)
canvas2.SetGrid()
canvas2.SetFillColor(0)
canvas2.GetFrame().SetBorderSize(12)

pad1 = r.TPad("pad1","",0,0,1,1)
pad1.SetGrid()
pad1.SetFillColor(0)
pad1.Draw()
pad1.cd()
fr1 = pad1.DrawFrame(0.0,1.5,45,4.)
fr1.SetXTitle("Transverse Momentum, p_{T} (GeV)")
fr1.SetYTitle("#scale[1.1]{#pi_{0}} Rejection Factor, R")
pad1.GetFrame().SetFillColor(10)
pad1.GetFrame().SetBorderSize(12)

gr1 = r.TGraphAsymmErrors(N,array('d',xs),array('d',y),array('d',ex),array('d',ex),array('d',eyl),array('d',eyh))
gr1.SetLineWidth(2)

gr1.SetLineColorAlpha(r.kBlack,0.5)
gr1.SetMarkerColor(r.kBlue)

gr1.SetMarkerStyle(21)

gr1.SetMarkerSize(1)

gr1.Draw("ACP")


gr1.GetYaxis().SetTitleOffset(1.4)
gr1.GetXaxis().SetTitleOffset(1.40)
gr1.SetTitle("")

xax = r.TPaveText()
gr1.GetXaxis().SetTitle("No. of variables, N")
gr1.GetYaxis().SetTitle("#scale[1.1]{#pi_{0}} Rejection Factor, R")

#maxh = gr1.GetMaximum()
#minh = gr2.GetMinimum()

gr1.SetMaximum(5.0)
gr1.SetMinimum(2.5)

canvas2.cd()
pad2 = r.TPad("overlay","",0,0,1,1)
pad2.SetFillStyle(4000)
pad2.SetFillColor(0)
pad2.SetFrameFillStyle(4000)
pad2.Draw("FA")
pad2.cd()

gr2 = r.TGraphAsymmErrors(N,array('d',xs),array('d',y2),array('d',ex),array('d',ex),array('d',ex),array('d',ex))
gr2.SetLineWidth(2)
gr2.SetLineColorAlpha(r.kOrange+1,0.5)
gr2.SetMarkerColor(r.kRed+1)

gr2.SetMarkerStyle(20)
gr2.SetMarkerSize(1)
rightmax = 1.1*gr2.GetMaximum()
rightmin = gr2.GetMinimum()/1.1
scale = r.gPad.GetUymax()/rightmax
#gr2.Scale(scale)
gr2.SetMaximum(400.0)
gr2.SetMinimum(50.0)
fr2 = pad2.DrawFrame(1,50,49,350)
fr2.GetXaxis().SetLabelOffset(99)
fr2.GetYaxis().SetTickLength(0)
fr2.GetXaxis().SetTickLength(0)
fr2.GetYaxis().SetLabelOffset(99)
gr2.Draw("CP")
axis = r.TGaxis(49,50,49,350,50,350,510,"+L")
axis.SetLineColor(r.kOrange+1)
axis.SetLabelColor(r.kOrange+1)
axis.SetTitleOffset(1.40)
axis.SetTitle("TMVA Run Time, t (Seconds)")
axis.Draw()

'''gr3 = r.TGraphAsymmErrors(N,array('d',x),array('d',y[2]),array('d',ex[2]),array('d',ex[2]),array('d',eyl[2]),array('d',eyh[2]))
gr3.SetLineWidth(2)
gr3.SetLineColorAlpha(r.kMagenta,0.5)
gr3.SetMarkerColor(r.kMagenta+1)

gr3.SetMarkerStyle(22)
gr3.SetMarkerSize(1)

#fr2 = pad2.DrawFrame(5.0,0.0,55.0,5.5)
gr3.Draw("LP")'''

'''gr4 = r.TGraphAsymmErrors(N,array('d',x),array('d',y[3][0]),array('d',ex[3][1]),array('d',ex[3][1]),array('d',eyl[3][2]),array('d',eyh[3][3]))
gr4.SetLineWidth(2)
gr4.SetLineColorAlpha(r.kMagenta,0.5)
gr4.SetMarkerColor(r.Magenta+1)

gr4.SetMarkerStyle(33)
gr4.SetMarkerSize(1)

#fr2 = pad2.DrawFrame(5.0,0.0,55.0,5.5)
gr4.Draw("LP")'''
#axis = r.TGaxis(5.0,0.0,55.0,6.0,"+L")
#axis.Draw()
leg = r.TLegend(0.60,0.75,0.84,0.88)
leg.SetFillColor(0)
leg.SetFillStyle(1001)
leg.SetLineColor(0)
leg.SetShadowColor(10)
leg.SetTextSize(0.04)
leg.SetTextFont(42)
leg.AddEntry(gr1,"#scale[1.1]{#pi_{0}} Rejection factor","lp")
leg.AddEntry(gr2,"Time to run TMVA","lp")
leg.Draw("same")

text = "90% #pm "+str(sigerr[0])+"% Signal Efficiency"
Text = r.TLatex()
Text.SetNDC(r.kTRUE)
Text.SetTextAlign(0);
Text.SetTextSize(0.035)
canvas2.RedrawAxis()

canvas2.Modified()
canvas2.Update()
canvas2.Print("Variable_Numbers.pdf")
