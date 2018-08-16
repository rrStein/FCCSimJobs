import numpy as n
import ROOT as r
import math as m
from array import array

Variables = ["emax_l01", "e2max_l01", "eocore_l01", "edmax_l01", "w9st_l01", "w41st_l01", "e_l21", "e_l1T"]
particles = ["pi0", "photon"]
TrainPath = "/afs/cern.ch/user/r/rastein/TMVAFacility/TMVA_BDT_"
TrainType = "pi0_vs_photon_half2nd3rdLayerSmallPhiAllLayers_"
TrainMethod = "BDT_"
Energy = "50"
unit = "GeV"

for var in Variables:
    if var == "emax_l01":
        name = "E_{max}"
    elif var == "e2max_l01":
        name = "E_{2ndmax}"
    elif var == "eocore_l01":
        name = "E_{ocore}"
    elif var == "edmax_l01":
        name = "E_{dmax}"
    elif var == "w9st_l01":
        name = "W_{9st}"
    elif var == "w41st_l01":
        name = "W41_{41st}"
    elif var == "e_l21":
        name = "E_{21}"
    elif var == "e_l1T":
        name = "E_{1T}"
    #if var == "w41st_l01":
        #file = TrainPath+TrainType+"70GeV.root"
        #Energy = "70"
        #unit ="GeV"
    #else:
    file = TrainPath+TrainType+Energy+unit+".root"
    infile = r.TFile.Open(file)
    dir_phot = TrainType+Energy+unit+"/InputVariables_Id/"+var+"__Signal_Id"
    dir_pi0 = TrainType+Energy+unit+"/InputVariables_Id/"+var+"__Background_Id"
    hist1 = infile.Get(dir_phot)
    hist2 = infile.Get(dir_pi0)


    logY = False
    if var == "edmax_l01": #or var == "e2max_l01":
        logY = True
    hist1.SetLineWidth(3)
    hist2.SetLineWidth(3)

    hist1.SetLineColor(r.kBlue+1)
    hist2.SetLineColor(r.kOrange+1)

    hist1.SetFillColorAlpha(r.kBlue, 0.5)
    hist2.SetFillColorAlpha(r.kOrange, 0.5)

    scale = hist1.GetEntries()
    #hist1.Scale(scale)
    #hist2.Scale(scale)

    canvas = r.TCanvas(var, var, 800, 600)
    canvas.SetLogy(logY)
    canvas.SetTicks(1,1)
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.08)
    r.gStyle.SetOptStat(0)

    leg = r.TLegend(0.50,0.62,0.80,0.86)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetShadowColor(10)
    leg.SetTextSize(0.04)
    leg.SetTextFont(42)
    leg.AddEntry(hist1,"#splitline{#scale[1.1]{#bf{#gamma}} (p_{T} = 50 GeV),}{mean = "+str(round(hist1.GetMean(),2))+" GeV}","l")
    leg.AddEntry(hist2,"#splitline{#scale[1.1]{#bf{#pi_{0}}} (p_{T} = 50 GeV),}{mean = "+str(round(hist2.GetMean(),2))+" GeV}","l")

    hist1.GetYaxis().SetTitleOffset(1.95)
    hist1.GetXaxis().SetTitleOffset(1.40)
    hist1.SetTitle("")
    hist1.GetXaxis().SetTitle("Energy (GeV)")
    hist1.GetYaxis().SetTitle("Entries")

    #hist1.Scale(lumi)
    #hist2.Scale(lumi)

    maxh1 = hist1.GetMaximum()
    minh1 = hist1.GetMinimum()
    maxh2 = hist2.GetMaximum()
    minh2 = hist2.GetMinimum()

    if logY:
        if maxh1 > maxh2:
            hist1.SetMaximum(1000*maxh1)
            hist2.SetMaximum(1000*maxh1)
        else:
            hist1.SetMaximum(1000*maxh2)
            hist2.SetMaximum(1000*maxh2)
    else:
        if maxh1 > maxh2:
            hist1.SetMaximum(1.2*maxh1)
            hist1.SetMinimum(0.)
            hist2.SetMaximum(1.2*maxh1)
            hist2.SetMinimum(0.)
        else:
            hist1.SetMaximum(1.2*maxh2)
            hist1.SetMinimum(0.)
            hist2.SetMaximum(1.2*maxh2)
            hist2.SetMinimum(0.)

    hist1.Draw()
    hist2.Draw("same")
    leg.Draw("same")

    Text = r.TLatex()
    Text.SetNDC(r.kTRUE)
    Text.SetTextAlign(0);
    Text.SetTextSize(0.055)

    Text.DrawLatex(0.65, 0.5, name)

    canvas.RedrawAxis()
    canvas.GetFrame().SetBorderSize( 12 )
    canvas.Modified()
    canvas.Update()
    canvas.Print(var+".pdf")
