import uproot
import ROOT
import awkward as ak
import numpy as np

fname = "/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ekauffma/JetMET1/crab_JetMET_2023C_v2/230620_114703/combined_histo_368382.root"

f = uproot.open(fname)
f.keys()

inFile = ROOT.TFile.Open( fname ,"READ")

##############################################################################
################################ jet et plots ################################
##############################################################################

# unprefirable jet et
c1 = ROOT.TCanvas( 'c1', 'c1', 100, 10, 600, 600 )
JetEtbx0_u = inFile.Get("demo/JetEtbx0_unprefirable")
JetEtbx0_u.SetTitle("Unprefirable")
JetEtbx0_u.SetLineColor(1)
JetEtbx0_u.GetXaxis().SetTitle("L1 Jet E_{T} (GeV)")
JetEtbx0_u.GetYaxis().SetTitle("a.u.")
JetEtbx0_u.SetStats(False)
JetEtbx0_u.SetMaximum(1e5)
JetEtbx0_u.SetMinimum(1e1)
JetEtbx0_u.Draw("")

JetEtbxm1_u = inFile.Get("demo/JetEtbxm1_unprefirable")
JetEtbxm1_u.SetStats(False)
JetEtbxm1_u.SetLineColor(2)
JetEtbxm1_u.Draw("same")

legend = ROOT.TLegend(0.6,0.7,0.9,0.8)
legend.AddEntry(JetEtbx0_u,"BX=0","L")
legend.AddEntry(JetEtbxm1_u,"BX=-1","L")
legend.SetBorderSize(0);
legend.SetTextSize(0.04);
legend.Draw()

c1.SetLogy();
c1.Draw()
c1.SaveAs(f"plots/jetet_unprefirable.png")
c1.Close()

# firstbunchintrain jet et
c = ROOT.TCanvas( 'c', 'c', 100, 10, 600, 600 )
JetEtbx0_f = inFile.Get("demo/JetEtbx0_firstbunch")
JetEtbx0_f.SetTitle("FirstBunchInTrain")
JetEtbx0_f.SetStats(False)
JetEtbx0_f.SetLineColor(1)
JetEtbx0_f.SetMaximum(1e5)
JetEtbx0_f.SetMinimum(1e1)
JetEtbx0_f.Draw("")

JetEtbxm1_f = inFile.Get("demo/JetEtbxm1_firstbunch")
JetEtbxm1_f.SetStats(False)
JetEtbxm1_f.SetLineColor(2)
JetEtbxm1_f.Draw("same")

legend2 = ROOT.TLegend(0.6,0.7,0.9,0.8)
legend2.AddEntry(JetEtbx0_f,"BX=0","L")
legend2.AddEntry(JetEtbxm1_f,"BX=-1","L")
legend2.SetBorderSize(0);
legend2.SetTextSize(0.04);
legend2.Draw()

c.SetLogy();
c.Draw()
c.SaveAs(f"plots/jetet_firstbunch.png")
c.Close()

# ratio jet et
c1_ratio = ROOT.TCanvas( 'c1_ratio', 'c1_ratio', 100, 10, 600, 600 )

JetEtbx0_u = inFile.Get("demo/JetEtbx0_unprefirable")
JetEtbxm1_u = inFile.Get("demo/JetEtbxm1_unprefirable")

JetEtcombined_u = JetEtbx0_u.Clone()
JetEtcombined_u.Add(JetEtbxm1_u)

JetEtratio_u = JetEtbxm1_u.Clone()
JetEtratio_u.Divide(JetEtcombined_u)

JetEtratio_u.SetLineColor(50)
JetEtratio_u.SetTitle("")
JetEtratio_u.GetXaxis().SetTitle("L1 Jet E_{T} (GeV)")
JetEtratio_u.GetYaxis().SetTitle("(bx=-1)/(bx=-1 or bx=0)")
JetEtratio_u.Draw("")

JetEtbx0_f = inFile.Get("demo/JetEtbx0_firstbunch")
JetEtbxm1_f = inFile.Get("demo/JetEtbxm1_firstbunch")

JetEtcombined_f = JetEtbx0_f.Clone()
JetEtcombined_f.Add(JetEtbxm1_f)

JetEtratio_f = JetEtbxm1_f.Clone()
JetEtratio_f.Divide(JetEtcombined_f)

JetEtratio_f.SetLineColor(9)
JetEtratio_f.Draw("same")

legend = ROOT.TLegend(0.4,0.7,0.9,0.8)
legend.AddEntry(JetEtratio_u,"UnprefirableEvent","L")
legend.AddEntry(JetEtratio_f,"FirstBunchInTrain","L")
legend.SetBorderSize(0);
legend.SetTextSize(0.04);
legend.Draw()

c1_ratio.Draw()
c1_ratio.SaveAs(f"plots/jetet_ratio.png")
c1_ratio.Close()


#############################################################################
############################# lead jet et plots #############################
#############################################################################

# unprefirable lead jet et
c2 = ROOT.TCanvas( 'c2', 'c2', 100, 10, 600, 600 )
LeadJetEtbx0 = inFile.Get("demo/LeadJetEtbx0_unprefirable")
LeadJetEtbx0.SetTitle("Unprefirable Event")
LeadJetEtbx0.SetLineColor(1)
LeadJetEtbx0.GetXaxis().SetTitle("L1 Leading Jet E_{T} (GeV)")
LeadJetEtbx0.GetYaxis().SetTitle("a.u.")
LeadJetEtbx0.SetStats(False)
LeadJetEtbx0.SetMaximum(1e5)
LeadJetEtbx0.Draw("")

LeadJetEtbxm1 = inFile.Get("demo/LeadJetEtbxm1_unprefirable")
LeadJetEtbxm1.SetStats(False)
LeadJetEtbxm1.SetLineColor(2)
LeadJetEtbxm1.Draw("same")

LeadJetEtbxm2 = inFile.Get("demo/LeadJetEtbxm2_unprefirable")
LeadJetEtbxm2.SetStats(False)
LeadJetEtbxm2.SetLineColor(4)
LeadJetEtbxm2.Draw("same")

legend = ROOT.TLegend(0.6,0.6,0.9,0.85)
legend.AddEntry(LeadJetEtbx0,"BX=0","L")
legend.AddEntry(LeadJetEtbxm1,"BX=-1","L")
legend.AddEntry(LeadJetEtbxm2,"BX=-2","L")
legend.SetBorderSize(0);
legend.SetTextSize(0.04);
legend.Draw()

c2.SetLogy();
c2.Draw()
c2.SaveAs(f"plots/leadjetet_unprefirable.png")
c2.Close()


# firstbunchintrain lead jet et
c20 = ROOT.TCanvas( 'c20', 'c20', 100, 10, 600, 600 )
LeadJetEtbx0 = inFile.Get("demo/LeadJetEtbx0_firstbunch")
LeadJetEtbx0.SetTitle("FirstBunchInTrain")
LeadJetEtbx0.SetLineColor(1)
LeadJetEtbx0.GetXaxis().SetTitle("L1 Leading Jet E_{T} (GeV)")
LeadJetEtbx0.GetYaxis().SetTitle("a.u.")
LeadJetEtbx0.SetStats(False)
LeadJetEtbx0.SetMaximum(1e5)
LeadJetEtbx0.Draw("")

LeadJetEtbxm1 = inFile.Get("demo/LeadJetEtbxm1_firstbunch")
LeadJetEtbxm1.SetStats(False)
LeadJetEtbxm1.SetLineColor(2)
LeadJetEtbxm1.Draw("same")

LeadJetEtbxm2 = inFile.Get("demo/LeadJetEtbxm2_firstbunch")
LeadJetEtbxm2.SetStats(False)
LeadJetEtbxm2.SetLineColor(4)
LeadJetEtbxm2.Draw("same")

legend = ROOT.TLegend(0.6,0.6,0.9,0.85)
legend.AddEntry(LeadJetEtbx0,"BX=0","L")
legend.AddEntry(LeadJetEtbxm1,"BX=-1","L")
legend.AddEntry(LeadJetEtbxm2,"BX=-2","L")
legend.SetBorderSize(0);
legend.SetTextSize(0.04);
legend.Draw()

c20.SetLogy();
c20.Draw()
c20.SaveAs(f"plots/leadjetet_firstbunch.png")
c20.Close()



# ratio lead jet et
c2_ratio = ROOT.TCanvas( 'c2_ratio', 'c2_ratio', 100, 10, 600, 600 )

LeadJetEtbx0_u = inFile.Get("demo/LeadJetEtbx0_unprefirable")
LeadJetEtbxm1_u = inFile.Get("demo/LeadJetEtbxm1_unprefirable")

LeadJetEtcombined_u = LeadJetEtbx0_u.Clone()
LeadJetEtcombined_u.Add(LeadJetEtbxm1_u)

LeadJetEtratio_u = LeadJetEtbxm1_u.Clone()
LeadJetEtratio_u.Divide(LeadJetEtcombined_u)

LeadJetEtratio_u.SetLineColor(50)
LeadJetEtratio_u.SetTitle("")
LeadJetEtratio_u.GetXaxis().SetTitle("L1 Lead Jet E_{T} (GeV)")
LeadJetEtratio_u.GetYaxis().SetTitle("(bx=-1)/(bx=-1 or bx=0)")
LeadJetEtratio_u.Draw("")

LeadJetEtbx0_f = inFile.Get("demo/LeadJetEtbx0_firstbunch")
LeadJetEtbxm1_f = inFile.Get("demo/LeadJetEtbxm1_firstbunch")

LeadJetEtcombined_f = LeadJetEtbx0_f.Clone()
LeadJetEtcombined_f.Add(LeadJetEtbxm1_f)

LeadJetEtratio_f = LeadJetEtbxm1_f.Clone()
LeadJetEtratio_f.Divide(LeadJetEtcombined_f)

LeadJetEtratio_f.SetLineColor(9)
LeadJetEtratio_f.Draw("same")

legend = ROOT.TLegend(0.4,0.7,0.9,0.8)
legend.AddEntry(LeadJetEtratio_u,"UnprefirableEvent","L")
legend.AddEntry(LeadJetEtratio_f,"FirstBunchInTrain","L")
legend.SetBorderSize(0);
legend.SetTextSize(0.04);
legend.Draw()

c2_ratio.Draw()
c2_ratio.SaveAs(f"plots/leadjetet_ratio.png")
c2_ratio.Close()


#############################################################################
############################ lead jet eta plots #############################
#############################################################################

# unprefirable lead jet eta
c1 = ROOT.TCanvas( 'c1', 'c1', 100, 10, 600, 600 )
LeadJetEtabx0 = inFile.Get("demo/LeadJetEtabx0_unprefirable")
LeadJetEtabx0.SetTitle("Unprefirable")
LeadJetEtabx0.SetLineColor(1)
LeadJetEtabx0.GetXaxis().SetTitle("L1 Leading Jet #eta")
LeadJetEtabx0.GetYaxis().SetTitle("a.u.")
LeadJetEtabx0.SetStats(False)
LeadJetEtabx0.SetMaximum(1e5)
LeadJetEtabx0.Draw("")

LeadJetEtabxm1 = inFile.Get("demo/LeadJetEtabxm1_unprefirable")
LeadJetEtabxm1.SetStats(False)
LeadJetEtabxm1.SetLineColor(2)
LeadJetEtabxm1.Draw("same")

legend = ROOT.TLegend(0.6,0.72,0.9,0.87)
legend.AddEntry(LeadJetEtabx0,"BX=0","L")
legend.AddEntry(LeadJetEtabxm1,"BX=-1","L")
legend.SetBorderSize(0);
legend.SetTextSize(0.04);
legend.Draw()

c1.SetLogy();
c1.Draw()
c1.SaveAs(f"plots/leadjeteta_unprefirable.png")
c1.Close()

# firstbunchintrain lead jet eta
c100 = ROOT.TCanvas( 'c100', 'c100', 100, 10, 600, 600 )
LeadJetEtabx0 = inFile.Get("demo/LeadJetEtabx0_firstbunch")
LeadJetEtabx0.SetTitle("FirstBunchInTrain")
LeadJetEtabx0.SetLineColor(1)
LeadJetEtabx0.GetXaxis().SetTitle("L1 Leading Jet #eta")
LeadJetEtabx0.GetYaxis().SetTitle("a.u.")
LeadJetEtabx0.SetStats(False)
LeadJetEtabx0.SetMaximum(1e5)
LeadJetEtabx0.Draw("")

LeadJetEtabxm1 = inFile.Get("demo/LeadJetEtabxm1_firstbunch")
LeadJetEtabxm1.SetStats(False)
LeadJetEtabxm1.SetLineColor(2)
LeadJetEtabxm1.Draw("same")

legend = ROOT.TLegend(0.6,0.72,0.9,0.87)
legend.AddEntry(LeadJetEtabx0,"BX=0","L")
legend.AddEntry(LeadJetEtabxm1,"BX=-1","L")
legend.SetBorderSize(0);
legend.SetTextSize(0.04);
legend.Draw()

c100.SetLogy();
c100.Draw()
c100.SaveAs(f"plots/leadjeteta_firstbunch.png")
c100.Close()

# ratio lead jet eta
c2_ratio = ROOT.TCanvas( 'c2_ratio', 'c2_ratio', 100, 10, 600, 600 )

LeadJetEtabx0_u = inFile.Get("demo/LeadJetEtabx0_unprefirable")
LeadJetEtabxm1_u = inFile.Get("demo/LeadJetEtabxm1_unprefirable")

LeadJetEtacombined_u = LeadJetEtabx0_u.Clone()
LeadJetEtacombined_u.Add(LeadJetEtabxm1_u)

LeadJetEtaratio_u = LeadJetEtabxm1_u.Clone()
LeadJetEtaratio_u.Divide(LeadJetEtacombined_u)

LeadJetEtaratio_u.SetLineColor(50)
LeadJetEtaratio_u.SetTitle("")
LeadJetEtaratio_u.GetXaxis().SetTitle("L1 Lead Jet #eta")
LeadJetEtaratio_u.GetYaxis().SetTitle("(bx=-1)/(bx=-1 or bx=0)")
LeadJetEtaratio_u.Draw("")

LeadJetEtabx0_f = inFile.Get("demo/LeadJetEtabx0_firstbunch")
LeadJetEtabxm1_f = inFile.Get("demo/LeadJetEtabxm1_firstbunch")

LeadJetEtacombined_f = LeadJetEtabx0_f.Clone()
LeadJetEtacombined_f.Add(LeadJetEtabxm1_f)

LeadJetEtaratio_f = LeadJetEtabxm1_f.Clone()
LeadJetEtaratio_f.Divide(LeadJetEtacombined_f)

LeadJetEtaratio_f.SetLineColor(9)
LeadJetEtaratio_f.Draw("same")

legend = ROOT.TLegend(0.4,0.75,0.9,0.85)
legend.AddEntry(LeadJetEtaratio_u,"UnprefirableEvent","L")
legend.AddEntry(LeadJetEtaratio_f,"FirstBunchInTrain","L")
legend.SetBorderSize(0);
legend.SetTextSize(0.04);
legend.Draw()

c2_ratio.Draw()
c2_ratio.SaveAs(f"plots/leadjeteta_ratio.png")
c2_ratio.Close()


############################################################################
############################## njet HT plot ################################
############################################################################
c8 = ROOT.TCanvas( 'c8', 'c8', 100, 10, 600, 600 )
nJetHT_u = inFile.Get("demo/nJetHT_unprefirable")
nJetHT_u.SetMarkerColor(2)
nJetHT_u.SetLineColor(2)
nJetHT_u.SetStats(False)
nJetHT_u.SetMaximum(90000)
nJetHT_u.GetXaxis().SetTitle("Jets sum passing HT threshold")
nJetHT_u.GetYaxis().SetTitle("a.u.")
nJetHT_u.Draw("E0")
nJetHT_f = inFile.Get("demo/nJetHT_firstbunch")
nJetHT_f.SetStats(False)
nJetHT_f.SetMarkerColor(4)
nJetHT_f.SetLineColor(4)
nJetHT_f.Draw("E0 same")

legend = ROOT.TLegend(0.5,0.72,0.9,0.87)
legend.AddEntry(nJetHT_u,"Unprefirable","L")
legend.AddEntry(nJetHT_f,"FirstBunchInTrain","L")
legend.SetBorderSize(0);
legend.SetTextSize(0.04);
legend.Draw()

c8.Draw()
c8.SaveAs(f"plots/njetht_compare.png")
c8.Close()

############################################################################
################################ HT plots ##################################
############################################################################

# unprefirable HT
c4 = ROOT.TCanvas( 'c4', 'c4', 100, 10, 600, 600 )
HTbx0 = inFile.Get("demo/HTbx0_unprefirable")
HTbx0.SetTitle("Unprefirable")
HTbx0.SetLineColor(1)
HTbx0.GetXaxis().SetTitle("L1 H_{T} (GeV)")
HTbx0.GetYaxis().SetTitle("a.u.")
HTbx0.SetStats(False)
HTbx0.SetMaximum(1e6)
HTbx0.Draw("")

HTbxm1 = inFile.Get("demo/HTbxm1_unprefirable")
HTbxm1.SetStats(False)
HTbxm1.SetLineColor(2)
HTbxm1.Draw("same")

legend = ROOT.TLegend(0.6,0.6,0.9,0.8)
legend.AddEntry(HTbx0,"BX=0","L")
legend.AddEntry(HTbxm1,"BX=-1","L")
legend.SetBorderSize(0);
legend.SetTextSize(0.04);
legend.Draw()

c4.SetLogy();
c4.Draw()
c4.SaveAs(f"plots/HT_unprefirable.png")
c4.Close()

# firstbunch HT
c40 = ROOT.TCanvas( 'c40', 'c40', 100, 10, 600, 600 )
HTbx0 = inFile.Get("demo/HTbx0_firstbunch")
HTbx0.SetTitle("FirstBunchInTrain")
HTbx0.SetLineColor(1)
HTbx0.GetXaxis().SetTitle("L1 H_{T} (GeV)")
HTbx0.GetYaxis().SetTitle("a.u.")
HTbx0.SetStats(False)
HTbx0.SetMaximum(1e6)
HTbx0.Draw("")

HTbxm1 = inFile.Get("demo/HTbxm1_firstbunch")
HTbxm1.SetStats(False)
HTbxm1.SetLineColor(2)
HTbxm1.Draw("same")

legend2 = ROOT.TLegend(0.6,0.6,0.9,0.8)
legend2.AddEntry(HTbx0,"BX=0","L")
legend2.AddEntry(HTbxm1,"BX=-1","L")
legend2.SetBorderSize(0);
legend2.SetTextSize(0.04);
legend2.Draw()

c40.SetLogy();
c40.Draw()
c40.SaveAs(f"plots/HT_firstbunch.png")
c40.Close()

# ratio ht
c4_ratio = ROOT.TCanvas( 'c4_ratio', 'c4_ratio', 100, 10, 600, 600 )

HTbx0_u = inFile.Get("demo/HTbx0_unprefirable")
HTbxm1_u = inFile.Get("demo/HTbxm1_unprefirable")

HTcombined_u = HTbx0_u.Clone()
HTcombined_u.Add(HTbxm1_u)

HTratio_u = HTbxm1_u.Clone()
HTratio_u.Divide(HTcombined_u)

HTratio_u.SetLineColor(50)
HTratio_u.SetTitle("")
HTratio_u.GetXaxis().SetTitle("Event H_{T}")
HTratio_u.GetYaxis().SetTitle("(bx=-1)/(bx=-1 or bx=0)")
HTratio_u.Draw("")

HTbx0_f = inFile.Get("demo/HTbx0_firstbunch")
HTbxm1_f = inFile.Get("demo/HTbxm1_firstbunch")

HTcombined_f = HTbx0_f.Clone()
HTcombined_f.Add(HTbxm1_f)

HTratio_f = HTbxm1_f.Clone()
HTratio_f.Divide(HTcombined_f)

HTratio_f.SetLineColor(9)
HTratio_f.Draw("same")

legend = ROOT.TLegend(0.4,0.75,0.9,0.85)
legend.AddEntry(HTratio_u,"UnprefirableEvent","L")
legend.AddEntry(HTratio_f,"FirstBunchInTrain","L")
legend.SetBorderSize(0);
legend.SetTextSize(0.04);
legend.Draw()

c4_ratio.Draw()
c4_ratio.SaveAs(f"plots/HT_ratio.png")
c4_ratio.Close()


############################################################################
############################# eta/phi plots ################################
############################################################################
c23 = ROOT.TCanvas( 'c23', 'c23', 100, 10, 600, 600 )
EtaPhi = inFile.Get("demo/EtaPhi_unprefirable")
EtaPhi.GetXaxis().SetTitle("#eta")
EtaPhi.GetYaxis().SetTitle("#phi")
EtaPhi.SetTitle("(Unprefirable) Jets Passing Pre-trigger Threshold")
EtaPhi.SetStats(False)
EtaPhi.Draw("")
c23.Draw()
c23.SaveAs(f"plots/etaphi_unprefirable.png")
c23.Close()

c24 = ROOT.TCanvas( 'c24', 'c24', 100, 10, 600, 600 )
EtaPhi = inFile.Get("demo/EtaPhi_firstbunch")
EtaPhi.GetXaxis().SetTitle("#eta")
EtaPhi.GetYaxis().SetTitle("#phi")
EtaPhi.SetTitle("(FirstBunchInTrain) Jets Passing Pre-trigger Threshold")
EtaPhi.SetStats(False)
EtaPhi.Draw("")
c24.Draw()
c24.SaveAs(f"plots/etaphi_firstbunch.png")
c24.Close()
