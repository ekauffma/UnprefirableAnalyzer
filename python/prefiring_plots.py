import argparse
import awkward as ak
import numpy as np
import ROOT
import uproot


def main(file_path, out_dir):
    
    inFile = ROOT.TFile.Open(file_path, "READ")
    
    ###################################################################################
    ################################ n jet per bx plot ################################
    ###################################################################################
    
    c = ROOT.TCanvas( 'c', 'c', 100, 10, 600, 600 )
    nJets_bx = inFile.Get("demo/nJets_bx")
    nJets_bx.SetLineColor(2)
    nJets_bx.SetMarkerColor(2)
    nJets_bx.SetMarkerStyle(47)
    nJets_bx.GetXaxis().SetTitle("BX")
    nJets_bx.GetYaxis().SetTitle("Number of Jets")
    nJets_bx.SetStats(False)
    nJets_bx.SetMaximum(1e8)
    nJets_bx.SetMinimum(1e1)
    nJets_bx.Draw("E")

    c.SetLogy()
    c.Draw()
    c.SaveAs(f"{out_dir}/njet_bx.png")
    c.Close()
    
    
    ##############################################################################
    ################################ jet pt plots ################################
    ##############################################################################
    
    c = ROOT.TCanvas( 'c', 'c', 100, 10, 600, 600 )
    JetEtbx0 = inFile.Get("demo/JetEt_bx0")
    JetEtbx0.SetTitle("Unprefirable OR FirstBunchInTrain")
    JetEtbx0.SetLineColor(2)
    JetEtbx0.SetMarkerColor(2)
    JetEtbx0.SetMarkerStyle(47)
    JetEtbx0.GetXaxis().SetTitle("Reco Jet E_{T} (GeV)")
    JetEtbx0.GetYaxis().SetTitle("a.u.")
    JetEtbx0.SetStats(False)
    JetEtbx0.SetMaximum(1e7)
    JetEtbx0.SetMinimum(1e-1)
    JetEtbx0.Draw("E")

    JetEtbxm1 = inFile.Get("demo/JetEt_bxm1")
    JetEtbxm1.SetStats(False)
    JetEtbxm1.SetLineColor(55)
    JetEtbxm1.SetMarkerColor(55)
    JetEtbxm1.SetMarkerStyle(22)
    JetEtbxm1.Draw("E same")

    JetEt_bx0_bxm1 = inFile.Get("demo/JetEt_bx0_bxm1")
    JetEt_bx0_bxm1.SetStats(False)
    JetEt_bx0_bxm1.SetLineColor(8)
    JetEt_bx0_bxm1.SetMarkerColor(8)
    JetEt_bx0_bxm1.SetMarkerStyle(33)
    JetEt_bx0_bxm1.Draw("E same")

    JetEt_bxm2 = inFile.Get("demo/JetEt_bxm2")
    JetEt_bxm2.SetStats(False)
    JetEt_bxm2.SetLineColor(6)
    JetEt_bxm2.SetMarkerColor(6)
    JetEt_bxm2.SetMarkerStyle(24)
    JetEt_bxm2.Draw("E same")


    legend = ROOT.TLegend(0.5,0.7,0.9,0.85)
    legend.AddEntry(JetEtbx0,"BX=0","PE")
    legend.AddEntry(JetEtbxm1,"BX=-1","PE")
    legend.AddEntry(JetEt_bx0_bxm1,"BX=0 or BX=-1","PE")
    legend.AddEntry(JetEt_bxm2,"BX=-2","PE")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c.SetGrid()
    c.SetLogy()
    c.Draw()
    c.SaveAs(f"{out_dir}/jetpt.png")
    c.Close()
    
    # ratio
    c_et_rat = ROOT.TCanvas( 'c_et_rat', 'c_et_rat', 10, 10, 800, 600 )

    JetEt_bxm1 = inFile.Get("demo/JetEt_bxm1")
    JetEt_bx0_bxm1 = inFile.Get("demo/JetEt_bx0_bxm1")

    JetEt_ratio = JetEt_bxm1.Clone()
    JetEt_ratio.Sumw2()
    JetEt_bx0_bxm1.Sumw2()
    JetEt_ratio.Divide(JetEt_bx0_bxm1)

    JetEt_ratio.SetLineColor(50)
    JetEt_ratio.SetMarkerColor(50)
    JetEt_ratio.SetMarkerStyle(8)
    JetEt_ratio.SetMarkerSize(0.5)
    JetEt_ratio.SetTitle("Unprefirable OR FirstBunchInTrain")
    JetEt_ratio.GetXaxis().SetTitle("Reco Jet E_{T} (GeV)")
    JetEt_ratio.GetYaxis().SetTitle("(bx=-1)/(bx=-1 or bx=0)")
    JetEt_ratio.SetMinimum(0)
    JetEt_ratio.SetMaximum(0.1)
    JetEt_ratio.Draw("E")

    c_et_rat.Draw()
    c_et_rat.SaveAs(f"{out_dir}/jetpt_ratio.png")
    c_et_rat.Close()
   
    # low resolution plot
    c = ROOT.TCanvas( 'c', 'c', 100, 10, 600, 600 )
    JetEtbx0 = inFile.Get("demo/JetPt_lowres_bx0")
    JetEtbx0.SetTitle("Jet E_{T} Distribution for Low Resolution (<0.5) Jets - Unprefirable OR FirstBunchInTrain")
    JetEtbx0.SetLineColor(2)
    JetEtbx0.SetMarkerColor(2)
    JetEtbx0.SetMarkerStyle(47)
    JetEtbx0.GetXaxis().SetTitle("Reco Jet E_{T} (GeV)")
    JetEtbx0.GetYaxis().SetTitle("a.u.")
    JetEtbx0.SetStats(False)
    JetEtbx0.SetMaximum(1e6)
    JetEtbx0.SetMinimum(1e-1)
    JetEtbx0.Draw("E")

    JetEtbxm1 = inFile.Get("demo/JetPt_lowres_bxm1")
    JetEtbxm1.SetStats(False)
    JetEtbxm1.SetLineColor(55)
    JetEtbxm1.SetMarkerColor(55)
    JetEtbxm1.SetMarkerStyle(22)
    JetEtbxm1.Draw("E same")

    JetEt_bxm2 = inFile.Get("demo/JetPt_lowres_bxm2")
    JetEt_bxm2.SetStats(False)
    JetEt_bxm2.SetLineColor(6)
    JetEt_bxm2.SetMarkerColor(6)
    JetEt_bxm2.SetMarkerStyle(24)
    JetEt_bxm2.Draw("E same")


    legend = ROOT.TLegend(0.5,0.7,0.9,0.85)
    legend.AddEntry(JetEtbx0,"BX=0","PE")
    legend.AddEntry(JetEtbxm1,"BX=-1","PE")
    legend.AddEntry(JetEt_bxm2,"BX=-2","PE")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c.SetGrid()
    c.SetLogy()
    c.Draw()
    c.SaveAs(f"{out_dir}/jetpt_lowres.png")
    c.Close()

    
    #############################################################################
    ############################### jet eta plots ###############################
    #############################################################################
    
    c1 = ROOT.TCanvas( 'c1', 'c1', 100, 10, 600, 600 )
    JetEtabx0 = inFile.Get("demo/JetEta_bx0")
    JetEtabx0.SetTitle("Unprefirable OR FirstBunchInTrain")
    JetEtabx0.SetLineColor(2)
    JetEtabx0.SetMarkerColor(2)
    JetEtabx0.SetMarkerStyle(47)
    JetEtabx0.GetXaxis().SetTitle("Reco Jet #eta")
    JetEtabx0.GetYaxis().SetTitle("a.u.")
    JetEtabx0.SetStats(False)
    JetEtabx0.SetMaximum(1e7)
    JetEtabx0.SetMinimum(1e-1)
    JetEtabx0.Draw("E")

    JetEtabxm1 = inFile.Get("demo/JetEta_bxm1")
    JetEtabxm1.SetStats(False)
    JetEtabxm1.SetLineColor(55)
    JetEtabxm1.SetMarkerColor(55)
    JetEtabxm1.SetMarkerStyle(22)
    JetEtabxm1.Draw("E same")

    JetEta_bx0_bxm1 = inFile.Get("demo/JetEta_bx0_bxm1")
    JetEta_bx0_bxm1.SetStats(False)
    JetEta_bx0_bxm1.SetLineColor(8)
    JetEta_bx0_bxm1.SetMarkerColor(8)
    JetEta_bx0_bxm1.SetMarkerStyle(33)
    JetEta_bx0_bxm1.Draw("E same")

    JetEta_bxm2 = inFile.Get("demo/JetEta_bxm2")
    JetEta_bxm2.SetStats(False)
    JetEta_bxm2.SetLineColor(6)
    JetEta_bxm2.SetMarkerColor(6)
    JetEta_bxm2.SetMarkerStyle(24)
    JetEta_bxm2.Draw("E same")

    legend = ROOT.TLegend(0.35,0.75,0.7,0.89)
    legend.AddEntry(JetEtabx0,"BX=0","PE")
    legend.AddEntry(JetEtabxm1,"BX=-1","PE")
    legend.AddEntry(JetEta_bx0_bxm1,"BX=0 or BX=-1","PE")
    legend.AddEntry(JetEta_bxm2,"BX=-2","PE")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c1.SetGrid()
    c1.SetLogy()
    c1.Draw()
    c1.SaveAs(f"{out_dir}/jeteta.png")
    c1.Close()

    # RATIO
    c_eta_rat = ROOT.TCanvas( 'c_eta_rat', 'c_eta_rat', 10, 10, 800, 600 )

    JetEta_bxm1 = inFile.Get("demo/JetEta_bxm1")
    JetEta_bx0_bxm1 = inFile.Get("demo/JetEta_bx0_bxm1")

    JetEta_ratio = JetEta_bxm1.Clone()
    
    JetEta_ratio.Sumw2()
    JetEta_bx0_bxm1.Sumw2()
    JetEta_ratio.Divide(JetEta_bx0_bxm1)

    JetEta_ratio.SetLineColor(50)
    JetEta_ratio.SetMarkerColor(50)
    JetEta_ratio.SetMarkerStyle(8)
    JetEta_ratio.SetMarkerSize(0.5)
    JetEta_ratio.SetTitle("Unprefirable OR FirstBunchInTrain")
    JetEta_ratio.GetXaxis().SetTitle("Reco Jet #eta")
    JetEta_ratio.GetYaxis().SetTitle("(bx=-1)/(bx=-1 or bx=0)")
    JetEta_ratio.SetMinimum(0)
    JetEta_ratio.SetMaximum(0.03)
    JetEta_ratio.Draw("E")

    c_eta_rat.Draw()
    c_eta_rat.SaveAs(f"{out_dir}/jeteta_ratio.png")
    c_eta_rat.Close()
   
    # RATIO SEPARATED PT
    c_eta_rat = ROOT.TCanvas( 'c_eta_rat', 'c_eta_rat', 10, 10, 800, 600 )

    JetEta_low_bxm1 = inFile.Get("demo/JetEta_lowpt_bxm1")
    JetEta_low_bx0_bxm1 = inFile.Get("demo/JetEta_lowpt_bx0_bxm1")

    JetEta_low_ratio = JetEta_low_bxm1.Clone()

    JetEta_low_ratio.Sumw2()
    JetEta_low_bx0_bxm1.Sumw2()
    JetEta_low_ratio.Divide(JetEta_low_bx0_bxm1)

    JetEta_low_ratio.SetLineColor(50)
    JetEta_low_ratio.SetMarkerColor(50)
    JetEta_low_ratio.SetMarkerStyle(8)
    JetEta_low_ratio.SetMarkerSize(0.5)
    JetEta_low_ratio.SetTitle("Unprefirable OR FirstBunchInTrain")
    JetEta_low_ratio.GetXaxis().SetTitle("Reco Jet #eta")
    JetEta_low_ratio.GetYaxis().SetTitle("(bx=-1)/(bx=-1 or bx=0)")
    JetEta_low_ratio.SetMinimum(0)
    JetEta_low_ratio.SetMaximum(0.1)
    JetEta_low_ratio.SetStats(False)
    JetEta_low_ratio.Draw("E")

    JetEta_med_bxm1 = inFile.Get("demo/JetEta_medpt_bxm1")
    JetEta_med_bx0_bxm1 = inFile.Get("demo/JetEta_medpt_bx0_bxm1")

    JetEta_med_ratio = JetEta_med_bxm1.Clone()

    JetEta_med_ratio.Sumw2()
    JetEta_med_bx0_bxm1.Sumw2()
    JetEta_med_ratio.Divide(JetEta_med_bx0_bxm1)

    JetEta_med_ratio.SetLineColor(54)
    JetEta_med_ratio.SetMarkerColor(54)
    JetEta_med_ratio.SetMarkerStyle(29)
    JetEta_med_ratio.SetMarkerSize(0.5)
    JetEta_med_ratio.SetStats(False);
    JetEta_med_ratio.Draw("E same")

    JetEta_high_bxm1 = inFile.Get("demo/JetEta_high_bxm1")
    JetEta_high_bx0_bxm1 = inFile.Get("demo/JetEta_high_bx0_bxm1")

    JetEta_high_ratio = JetEta_high_bxm1.Clone()

    JetEta_high_ratio.Sumw2()
    JetEta_high_bx0_bxm1.Sumw2()
    JetEta_high_ratio.Divide(JetEta_high_bx0_bxm1)

    JetEta_high_ratio.SetLineColor(108)
    JetEta_high_ratio.SetMarkerColor(108)
    JetEta_high_ratio.SetMarkerStyle(21)
    JetEta_high_ratio.SetMarkerSize(0.5)
    JetEta_high_ratio.SetStats(False);
    JetEta_high_ratio.Draw("E same")

    legend = ROOT.TLegend(0.3,0.75,0.9,0.88)
    legend.AddEntry(JetEta_low_ratio,"L1 p_T <= 15 GeV","L")
    legend.AddEntry(JetEta_med_ratio,"15 GeV < L1 p_T <= 30 GeV","L")
    legend.AddEntry(JetEta_high_ratio,"L1 p_T > 30 GeV","L")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c_eta_rat.Draw()
    c_eta_rat.SaveAs(f"{out_dir}/jeteta_ratio_ptsep.png")
    c_eta_rat.Close() 
   
    ############################################################################
    ############################# eta/phi plots ################################
    ############################################################################
    
    c_ep_0 = ROOT.TCanvas( 'c_ep_0', 'c_ep_0', 100, 10, 600, 600 )
    EtaPhi = inFile.Get("demo/JetEtaPhi_bx0")
    EtaPhi.GetXaxis().SetTitle("#eta")
    EtaPhi.GetYaxis().SetTitle("#phi")
    EtaPhi.SetTitle("Unprefirable OR FirstBunchInTrain (Offline) Jets BX=0 (Offline Jet pT>30 GeV)")
    EtaPhi.SetStats(False)
    EtaPhi.Draw("colz")
    c_ep_0.Draw()
    c_ep_0.SaveAs(f"{out_dir}/jetetaphi_bx0_offline.png")
    c_ep_0.Close()

    c_ep_m1 = ROOT.TCanvas( 'c_ep_m1', 'c_ep_m1', 100, 10, 600, 600 )
    EtaPhi = inFile.Get("demo/JetEtaPhi_bxm1")
    EtaPhi.GetXaxis().SetTitle("#eta")
    EtaPhi.GetYaxis().SetTitle("#phi")
    EtaPhi.SetTitle("Unprefirable OR FirstBunchInTrain (Offline) Jets BX=-1 (Offline Jet pT>30 GeV)")
    EtaPhi.SetStats(False)
    EtaPhi.Draw("colz")
    c_ep_m1.Draw()
    c_ep_m1.SaveAs(f"{out_dir}/jetetaphi_bxm1_offline.png")
    c_ep_m1.Close()

    c_ep_0_on = ROOT.TCanvas( 'c_ep_0_on', 'c_ep_0_on', 100, 10, 600, 600 )
    EtaPhi = inFile.Get("demo/JetEtaPhi_bx0_online")
    EtaPhi.GetXaxis().SetTitle("#eta")
    EtaPhi.GetYaxis().SetTitle("#phi")
    EtaPhi.SetTitle("Unprefirable OR FirstBunchInTrain (Online) Jets BX=0 (Offline Jet pT>30 GeV)")
    EtaPhi.SetStats(False)
    EtaPhi.Draw("colz")
    c_ep_0_on.Draw()
    c_ep_0_on.SaveAs(f"{out_dir}/jetetaphi_bx0_online.png")
    c_ep_0_on.Close()

    c_ep_m1_on = ROOT.TCanvas( 'c_ep_m1_on', 'c_ep_m1_on', 100, 10, 600, 600 )
    EtaPhi = inFile.Get("demo/JetEtaPhi_bxm1_online")
    EtaPhi.GetXaxis().SetTitle("#eta")
    EtaPhi.GetYaxis().SetTitle("#phi")
    EtaPhi.SetTitle("Unprefirable OR FirstBunchInTrain (Online) Jets BX=-1 (Offline Jet pT>30 GeV)")
    EtaPhi.SetStats(False)
    EtaPhi.Draw("colz")
    c_ep_m1_on.Draw()
    c_ep_m1_on.SaveAs(f"{out_dir}/jetetaphi_bxm1_online.png")
    c_ep_m1_on.Close()
    
    ############################################################################
    ############################# pt/eta plots ################################
    ############################################################################

    c_pe_0 = ROOT.TCanvas( 'c_pe_0', 'c_pe_0', 100, 10, 600, 600 )
    PtEta = inFile.Get("demo/JetPtEta_bx0")
    PtEta.GetXaxis().SetTitle("p_T (GeV)")
    PtEta.GetYaxis().SetTitle("#eta")
    PtEta.SetTitle("Unprefirable OR FirstBunchInTrain (Offline) Jets BX=0 (Offline Jet pT>30 GeV)")
    PtEta.SetStats(False)
    PtEta.Draw("colz")
    c_pe_0.Draw()
    c_pe_0.SaveAs(f"{out_dir}/jetpteta_bx0_offline.png")
    c_pe_0.Close()
 
    c_pe_m1 = ROOT.TCanvas( 'c_pe_m1', 'c_pe_m1', 100, 10, 600, 600 )
    PtEta = inFile.Get("demo/JetPtEta_bxm1")
    PtEta.GetXaxis().SetTitle("p_T (GeV)")
    PtEta.GetYaxis().SetTitle("#eta")
    PtEta.SetTitle("Unprefirable OR FirstBunchInTrain (Offline) Jets BX=-1 (Offline Jet pT>30 GeV)")
    PtEta.SetStats(False)
    PtEta.Draw("colz")
    c_pe_m1.Draw()
    c_pe_m1.SaveAs(f"{out_dir}/jetpteta_bxm1_offline.png")
    c_pe_m1.Close()

    ################################################################################
    ################################# pt res plots #################################
    ################################################################################
    
    c1 = ROOT.TCanvas( 'c1', 'c1', 100, 10, 1000, 600 )
    JetEResbx0 = inFile.Get("demo/JetERes_bx0")
    JetEResbx0.SetTitle("Unprefirable OR FirstBunchInTrain")
    JetEResbx0.SetLineColor(2)
    JetEResbx0.SetMarkerColor(2)
    JetEResbx0.SetMarkerStyle(47)
    JetEResbx0.GetXaxis().SetTitle("Online Jet p_{T} / Offline Jet p_{T}")
    JetEResbx0.GetYaxis().SetTitle("a.u.")
    JetEResbx0.SetStats(False)
    JetEResbx0.SetMaximum(1e7)
    JetEResbx0.Draw("E")

    JetEResbxm1 = inFile.Get("demo/JetERes_bxm1")
    JetEResbxm1.SetStats(False)
    JetEResbxm1.SetLineColor(55)
    JetEResbxm1.SetMarkerColor(55)
    JetEResbxm1.SetMarkerStyle(22)
    JetEResbxm1.Draw("E same")

    JetEResbxm2 = inFile.Get("demo/JetERes_bxm2")
    JetEResbxm2.SetStats(False)
    JetEResbxm2.SetLineColor(6)
    JetEResbxm2.SetMarkerColor(6)
    JetEResbxm2.SetMarkerStyle(23)
    JetEResbxm2.Draw("E same")

    JetEResbx1 = inFile.Get("demo/JetERes_bx1")
    JetEResbx1.SetStats(False)
    JetEResbx1.SetLineColor(8)
    JetEResbx1.SetMarkerColor(8)
    JetEResbx1.SetMarkerStyle(42)
    JetEResbx1.Draw("E same")

    JetEResbx2 = inFile.Get("demo/JetERes_bx2")
    JetEResbx2.SetStats(False)
    JetEResbx2.SetLineColor(95)
    JetEResbx2.SetMarkerColor(95)
    JetEResbx2.SetMarkerStyle(43)
    JetEResbx2.Draw("E same")

    legend = ROOT.TLegend(0.7,0.7,0.9,0.85)
    legend.AddEntry(JetEResbx0,"BX=0","PE")
    legend.AddEntry(JetEResbxm1,"BX=-1","PE")
    legend.AddEntry(JetEResbxm2,"BX=-2","PE")
    legend.AddEntry(JetEResbx1,"BX=1","PE")
    legend.AddEntry(JetEResbx2,"BX=2","PE")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c1.SetLogy()
    c1.Draw()
    c1.SaveAs(f"{out_dir}/jetptres.png")
    c1.Close()

    
if __name__ == "__main__":
    
    # parse arguments
    parser = argparse.ArgumentParser(description="This program creates plots to study jet pre-firing in the L1 Trigger.")
    parser.add_argument("-p", "--file_path", help="path to input ROOT file containing histograms")
    parser.add_argument("-o", "--out_dir", help="path to directory to save plots", default=".")
    
    args = parser.parse_args()
    
    main(args.file_path, args.out_dir)
