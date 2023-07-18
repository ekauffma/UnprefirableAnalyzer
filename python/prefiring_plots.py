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
    nJets_bx_unprefirable = inFile.Get("demo/nJets_bx_unprefirable")
    nJets_bx_unprefirable.SetLineColor(2)
    nJets_bx_unprefirable.SetMarkerColor(2)
    nJets_bx_unprefirable.SetMarkerStyle(47)
    nJets_bx_unprefirable.GetXaxis().SetTitle("BX")
    nJets_bx_unprefirable.GetYaxis().SetTitle("Number of Jets")
    nJets_bx_unprefirable.SetStats(False)
    nJets_bx_unprefirable.SetMaximum(2000)
    nJets_bx_unprefirable.Draw("E")

    nJets_bx_firstbunch = inFile.Get("demo/nJets_bx_firstbunch")
    nJets_bx_firstbunch.SetLineColor(8)
    nJets_bx_firstbunch.SetMarkerColor(8)
    nJets_bx_firstbunch.SetMarkerStyle(40)
    nJets_bx_firstbunch.SetStats(False)
    nJets_bx_firstbunch.Draw("E same")

    legend = ROOT.TLegend(0.5,0.8,0.9,0.88)
    legend.AddEntry(nJets_bx_unprefirable,"Unprefirable","PE")
    legend.AddEntry(nJets_bx_firstbunch,"FirstBunchInTrain","PE")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c.Draw()
    c.SaveAs(f"{out_dir}/njet_bx.png")
    c.Close()
    
    
    ##############################################################################
    ################################ jet pt plots ################################
    ##############################################################################
    
    # UNPREFIRABLE
    c = ROOT.TCanvas( 'c', 'c', 100, 10, 600, 600 )
    JetEtbx0_u = inFile.Get("demo/JetEt_bx0_unprefirable")
    JetEtbx0_u.SetTitle("Unprefirable")
    JetEtbx0_u.SetLineColor(2)
    JetEtbx0_u.SetMarkerColor(2)
    JetEtbx0_u.SetMarkerStyle(47)
    JetEtbx0_u.GetXaxis().SetTitle("Reco Jet E_{T} (GeV)")
    JetEtbx0_u.GetYaxis().SetTitle("a.u.")
    JetEtbx0_u.SetStats(False)
    JetEtbx0_u.SetMaximum(1e3)
    JetEtbx0_u.SetMinimum(1e-1)
    JetEtbx0_u.Draw("E")

    JetEtbxm1_u = inFile.Get("demo/JetEt_bxm1_unprefirable")
    JetEtbxm1_u.SetStats(False)
    JetEtbxm1_u.SetLineColor(55)
    JetEtbxm1_u.SetMarkerColor(55)
    JetEtbxm1_u.SetMarkerStyle(22)
    JetEtbxm1_u.Draw("E same")

    JetEt_bx0_bxm1_u = inFile.Get("demo/JetEt_bx0_bxm1_unprefirable")
    JetEt_bx0_bxm1_u.SetStats(False)
    JetEt_bx0_bxm1_u.SetLineColor(8)
    JetEt_bx0_bxm1_u.SetMarkerColor(8)
    JetEt_bx0_bxm1_u.SetMarkerStyle(33)
    JetEt_bx0_bxm1_u.Draw("E same")

    JetEt_bxm2_u = inFile.Get("demo/JetEt_bxm2_unprefirable")
    JetEt_bxm2_u.SetStats(False)
    JetEt_bxm2_u.SetLineColor(6)
    JetEt_bxm2_u.SetMarkerColor(6)
    JetEt_bxm2_u.SetMarkerStyle(24)
    JetEt_bxm2_u.Draw("E same")


    legend = ROOT.TLegend(0.5,0.7,0.9,0.85)
    legend.AddEntry(JetEtbx0_u,"BX=0","PE")
    legend.AddEntry(JetEtbxm1_u,"BX=-1","PE")
    legend.AddEntry(JetEt_bx0_bxm1_u,"BX=0 or BX=-1","PE")
    legend.AddEntry(JetEt_bxm2_u,"BX=-2","PE")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c.SetGrid()
    c.SetLogy()
    c.Draw()
    c.SaveAs(f"{out_dir}/jetpt_unprefirable.png")
    c.Close()
    
    # FIRSTBUNCHINTRAIN
    c = ROOT.TCanvas( 'c', 'c', 100, 10, 600, 600 )
    JetEtbx0_f = inFile.Get("demo/JetEt_bx0_firstbunch")
    JetEtbx0_f.SetTitle("FirstBunchInTrain")
    JetEtbx0_f.GetXaxis().SetTitle("Reco Jet E_{T} (GeV)")
    JetEtbx0_f.GetYaxis().SetTitle("a.u.")
    JetEtbx0_f.SetStats(False)
    JetEtbx0_f.SetLineColor(2)
    JetEtbx0_f.SetMarkerColor(2)
    JetEtbx0_f.SetMarkerStyle(47)
    JetEtbx0_f.SetMaximum(1e4)
    JetEtbx0_f.SetMinimum(1e-1)
    JetEtbx0_f.Draw("E")

    JetEtbxm1_f = inFile.Get("demo/JetEt_bxm1_firstbunch")
    JetEtbxm1_f.SetStats(False)
    JetEtbxm1_f.SetLineColor(55)
    JetEtbxm1_f.SetMarkerColor(55)
    JetEtbxm1_f.SetMarkerStyle(22)
    JetEtbxm1_f.Draw("E same")

    JetEt_bx0_bxm1_f = inFile.Get("demo/JetEt_bx0_bxm1_firstbunch")
    JetEt_bx0_bxm1_f.SetStats(False)
    JetEt_bx0_bxm1_f.SetLineColor(8)
    JetEt_bx0_bxm1_f.SetMarkerColor(8)
    JetEt_bx0_bxm1_f.SetMarkerStyle(33)
    JetEt_bx0_bxm1_f.Draw("E same")

    JetEt_bxm2_f = inFile.Get("demo/JetEt_bxm2_firstbunch")
    JetEt_bxm2_f.SetStats(False)
    JetEt_bxm2_f.SetLineColor(6)
    JetEt_bxm2_f.SetMarkerColor(6)
    JetEt_bxm2_f.SetMarkerStyle(24)
    JetEt_bxm2_f.Draw("E same")

    legend = ROOT.TLegend(0.5,0.7,0.9,0.85)
    legend.AddEntry(JetEtbx0_f,"BX=0","PE")
    legend.AddEntry(JetEtbxm1_f,"BX=-1","PE")
    legend.AddEntry(JetEt_bx0_bxm1_f,"BX=0 or BX=-1","PE")
    legend.AddEntry(JetEt_bxm2_f,"BX=-2","PE")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c.SetGrid()
    c.SetLogy()
    c.Draw()
    c.SaveAs(f"{out_dir}/jetpt_firstbunch.png")
    c.Close()
    
    # ratio
    c_et_rat = ROOT.TCanvas( 'c_et_rat', 'c_et_rat', 100, 10, 600, 600 )

    JetEt_bxm1_u = inFile.Get("demo/JetEt_bxm1_unprefirable")
    JetEt_bx0_bxm1_u = inFile.Get("demo/JetEt_bx0_bxm1_unprefirable")

    JetEt_ratio_u = JetEt_bxm1_u.Clone()
    JetEt_ratio_u.Divide(JetEt_bx0_bxm1_u)

    JetEt_ratio_u.SetLineColor(50)
    JetEt_ratio_u.SetMarkerColor(50)
    JetEt_ratio_u.SetMarkerStyle(8)
    JetEt_ratio_u.SetMarkerSize(0.5)
    JetEt_ratio_u.SetTitle("")
    JetEt_ratio_u.GetXaxis().SetTitle("Reco Jet E_{T} (GeV)")
    JetEt_ratio_u.GetYaxis().SetTitle("(bx=-1)/(bx=-1 or bx=0)")
    JetEt_ratio_u.SetMinimum(0)
    JetEt_ratio_u.SetMaximum(0.6)
    JetEt_ratio_u.Draw("E")

    JetEt_bxm1_f = inFile.Get("demo/JetEt_bxm1_firstbunch")
    JetEt_bx0_bxm1_f = inFile.Get("demo/JetEt_bx0_bxm1_firstbunch")

    JetEt_ratio_f = JetEt_bxm1_f.Clone()
    JetEt_ratio_f.Divide(JetEt_bx0_bxm1_f)

    JetEt_ratio_f.SetLineColor(9)
    JetEt_ratio_f.SetMarkerColor(9)
    JetEt_ratio_f.SetMarkerStyle(8)
    JetEt_ratio_f.SetMarkerSize(0.5)
    JetEt_ratio_f.Draw("E same")

    legend = ROOT.TLegend(0.5,0.8,0.9,0.88)
    legend.AddEntry(JetEt_ratio_u,"UnprefirableEvent","L")
    legend.AddEntry(JetEt_ratio_f,"FirstBunchInTrain","L")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c_et_rat.Draw()
    c_et_rat.SaveAs(f"{out_dir}/jetpt_ratio.png")
    c_et_rat.Close()
    
    
    #############################################################################
    ############################### jet eta plots ###############################
    #############################################################################
    
    # UNPREFIRABLE
    c1 = ROOT.TCanvas( 'c1', 'c1', 100, 10, 600, 600 )
    JetEtabx0_u = inFile.Get("demo/JetEta_bx0_unprefirable")
    JetEtabx0_u.SetTitle("Unprefirable")
    JetEtabx0_u.SetLineColor(2)
    JetEtabx0_u.SetMarkerColor(2)
    JetEtabx0_u.SetMarkerStyle(47)
    JetEtabx0_u.GetXaxis().SetTitle("Reco Jet #eta (GeV)")
    JetEtabx0_u.GetYaxis().SetTitle("a.u.")
    JetEtabx0_u.SetStats(False)
    JetEtabx0_u.SetMaximum(1e3)
    JetEtabx0_u.SetMinimum(1e-1)
    JetEtabx0_u.Draw("E")

    JetEtabxm1_u = inFile.Get("demo/JetEta_bxm1_unprefirable")
    JetEtabxm1_u.SetStats(False)
    JetEtabxm1_u.SetLineColor(55)
    JetEtabxm1_u.SetMarkerColor(55)
    JetEtabxm1_u.SetMarkerStyle(22)
    JetEtabxm1_u.Draw("E same")

    JetEta_bx0_bxm1_u = inFile.Get("demo/JetEta_bx0_bxm1_unprefirable")
    JetEta_bx0_bxm1_u.SetStats(False)
    JetEta_bx0_bxm1_u.SetLineColor(8)
    JetEta_bx0_bxm1_u.SetMarkerColor(8)
    JetEta_bx0_bxm1_u.SetMarkerStyle(33)
    JetEta_bx0_bxm1_u.Draw("E same")

    JetEta_bxm2_u = inFile.Get("demo/JetEta_bxm2_unprefirable")
    JetEta_bxm2_u.SetStats(False)
    JetEta_bxm2_u.SetLineColor(6)
    JetEta_bxm2_u.SetMarkerColor(6)
    JetEta_bxm2_u.SetMarkerStyle(24)
    JetEta_bxm2_u.Draw("E same")

    legend = ROOT.TLegend(0.35,0.7,0.7,0.88)
    legend.AddEntry(JetEtabx0_u,"BX=0","PE")
    legend.AddEntry(JetEtabxm1_u,"BX=-1","PE")
    legend.AddEntry(JetEta_bx0_bxm1_u,"BX=0 or BX=-1","PE")
    legend.AddEntry(JetEta_bxm2_u,"BX=-2","PE")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c1.SetGrid()
    c1.SetLogy()
    c1.Draw()
    c1.SaveAs(f"{out_dir}/jeteta_unprefirable.png")
    c1.Close()

    # FIRSTBUNCHINTRAIN
    c = ROOT.TCanvas( 'c', 'c', 100, 10, 600, 600 )
    JetEtabx0_f = inFile.Get("demo/JetEta_bx0_firstbunch")
    JetEtabx0_f.SetTitle("FirstBunchInTrain")
    JetEtabx0_f.GetXaxis().SetTitle("Reco Jet #eta (GeV)")
    JetEtabx0_f.GetYaxis().SetTitle("a.u.")
    JetEtabx0_f.SetStats(False)
    JetEtabx0_f.SetLineColor(2)
    JetEtabx0_f.SetMarkerColor(2)
    JetEtabx0_f.SetMarkerStyle(47)
    JetEtabx0_f.SetMaximum(5e4)
    JetEtabx0_f.SetMinimum(1e-1)
    JetEtabx0_f.Draw("E")

    JetEtabxm1_f = inFile.Get("demo/JetEta_bxm1_firstbunch")
    JetEtabxm1_f.SetStats(False)
    JetEtabxm1_f.SetLineColor(55)
    JetEtabxm1_f.SetMarkerColor(55)
    JetEtabxm1_f.SetMarkerStyle(22)
    JetEtabxm1_f.Draw("E same")

    JetEta_bx0_bxm1_f = inFile.Get("demo/JetEta_bx0_bxm1_firstbunch")
    JetEta_bx0_bxm1_f.SetStats(False)
    JetEta_bx0_bxm1_f.SetLineColor(8)
    JetEta_bx0_bxm1_f.SetMarkerColor(8)
    JetEta_bx0_bxm1_f.SetMarkerStyle(33)
    JetEta_bx0_bxm1_f.Draw("E same")

    JetEta_bxm2_f = inFile.Get("demo/JetEta_bxm2_firstbunch")
    JetEta_bxm2_f.SetStats(False)
    JetEta_bxm2_f.SetLineColor(6)
    JetEta_bxm2_f.SetMarkerColor(6)
    JetEta_bxm2_f.SetMarkerStyle(24)
    JetEta_bxm2_f.Draw("E same")

    legend2 = ROOT.TLegend(0.3,0.7,0.7,0.85)
    legend2.AddEntry(JetEtabx0_f,"BX=0","PE")
    legend2.AddEntry(JetEtabxm1_f,"BX=-1","PE")
    legend2.AddEntry(JetEta_bx0_bxm1_f,"BX=0 or BX=-1","PE")
    legend2.AddEntry(JetEta_bxm2_f,"BX=-2","PE")
    legend2.SetBorderSize(0)
    legend2.SetTextSize(0.04)
    legend2.Draw()

    c.SetGrid()
    c.SetLogy()
    c.Draw()
    c.SaveAs(f"{out_dir}/jeteta_firstbunchintrain.png")
    c.Close()
    
    # RATIO
    c_eta_rat = ROOT.TCanvas( 'c_eta_rat', 'c_eta_rat', 100, 10, 600, 600 )

    JetEta_bxm1_u = inFile.Get("demo/JetEta_bxm1_unprefirable")
    JetEta_bx0_bxm1_u = inFile.Get("demo/JetEta_bx0_bxm1_unprefirable")

    JetEta_ratio_u = JetEta_bxm1_u.Clone()
    JetEta_ratio_u.Divide(JetEta_bx0_bxm1_u)

    JetEta_ratio_u.SetLineColor(50)
    JetEta_ratio_u.SetMarkerColor(50)
    JetEta_ratio_u.SetMarkerStyle(8)
    JetEta_ratio_u.SetMarkerSize(0.5)
    JetEta_ratio_u.SetTitle("")
    JetEta_ratio_u.GetXaxis().SetTitle("Reco Jet #eta (GeV)")
    JetEta_ratio_u.GetYaxis().SetTitle("(bx=-1)/(bx=-1 or bx=0)")
    JetEta_ratio_u.SetMinimum(0)
    JetEta_ratio_u.SetMaximum(0.6)
    JetEta_ratio_u.Draw("E")

    JetEta_bxm1_f = inFile.Get("demo/JetEta_bxm1_firstbunch")
    JetEta_bx0_bxm1_f = inFile.Get("demo/JetEta_bx0_bxm1_firstbunch")

    JetEta_ratio_f = JetEta_bxm1_f.Clone()
    JetEta_ratio_f.Divide(JetEta_bx0_bxm1_f)

    JetEta_ratio_f.SetLineColor(9)
    JetEta_ratio_f.SetMarkerColor(9)
    JetEta_ratio_f.SetMarkerStyle(8)
    JetEta_ratio_f.SetMarkerSize(0.5)
    JetEta_ratio_f.Draw("E same")

    legend = ROOT.TLegend(0.5,0.8,0.9,0.88)
    legend.AddEntry(JetEta_ratio_u,"UnprefirableEvent","L")
    legend.AddEntry(JetEta_ratio_f,"FirstBunchInTrain","L")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c_eta_rat.Draw()
    c_eta_rat.SaveAs(f"{out_dir}/jeteta_ratio.png")
    c_eta_rat.Close()
    
    
    ############################################################################
    ############################# eta/phi plots ################################
    ############################################################################
    
    # UNPREFIRABLE
    c_ep_0_u = ROOT.TCanvas( 'c_ep_0_u', 'c_ep_0_u', 100, 10, 600, 600 )
    EtaPhi = inFile.Get("demo/JetEtaPhi_bx0_unprefirable")
    EtaPhi.GetXaxis().SetTitle("#eta")
    EtaPhi.GetYaxis().SetTitle("#phi")
    EtaPhi.SetTitle("Unprefirable (Offline) Jets BX=0 (Offline Jet pT>30 GeV)")
    EtaPhi.SetStats(False)
    EtaPhi.Draw("colz")
    c_ep_0_u.Draw()
    c_ep_0_u.SaveAs(f"{out_dir}/jetetaphi_bx0_offline_unprefirable.png")
    c_ep_0_u.Close()

    c_ep_m1_u = ROOT.TCanvas( 'c_ep_m1_u', 'c_ep_m1_u', 100, 10, 600, 600 )
    EtaPhi = inFile.Get("demo/JetEtaPhi_bxm1_unprefirable")
    EtaPhi.GetXaxis().SetTitle("#eta")
    EtaPhi.GetYaxis().SetTitle("#phi")
    EtaPhi.SetTitle("Unprefirable (Offline) Jets BX=-1 (Offline Jet pT>30 GeV)")
    EtaPhi.SetStats(False)
    EtaPhi.Draw("colz")
    c_ep_m1_u.Draw()
    c_ep_m1_u.SaveAs(f"{out_dir}/jetetaphi_bxm1_offline_unprefirable.png")
    c_ep_m1_u.Close()

    c_ep_0_u_on = ROOT.TCanvas( 'c_ep_0_u_on', 'c_ep_0_u_on', 100, 10, 600, 600 )
    EtaPhi = inFile.Get("demo/JetEtaPhi_bx0_online_unprefirable")
    EtaPhi.GetXaxis().SetTitle("#eta")
    EtaPhi.GetYaxis().SetTitle("#phi")
    EtaPhi.SetTitle("Unprefirable (Online) Jets BX=0 (Offline Jet pT>30 GeV)")
    EtaPhi.SetStats(False)
    EtaPhi.Draw("colz")
    c_ep_0_u_on.Draw()
    c_ep_0_u_on.SaveAs(f"{out_dir}/jetetaphi_bx0_online_unprefirable.png")
    c_ep_0_u_on.Close()

    c_ep_m1_u_on = ROOT.TCanvas( 'c_ep_m1_u_on', 'c_ep_m1_u_on', 100, 10, 600, 600 )
    EtaPhi = inFile.Get("demo/JetEtaPhi_bxm1_online_unprefirable")
    EtaPhi.GetXaxis().SetTitle("#eta")
    EtaPhi.GetYaxis().SetTitle("#phi")
    EtaPhi.SetTitle("Unprefirable (Online) Jets BX=-1 (Offline Jet pT>30 GeV)")
    EtaPhi.SetStats(False)
    EtaPhi.Draw("colz")
    c_ep_m1_u_on.Draw()
    c_ep_m1_u_on.SaveAs(f"{out_dir}/jetetaphi_bxm1_online_unprefirable.png")
    c_ep_m1_u_on.Close()
    
    # FIRSTBUNCHINTRAIN
    c_ep_0_f = ROOT.TCanvas( 'c_ep_0_f', 'c_ep_0_f', 100, 10, 600, 600 )
    EtaPhi = inFile.Get("demo/JetEtaPhi_bx0_firstbunch")
    EtaPhi.GetXaxis().SetTitle("#eta")
    EtaPhi.GetYaxis().SetTitle("#phi")
    EtaPhi.SetTitle("FirstBunchInTrain (Offline) Jets BX=0 (Offline Jet pT>30 GeV)")
    EtaPhi.SetStats(False)
    EtaPhi.Draw("colz")
    c_ep_0_f.Draw()
    c_ep_0_f.SaveAs(f"{out_dir}/jetetaphi_bx0_offline_firstbunch.png")
    c_ep_0_f.Close()

    c_ep_m1_f = ROOT.TCanvas( 'c_ep_m1_f', 'c_ep_m1_f', 100, 10, 600, 600 )
    EtaPhi = inFile.Get("demo/JetEtaPhi_bxm1_firstbunch")
    EtaPhi.GetXaxis().SetTitle("#eta")
    EtaPhi.GetYaxis().SetTitle("#phi")
    EtaPhi.SetTitle("FirstBunchInTrain (Offline) Jets BX=-1 (Offline Jet pT>30 GeV)")
    EtaPhi.SetStats(False)
    EtaPhi.Draw("colz")
    c_ep_m1_f.Draw()
    c_ep_m1_f.SaveAs(f"{out_dir}/jetetaphi_bxm1_offline_firstbunch.png")
    c_ep_m1_f.Close()

    c_ep_0_f_on = ROOT.TCanvas( 'c_ep_0_f_on', 'c_ep_0_f_on', 100, 10, 600, 600 )
    EtaPhi = inFile.Get("demo/JetEtaPhi_bx0_online_firstbunch")
    EtaPhi.GetXaxis().SetTitle("#eta")
    EtaPhi.GetYaxis().SetTitle("#phi")
    EtaPhi.SetTitle("FirstBunchInTrain (Online) Jets BX=0 (Offline Jet pT>30 GeV)")
    EtaPhi.SetStats(False)
    EtaPhi.Draw("colz")
    c_ep_0_f_on.Draw()
    c_ep_0_f_on.SaveAs(f"{out_dir}/jetetaphi_bx0_online_firstbunch.png")
    c_ep_0_f_on.Close()

    c_ep_m1_f_on = ROOT.TCanvas( 'c_ep_m1_f_on', 'c_ep_m1_f_on', 100, 10, 600, 600 )
    EtaPhi = inFile.Get("demo/JetEtaPhi_bxm1_online_firstbunch")
    EtaPhi.GetXaxis().SetTitle("#eta")
    EtaPhi.GetYaxis().SetTitle("#phi")
    EtaPhi.SetTitle("FirstBunchInTrain (Online) Jets BX=-1 (Offline Jet pT>30 GeV)")
    EtaPhi.SetStats(False)
    EtaPhi.Draw("colz")
    c_ep_m1_f_on.Draw()
    c_ep_m1_f_on.SaveAs(f"{out_dir}/jetetaphi_bxm1_online_firstbunch.png")
    c_ep_m1_f_on.Close()
    
    
    ################################################################################
    ################################# pt res plots #################################
    ################################################################################
    
    # unprefirable
    c1 = ROOT.TCanvas( 'c1', 'c1', 100, 10, 1000, 600 )
    JetEResbx0_u = inFile.Get("demo/JetERes_bx0_unprefirable")
    JetEResbx0_u.SetTitle("Unprefirable")
    JetEResbx0_u.SetLineColor(2)
    JetEResbx0_u.SetMarkerColor(2)
    JetEResbx0_u.SetMarkerStyle(47)
    JetEResbx0_u.GetXaxis().SetTitle("Online Jet p_{T} / Offline Jet p_{T}")
    JetEResbx0_u.GetYaxis().SetTitle("a.u.")
    JetEResbx0_u.SetStats(False)
    JetEResbx0_u.SetMaximum(1e2)
    JetEResbx0_u.Draw("E")

    JetEResbxm1_u = inFile.Get("demo/JetERes_bxm1_unprefirable")
    JetEResbxm1_u.SetStats(False)
    JetEResbxm1_u.SetLineColor(55)
    JetEResbxm1_u.SetMarkerColor(55)
    JetEResbxm1_u.SetMarkerStyle(22)
    JetEResbxm1_u.Draw("E same")

    JetEResbxm2_u = inFile.Get("demo/JetERes_bxm2_unprefirable")
    JetEResbxm2_u.SetStats(False)
    JetEResbxm2_u.SetLineColor(6)
    JetEResbxm2_u.SetMarkerColor(6)
    JetEResbxm2_u.SetMarkerStyle(23)
    JetEResbxm2_u.Draw("E same")

    JetEResbx1_u = inFile.Get("demo/JetERes_bx1_unprefirable")
    JetEResbx1_u.SetStats(False)
    JetEResbx1_u.SetLineColor(8)
    JetEResbx1_u.SetMarkerColor(8)
    JetEResbx1_u.SetMarkerStyle(42)
    JetEResbx1_u.Draw("E same")

    JetEResbx2_u = inFile.Get("demo/JetERes_bx2_unprefirable")
    JetEResbx2_u.SetStats(False)
    JetEResbx2_u.SetLineColor(95)
    JetEResbx2_u.SetMarkerColor(95)
    JetEResbx2_u.SetMarkerStyle(43)
    JetEResbx2_u.Draw("E same")

    legend = ROOT.TLegend(0.7,0.7,0.9,0.85)
    legend.AddEntry(JetEResbx0_u,"BX=0","PE")
    legend.AddEntry(JetEResbxm1_u,"BX=-1","PE")
    legend.AddEntry(JetEResbxm2_u,"BX=-2","PE")
    legend.AddEntry(JetEResbx1_u,"BX=1","PE")
    legend.AddEntry(JetEResbx2_u,"BX=2","PE")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c1.SetLogy()
    c1.Draw()
    c1.SaveAs(f"{out_dir}/jetptres_unprefirable.png")
    c1.Close()

    # firstbunchintrain
    c2 = ROOT.TCanvas( 'c2', 'c2', 100, 10, 1000, 600 )
    JetEResbx0_f = inFile.Get("demo/JetERes_bx0_firstbunch")
    JetEResbx0_f.SetTitle("FirstBunchInTrain")
    JetEResbx0_f.SetLineColor(2)
    JetEResbx0_f.SetMarkerColor(2)
    JetEResbx0_f.SetMarkerStyle(47)
    JetEResbx0_f.GetXaxis().SetTitle("Offline Jet p_{T} / Online Jet p_{T}")
    JetEResbx0_f.GetYaxis().SetTitle("a.u.")
    JetEResbx0_f.SetStats(False)
    JetEResbx0_f.SetMaximum(5e3)
    JetEResbx0_f.Draw("E")

    JetEResbxm1_f = inFile.Get("demo/JetERes_bxm1_firstbunch")
    JetEResbxm1_f.SetStats(False)
    JetEResbxm1_f.SetLineColor(55)
    JetEResbxm1_f.SetMarkerColor(55)
    JetEResbxm1_f.SetMarkerStyle(22)
    JetEResbxm1_f.Draw("E same")

    JetEResbxm2_f = inFile.Get("demo/JetERes_bxm2_firstbunch")
    JetEResbxm2_f.SetStats(False)
    JetEResbxm2_f.SetLineColor(6)
    JetEResbxm2_f.SetMarkerColor(6)
    JetEResbxm2_f.SetMarkerStyle(23)
    JetEResbxm2_f.Draw("E same")

    JetEResbx1_f = inFile.Get("demo/JetERes_bx1_firstbunch")
    JetEResbx1_f.SetStats(False)
    JetEResbx1_f.SetLineColor(8)
    JetEResbx1_f.SetMarkerColor(8)
    JetEResbx1_f.SetMarkerStyle(42)
    JetEResbx1_f.Draw("E same")

    JetEResbx2_f = inFile.Get("demo/JetERes_bx2_firstbunch")
    JetEResbx2_f.SetStats(False)
    JetEResbx2_f.SetLineColor(95)
    JetEResbx2_f.SetMarkerColor(95)
    JetEResbx2_f.SetMarkerStyle(43)
    JetEResbx2_f.Draw("E same")

    legend2 = ROOT.TLegend(0.7,0.7,0.9,0.85)
    legend2.AddEntry(JetEResbx0_f,"BX=0","PE")
    legend2.AddEntry(JetEResbxm1_f,"BX=-1","PE")
    legend2.AddEntry(JetEResbxm2_f,"BX=-2","PE")
    legend2.AddEntry(JetEResbx1_f,"BX=1","PE")
    legend2.AddEntry(JetEResbx2_f,"BX=2","PE")
    legend2.SetBorderSize(0)
    legend2.SetTextSize(0.04)
    legend2.Draw()

    c2.SetLogy()
    c2.Draw()
    c2.SaveAs(f"{out_dir}/jetptres_firstbunch.png")
    c2.Close()
    
    
if __name__ == "__main__":
    
    # parse arguments
    parser = argparse.ArgumentParser(description="This program creates plots to study jet pre-firing in the L1 Trigger.")
    parser.add_argument("-p", "--file_path", help="path to input ROOT file containing histograms")
    parser.add_argument("-o", "--out_dir", help="path to directory to save plots", default=".")
    
    args = parser.parse_args()
    
    main(args.file_path, args.out_dir)
