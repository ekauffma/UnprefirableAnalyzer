import argparse
import awkward as ak
import numpy as np
import ROOT
import uproot


def main(file_path, out_dir):
    
    inFile = ROOT.TFile.Open(file_path, "READ")

    f = uproot.open(file_path)
    print(f.keys())
 
    ###################################################################################
    ########################## dijet mass pass fail plot ##############################
    ###################################################################################

    c_mjj = ROOT.TCanvas( 'c_mjj', 'c_mjj', 10, 10, 800, 600 )

    mjj_pass = inFile.Get("demo/mjj_passL1FinalORinBXm1_unprefirable")
    mjj_fail = inFile.Get("demo/mjj_failL1FinalORinBXm1_unprefirable")

    mjj_demoniminator = mjj_fail.Clone()
    mjj_demoniminator.Add(mjj_pass)
    mjj_demoniminator.Sumw2()
    mjj_ratio = mjj_pass.Clone()
    mjj_ratio.Sumw2()
    mjj_ratio.Divide(mjj_demoniminator)

    mjj_ratio.SetLineColor(50)
    mjj_ratio.SetMarkerColor(50)
    mjj_ratio.SetMarkerStyle(8)
    mjj_ratio.SetMarkerSize(1)
    mjj_ratio.SetTitle("")
    mjj_ratio.GetXaxis().SetTitle("Dijet Mass (GeV)")
    mjj_ratio.GetYaxis().SetTitle("Fraction of events passing L1FinalOR in BX=-1")
    mjj_ratio.SetMinimum(0)
    mjj_ratio.SetMaximum(0.5)
    mjj_ratio.SetStats(False)
    mjj_ratio.Draw("E")

    c_mjj.SetGrid()
    c_mjj.Draw()
    c_mjj.SaveAs(f"{out_dir}/mjj_ratio.png")
    c_mjj.Close()

    ###################################################################################
    ######################### dijet deltaR pass fail plot #############################
    ###################################################################################

    c_deltaR = ROOT.TCanvas( 'c_deltaR', 'c_deltaR', 0, 10, 700, 500 )

    deltaR_pass = inFile.Get("demo/deltaR_passL1FinalORinBXm1_unprefirable")
    deltaR_fail = inFile.Get("demo/deltaR_failL1FinalORinBXm1_unprefirable")

    deltaR_demoniminator = deltaR_fail.Clone()
    deltaR_demoniminator.Add(deltaR_pass)
    deltaR_demoniminator.Sumw2()
    deltaR_ratio = deltaR_pass.Clone()
    deltaR_ratio.Sumw2()
    deltaR_ratio.Divide(deltaR_demoniminator)

    deltaR_ratio.SetLineColor(50)
    deltaR_ratio.SetMarkerColor(50)
    deltaR_ratio.SetMarkerStyle(8)
    deltaR_ratio.SetMarkerSize(1)
    deltaR_ratio.SetTitle("")
    deltaR_ratio.GetXaxis().SetTitle("Dijet $\Delta R$")
    deltaR_ratio.GetYaxis().SetTitle("Fraction of events passing L1FinalOR in BX=-1")
    deltaR_ratio.GetXaxis().SetRangeUser(2.6,3.5)
    deltaR_ratio.SetMinimum(0)
    deltaR_ratio.SetMaximum(0.03)
    deltaR_ratio.SetStats(False)
    deltaR_ratio.Draw("E")

    c_deltaR.SetGrid()
    c_deltaR.Draw()
    c_deltaR.SaveAs(f"{out_dir}/deltaR_ratio.png")
    c_deltaR.Close()

    ###################################################################################
    ######################## dijet deltaPhi pass fail plot ############################
    ###################################################################################

    c_deltaPhi = ROOT.TCanvas( 'c_deltaPhi', 'c_deltaPhi', 0, 10, 700, 500)

    deltaPhi_pass = inFile.Get("demo/deltaPhi_passL1FinalORinBXm1_unprefirable")
    deltaPhi_fail = inFile.Get("demo/deltaPhi_failL1FinalORinBXm1_unprefirable")

    deltaPhi_demoniminator = deltaPhi_fail.Clone()
    deltaPhi_demoniminator.Add(deltaPhi_pass)
    deltaPhi_demoniminator.Sumw2()
    deltaPhi_ratio = deltaPhi_pass.Clone()
    deltaPhi_ratio.Sumw2()
    deltaPhi_ratio.Divide(deltaPhi_demoniminator)

    deltaPhi_ratio.SetLineColor(50)
    deltaPhi_ratio.SetMarkerColor(50)
    deltaPhi_ratio.SetMarkerStyle(8)
    deltaPhi_ratio.SetMarkerSize(1)
    deltaPhi_ratio.SetTitle("")
    deltaPhi_ratio.GetXaxis().SetTitle("Dijet $\Delta \phi$")
    deltaPhi_ratio.GetYaxis().SetTitle("Fraction of events passing L1FinalOR in BX=-1")
    deltaPhi_ratio.GetXaxis().SetRangeUser(2.68,3.14)
    deltaPhi_ratio.SetMinimum(0)
    deltaPhi_ratio.SetMaximum(0.03)
    deltaPhi_ratio.SetStats(False)
    deltaPhi_ratio.Draw("E")

    c_deltaPhi.SetGrid()
    c_deltaPhi.Draw()
    c_deltaPhi.SaveAs(f"{out_dir}/deltaPhi_ratio.png")
    c_deltaPhi.Close()

    ###################################################################################
    ######################## dijet deltaEta pass fail plot ############################
    ###################################################################################

    #c_deltaEta = ROOT.TCanvas( 'c_deltaEta', 'c_deltaEta', 10, 10, 800, 600 )

    #deltaEta_pass = inFile.Get("demo/deltaEta_passL1FinalORinBXm1_unprefirable")
    #deltaEta_fail = inFile.Get("demo/deltaEta_failL1FinalORinBXm1_unprefirable")

    #deltaEta_demoniminator = deltaEta_fail.Clone()
    #deltaEta_demoniminator.Add(deltaEta_pass)
    #deltaEta_demoniminator.Sumw2()
    #deltaEta_ratio = deltaEta_pass.Clone()
    #deltaEta_ratio.Sumw2()
    #deltaEta_ratio.Divide(deltaEta_demoniminator)

    #deltaEta_ratio.SetLineColor(50)
    #deltaEta_ratio.SetMarkerColor(50)
    #deltaEta_ratio.SetMarkerStyle(8)
    #deltaEta_ratio.SetMarkerSize(0.5)
    #deltaEta_ratio.SetTitle("")
    #deltaEta_ratio.GetXaxis().SetTitle("Dijet $\Delta \eta$")
    #deltaEta_ratio.GetYaxis().SetTitle("Fraction of events passing L1FinalOR in BX=-1")
    #deltaEta_ratio.SetMinimum(0)
    #deltaEta_ratio.SetMaximum(0.05)
    #deltaEta_ratio.SetStats(False)
    #deltaEta_ratio.Draw("E")

    #c_deltaEta.Draw()
    #c_deltaEta.SaveAs(f"{out_dir}/deltaEta_ratio.png")
    #c_deltaEta.Close()
  
    ###################################################################################
    ##################### dijet leading jet pt pass fail plot #########################
    ###################################################################################   

    c_pt0 = ROOT.TCanvas( 'c_pt0', 'c_pt0', 10, 10, 800, 600 )

    pt0_pass = inFile.Get("demo/pt0_passL1FinalORinBXm1_unprefirable")
    pt0_fail = inFile.Get("demo/pt0_failL1FinalORinBXm1_unprefirable")

    pt0_demoniminator = pt0_fail.Clone()
    pt0_demoniminator.Add(pt0_pass)
    pt0_demoniminator.Sumw2()
    pt0_ratio = pt0_pass.Clone()
    pt0_ratio.Sumw2()
    pt0_ratio.Divide(pt0_demoniminator)

    pt0_ratio.SetLineColor(50)
    pt0_ratio.SetMarkerColor(50)
    pt0_ratio.SetMarkerStyle(8)
    pt0_ratio.SetMarkerSize(1)
    pt0_ratio.SetTitle("")
    pt0_ratio.GetXaxis().SetTitle("Leading Jet pT (GeV)")
    pt0_ratio.GetYaxis().SetTitle("Fraction of events passing L1FinalOR in BX=-1")
    pt0_ratio.SetMinimum(0)
    pt0_ratio.SetMaximum(0.7)
    pt0_ratio.SetStats(False)
    pt0_ratio.Draw("E")

    c_pt0.SetGrid()
    c_pt0.Draw()
    c_pt0.SaveAs(f"{out_dir}/pt0_ratio.png")
    c_pt0.Close()

    ###################################################################################
    ################### dijet subleading jet pt pass fail plot ########################
    ###################################################################################   

    c_pt1 = ROOT.TCanvas( 'c_pt1', 'c_pt1', 10, 10, 800, 600 )

    pt1_pass = inFile.Get("demo/pt1_passL1FinalORinBXm1_unprefirable")
    pt1_fail = inFile.Get("demo/pt1_failL1FinalORinBXm1_unprefirable")

    pt1_demoniminator = pt1_fail.Clone()
    pt1_demoniminator.Add(pt1_pass)
    pt1_demoniminator.Sumw2()
    pt1_ratio = pt1_pass.Clone()
    pt1_ratio.Sumw2()
    pt1_ratio.Divide(pt1_demoniminator)

    pt1_ratio.SetLineColor(50)
    pt1_ratio.SetMarkerColor(50)
    pt1_ratio.SetMarkerStyle(8)
    pt1_ratio.SetMarkerSize(1)
    pt1_ratio.SetTitle("")
    pt1_ratio.GetXaxis().SetTitle("Subleading Jet pT (GeV)")
    pt1_ratio.GetYaxis().SetTitle("Fraction of events passing L1FinalOR in BX=-1")
    pt1_ratio.SetMinimum(0)
    pt1_ratio.SetMaximum(0.7)
    pt1_ratio.SetStats(False)
    pt1_ratio.Draw("E")

    c_pt1.SetGrid()
    c_pt1.Draw()
    c_pt1.SaveAs(f"{out_dir}/pt1_ratio.png")
    c_pt1.Close()

    ###################################################################################
    #################### dijet leading jet eta pass fail plot #########################
    ###################################################################################   

    c_eta0 = ROOT.TCanvas( 'c_eta0', 'c_eta0',0, 10, 700, 500)

    eta0_pass = inFile.Get("demo/eta0_passL1FinalORinBXm1_unprefirable")
    eta0_fail = inFile.Get("demo/eta0_failL1FinalORinBXm1_unprefirable")

    eta0_demoniminator = eta0_fail.Clone()
    eta0_demoniminator.Add(eta0_pass)
    eta0_demoniminator.Sumw2()
    eta0_ratio = eta0_pass.Clone()
    eta0_ratio.Sumw2()
    eta0_ratio.Divide(eta0_demoniminator)

    eta0_ratio.SetLineColor(50)
    eta0_ratio.SetMarkerColor(50)
    eta0_ratio.SetMarkerStyle(8)
    eta0_ratio.SetMarkerSize(1)
    eta0_ratio.SetTitle("")
    eta0_ratio.GetXaxis().SetTitle("Leading Jet $\eta$")
    eta0_ratio.GetYaxis().SetTitle("Fraction of events passing L1FinalOR in BX=-1")
    eta0_ratio.SetMinimum(0)
    eta0_ratio.SetMaximum(0.03)
    eta0_ratio.SetStats(False)
    eta0_ratio.Draw("E")

    c_eta0.SetGrid()
    c_eta0.Draw()
    c_eta0.SaveAs(f"{out_dir}/eta0_ratio.png")
    c_eta0.Close()

    ###################################################################################
    ################### dijet subleading jet eta pass fail plot #######################
    ###################################################################################   

    c_eta1 = ROOT.TCanvas( 'c_eta1', 'c_eta1', 0, 10, 700, 500)

    eta1_pass = inFile.Get("demo/eta1_passL1FinalORinBXm1_unprefirable")
    eta1_fail = inFile.Get("demo/eta1_failL1FinalORinBXm1_unprefirable")

    eta1_demoniminator = eta1_fail.Clone()
    eta1_demoniminator.Add(eta1_pass)
    eta1_demoniminator.Sumw2()
    eta1_ratio = eta1_pass.Clone()
    eta1_ratio.Sumw2()
    eta1_ratio.Divide(eta1_demoniminator)

    eta1_ratio.SetLineColor(50)
    eta1_ratio.SetMarkerColor(50)
    eta1_ratio.SetMarkerStyle(8)
    eta1_ratio.SetMarkerSize(1)
    eta1_ratio.SetTitle("")
    eta1_ratio.GetXaxis().SetTitle("Subleading Jet $\eta$")
    eta1_ratio.GetYaxis().SetTitle("Fraction of events passing L1FinalOR in BX=-1")
    eta1_ratio.SetMinimum(0)
    eta1_ratio.SetMaximum(0.03)
    eta1_ratio.SetStats(False)
    eta1_ratio.Draw("E")

    c_eta1.SetGrid()
    c_eta1.Draw()
    c_eta1.SaveAs(f"{out_dir}/eta1_ratio.png")
    c_eta1.Close()

    ###################################################################################
    #################### dijet leading jet phi pass fail plot #########################
    ###################################################################################   

    c_phi0 = ROOT.TCanvas( 'c_phi0', 'c_phi0', 0, 10, 700, 500 )

    phi0_pass = inFile.Get("demo/phi0_passL1FinalORinBXm1_unprefirable")
    phi0_fail = inFile.Get("demo/phi0_failL1FinalORinBXm1_unprefirable")

    phi0_demoniminator = phi0_fail.Clone()
    phi0_demoniminator.Add(phi0_pass)
    phi0_demoniminator.Sumw2()
    phi0_ratio = phi0_pass.Clone()
    phi0_ratio.Sumw2()
    phi0_ratio.Divide(phi0_demoniminator)

    phi0_ratio.SetLineColor(50)
    phi0_ratio.SetMarkerColor(50)
    phi0_ratio.SetMarkerStyle(8)
    phi0_ratio.SetMarkerSize(1)
    phi0_ratio.SetTitle("")
    phi0_ratio.GetXaxis().SetTitle("Leading Jet $\eta$")
    phi0_ratio.GetYaxis().SetTitle("Fraction of events passing L1FinalOR in BX=-1")
    phi0_ratio.SetMinimum(0)
    phi0_ratio.SetMaximum(0.03)
    phi0_ratio.SetStats(False)
    phi0_ratio.Draw("E")

    c_phi0.SetGrid()
    c_phi0.Draw()
    c_phi0.SaveAs(f"{out_dir}/phi0_ratio.png")
    c_phi0.Close()

    ###################################################################################
    #################### dijet leading jet phi pass fail plot #########################
    ###################################################################################   

    c_phi1 = ROOT.TCanvas( 'c_phi1', 'c_phi1', 0, 10, 700, 500)

    phi1_pass = inFile.Get("demo/phi1_passL1FinalORinBXm1_unprefirable")
    phi1_fail = inFile.Get("demo/phi1_failL1FinalORinBXm1_unprefirable")

    phi1_demoniminator = phi1_fail.Clone()
    phi1_demoniminator.Add(phi1_pass)
    phi1_demoniminator.Sumw2()
    phi1_ratio = phi1_pass.Clone()
    phi1_ratio.Sumw2()
    phi1_ratio.Divide(phi1_demoniminator)

    phi1_ratio.SetLineColor(50)
    phi1_ratio.SetMarkerColor(50)
    phi1_ratio.SetMarkerStyle(8)
    phi1_ratio.SetMarkerSize(1)
    phi1_ratio.SetTitle("")
    phi1_ratio.GetXaxis().SetTitle("Leading Jet $\phi$")
    phi1_ratio.GetYaxis().SetTitle("Fraction of events passing L1FinalOR in BX=-1")
    phi1_ratio.SetMinimum(0)
    phi1_ratio.SetMaximum(0.03)
    phi1_ratio.SetStats(False)
    phi1_ratio.Draw("E")

    c_phi1.SetGrid()
    c_phi1.Draw()
    c_phi1.SaveAs(f"{out_dir}/phi1_ratio.png")
    c_phi1.Close()
   
    # UNPREFIRABLE RATIO SEPARATED PT
    c_eta_rat_u = ROOT.TCanvas( 'c_eta_rat_u', 'c_eta_rat_u', 10, 10, 800, 600 )

    JetEta_low_bxm1_u = inFile.Get("demo/JetEta_lowpt_bxm1_unprefirable")
    JetEta_low_bx0_bxm1_u = inFile.Get("demo/JetEta_lowpt_bx0_bxm1_unprefirable")

    JetEta_low_ratio_u = JetEta_low_bxm1_u.Clone()

    JetEta_low_ratio_u.Sumw2()
    JetEta_low_bx0_bxm1_u.Sumw2()
    JetEta_low_ratio_u.Divide(JetEta_low_bx0_bxm1_u)

    JetEta_low_ratio_u.SetLineColor(50)
    JetEta_low_ratio_u.SetMarkerColor(50)
    JetEta_low_ratio_u.SetMarkerStyle(8)
    JetEta_low_ratio_u.SetMarkerSize(0.5)
    JetEta_low_ratio_u.SetTitle("Unprefirable")
    JetEta_low_ratio_u.GetXaxis().SetTitle("Reco Jet #eta")
    JetEta_low_ratio_u.GetYaxis().SetTitle("(bx=-1)/(bx=-1 or bx=0)")
    JetEta_low_ratio_u.SetMinimum(0)
    JetEta_low_ratio_u.SetMaximum(0.25)
    JetEta_low_ratio_u.SetStats(False)
    JetEta_low_ratio_u.Draw("E")

    JetEta_med_bxm1_u = inFile.Get("demo/JetEta_medpt_bxm1_unprefirable")
    JetEta_med_bx0_bxm1_u = inFile.Get("demo/JetEta_medpt_bx0_bxm1_unprefirable")

    JetEta_med_ratio_u = JetEta_med_bxm1_u.Clone()

    JetEta_med_ratio_u.Sumw2()
    JetEta_med_bx0_bxm1_u.Sumw2()
    JetEta_med_ratio_u.Divide(JetEta_med_bx0_bxm1_u)

    JetEta_med_ratio_u.SetLineColor(54)
    JetEta_med_ratio_u.SetMarkerColor(54)
    JetEta_med_ratio_u.SetMarkerStyle(29)
    JetEta_med_ratio_u.SetMarkerSize(0.5)
    JetEta_med_ratio_u.SetStats(False);
    JetEta_med_ratio_u.Draw("E same")

    JetEta_high_bxm1_u = inFile.Get("demo/JetEta_high_bxm1_unprefirable")
    JetEta_high_bx0_bxm1_u = inFile.Get("demo/JetEta_high_bx0_bxm1_unprefirable")

    JetEta_high_ratio_u = JetEta_high_bxm1_u.Clone()

    JetEta_high_ratio_u.Sumw2()
    JetEta_high_bx0_bxm1_u.Sumw2()
    JetEta_high_ratio_u.Divide(JetEta_high_bx0_bxm1_u)

    JetEta_high_ratio_u.SetLineColor(108)
    JetEta_high_ratio_u.SetMarkerColor(108)
    JetEta_high_ratio_u.SetMarkerStyle(21)
    JetEta_high_ratio_u.SetMarkerSize(0.5)
    JetEta_high_ratio_u.SetStats(False);
    JetEta_high_ratio_u.Draw("E same")

    legend = ROOT.TLegend(0.3,0.75,0.6,0.88)
    legend.AddEntry(JetEta_low_ratio_u,"L1 p_T > 15 GeV","L")
    legend.AddEntry(JetEta_med_ratio_u,"L1 p_T > 40 GeV","L")
    legend.AddEntry(JetEta_high_ratio_u,"L1 p_T > 180 GeV","L")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c_eta_rat_u.Draw()
    c_eta_rat_u.SaveAs(f"{out_dir}/jeteta_ratio_ptsep_unprefirable.png")
    c_eta_rat_u.Close() 
   
    # FIRSTBUNCHINTRAIN RATIO SEPARATED PT
    c_eta_rat_f = ROOT.TCanvas( 'c_eta_rat_f', 'c_eta_rat_f', 10, 10, 800, 600 )

    JetEta_low_bxm1_f = inFile.Get("demo/JetEta_lowpt_bxm1_firstbunch")
    JetEta_low_bx0_bxm1_f = inFile.Get("demo/JetEta_lowpt_bx0_bxm1_firstbunch")

    JetEta_low_ratio_f = JetEta_low_bxm1_f.Clone()

    JetEta_low_ratio_f.Sumw2()
    JetEta_low_bx0_bxm1_f.Sumw2()
    JetEta_low_ratio_f.Divide(JetEta_low_bx0_bxm1_f)

    JetEta_low_ratio_f.SetLineColor(50)
    JetEta_low_ratio_f.SetMarkerColor(50)
    JetEta_low_ratio_f.SetMarkerStyle(8)
    JetEta_low_ratio_f.SetMarkerSize(0.5)
    JetEta_low_ratio_f.SetTitle("FirstBunchInTrain")
    JetEta_low_ratio_f.GetXaxis().SetTitle("Reco Jet #eta")
    JetEta_low_ratio_f.GetYaxis().SetTitle("(bx=-1)/(bx=-1 or bx=0)")
    JetEta_low_ratio_f.SetMinimum(0)
    JetEta_low_ratio_f.SetMaximum(0.037)
    JetEta_low_ratio_f.SetStats(False)
    JetEta_low_ratio_f.Draw("E")

    JetEta_med_bxm1_f = inFile.Get("demo/JetEta_medpt_bxm1_firstbunch")
    JetEta_med_bx0_bxm1_f = inFile.Get("demo/JetEta_medpt_bx0_bxm1_firstbunch")

    JetEta_med_ratio_f = JetEta_med_bxm1_f.Clone()

    JetEta_med_ratio_f.Sumw2()
    JetEta_med_bx0_bxm1_f.Sumw2()
    JetEta_med_ratio_f.Divide(JetEta_med_bx0_bxm1_f)

    JetEta_med_ratio_f.SetLineColor(54)
    JetEta_med_ratio_f.SetMarkerColor(54)
    JetEta_med_ratio_f.SetMarkerStyle(29)
    JetEta_med_ratio_f.SetMarkerSize(0.5)
    JetEta_med_ratio_f.SetStats(False)
    JetEta_med_ratio_f.Draw("E same")

    JetEta_high_bxm1_f = inFile.Get("demo/JetEta_high_bxm1_firstbunch")
    JetEta_high_bx0_bxm1_f = inFile.Get("demo/JetEta_high_bx0_bxm1_firstbunch")

    JetEta_high_ratio_f = JetEta_high_bxm1_f.Clone()

    JetEta_high_ratio_f.Sumw2()
    JetEta_high_bx0_bxm1_f.Sumw2()
    JetEta_high_ratio_f.Divide(JetEta_high_bx0_bxm1_f)

    JetEta_high_ratio_f.SetLineColor(108)
    JetEta_high_ratio_f.SetMarkerColor(108)
    JetEta_high_ratio_f.SetMarkerStyle(21)
    JetEta_high_ratio_f.SetMarkerSize(0.5)
    JetEta_high_ratio_f.SetStats(False)
    JetEta_high_ratio_f.Draw("E same")

    legend = ROOT.TLegend(0.3,0.75,0.6,0.88)
    legend.AddEntry(JetEta_low_ratio_f,"L1 p_T > 15 GeV","L")
    legend.AddEntry(JetEta_med_ratio_f,"L1 p_T > 40 GeV","L")
    legend.AddEntry(JetEta_high_ratio_f,"L1 p_T > 180 GeV","L")
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.Draw()

    c_eta_rat_f.Draw()
    c_eta_rat_f.SaveAs(f"{out_dir}/jeteta_ratio_ptsep_firstbunch.png")
    c_eta_rat_f.Close()
 
 
    # RATIO
    c_etaphi_rat_u = ROOT.TCanvas( 'c_etaphi_rat_u', 'c_etaphi_rat_u', 10, 10, 800, 600 )

    JetEtaPhi_bxm1_u = inFile.Get("demo/JetEtaPhi_bxm1_unprefirable")
    JetEtaPhi_bx0_bxm1_u = inFile.Get("demo/JetEtaPhi_bx0_bxm1_unprefirable")
    JetEtaPhi_ratio_u = JetEtaPhi_bxm1_u.Clone()

    JetEtaPhi_ratio_u.Sumw2()
    JetEtaPhi_bx0_bxm1_u.Sumw2()
    JetEtaPhi_ratio_u.Divide(JetEtaPhi_bx0_bxm1_u)
    JetEtaPhi_ratio_u.SetStats(False)
    JetEtaPhi_ratio_u.SetTitle("Unprefirable (bx=-1)/(bx=0 or bx=-1) for Jets with pT>180 GeV")
    JetEtaPhi_ratio_u.GetXaxis().SetTitle("Reco Jet #eta")
    JetEtaPhi_ratio_u.GetYaxis().SetTitle("Reco Jet #phi")
    JetEtaPhi_ratio_u.Rebin2D()
    JetEtaPhi_ratio_u.Draw("colz")

    c_etaphi_rat_u.Draw()
    c_etaphi_rat_u.SaveAs(f"{out_dir}/jetetaphi_ratio_unprefirable.png")
    c_etaphi_rat_u.Close() 

    c_etaphi_rat_f = ROOT.TCanvas( 'c_etaphi_rat_f', 'c_etaphi_rat_f', 10, 10, 800, 600 )

    JetEtaPhi_bxm1_f = inFile.Get("demo/JetEtaPhi_bxm1_firstbunch")
    JetEtaPhi_bx0_bxm1_f = inFile.Get("demo/JetEtaPhi_bx0_bxm1_firstbunch")
    JetEtaPhi_ratio_f = JetEtaPhi_bxm1_f.Clone()

    JetEtaPhi_ratio_f.Sumw2()
    JetEtaPhi_bx0_bxm1_f.Sumw2()
    JetEtaPhi_ratio_f.Divide(JetEtaPhi_bx0_bxm1_f)
    JetEtaPhi_ratio_f.SetStats(False)
    JetEtaPhi_ratio_f.SetTitle("FirstBunchInTrain (bx=-1)/(bx=0 or bx=-1) for Jets with pT>180 GeV")
    JetEtaPhi_ratio_f.GetXaxis().SetTitle("Reco Jet #eta")
    JetEtaPhi_ratio_f.GetYaxis().SetTitle("Reco Jet #phi")
    JetEtaPhi_ratio_f.Draw("colz")

    c_etaphi_rat_f.Draw()
    c_etaphi_rat_f.SaveAs(f"{out_dir}/jetetaphi_ratio_firstbunch.png")
    c_etaphi_rat_f.Close()
    
    
if __name__ == "__main__":
    
    # parse arguments
    parser = argparse.ArgumentParser(description="This program creates plots to study jet pre-firing in the L1 Trigger.")
    parser.add_argument("-p", "--file_path", help="path to input ROOT file containing histograms")
    parser.add_argument("-o", "--out_dir", help="path to directory to save plots", default=".")
    
    args = parser.parse_args()
    
    main(args.file_path, args.out_dir)
