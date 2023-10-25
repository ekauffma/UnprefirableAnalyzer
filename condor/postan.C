#define postan_cxx
#include "postan.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

float phiAbs(float phi1, float phi2){
  float dp = abs(phi1-phi2);
  if (dp > M_PI){
    dp -= 2*M_PI;
  }
  return abs(dp);
}

int main(int argc, char *argv[])
{

  if(argc > 1)
    {
      postan t(argv[1], argv[2]);
      t.Loop();
    }
  return 0;
}

using namespace std;

void postan::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      // if (Cut(ientry) < 0) continue;

      //int njet = 0;
      //for(unsigned int i = 0; i < reco_jets->size(); i++){
      //   if((*reco_jetId)[i] && njet == 0) {
      //      TLorentzVector* j1 = (TLorentzVector*) (*reco_jets)[i];
      //      njet++;
      //   }
      //   else if((*reco_jetId)[i] && njet == 1) {
      //      TLorentzVector* j2 = (TLorentzVector*) (*reco_jets)[i];
      //      njet++;
      //   }
      //   else if (njet == 2) break;
      //}
   
      //if(njet == 2){
      //   if(((j1+j2)->M() > 500) && (phiAbs(j1->Phi(), j2->Phi()) > 2.7) && (abs(j1->Eta() - j2->Eta()) < 1.3)){
            if(isUnprefirable == 1){
               for(unsigned int i = 0; i < reco_jets->size(); i++){
                  if(match_l1_bx0->size() > 0){
                     if((*match_l1_bx0)[i].Pt() > 180.) { 
                        h_jetet_bx0_u->Fill((*reco_jets)[i].Pt()); 
                        h_jeteta_bx0_u->Fill((*reco_jets)[i].Eta()); 
                        h_jetetaphi_bx0_u->Fill((*reco_jets)[i].Eta(), (*reco_jets)[i].Phi()); 
                     }
                  }
                  if(match_l1_bxm1->size() > 0){
                     if((*match_l1_bxm1)[i].Pt() > 180.) { 
                        h_jetet_bxm1_u->Fill((*reco_jets)[i].Pt()); 
                        h_jeteta_bxm1_u->Fill((*reco_jets)[i].Eta()); 
                        h_jetetaphi_bxm1_u->Fill((*reco_jets)[i].Eta(), (*reco_jets)[i].Phi()); 
                     }
                  }
                  if(match_l1_bx0->size() > 0 && match_l1_bxm1->size() > 0){
                     if((*match_l1_bx0)[i].Pt() > 180. || (*match_l1_bxm1)[i].Pt() > 180.) { 
                        h_jetet_bx0_bxm1_u->Fill((*reco_jets)[i].Pt()); 
                        h_jeteta_bx0_bxm1_u->Fill((*reco_jets)[i].Eta());  
                        h_jetetaphi_bx0_bxm1_u->Fill((*reco_jets)[i].Eta(), (*reco_jets)[i].Phi()); 
                     }
                  }
                  if(match_l1_bxm2->size() > 0){
                     if((*match_l1_bxm2)[i].Pt() > 180.) { 
                        h_jetet_bxm2_u->Fill((*reco_jets)[i].Pt()); 
                        h_jeteta_bxm2_u->Fill((*reco_jets)[i].Eta()); 
                        h_jetetaphi_bxm2_u->Fill((*reco_jets)[i].Eta(), (*reco_jets)[i].Phi()); 
                     }
                  }
                  if(match_l1_bx0->size() > 0 && match_l1_bxm2->size() > 0){
                     if((*match_l1_bx0)[i].Pt() > 180. || (*match_l1_bxm2)[i].Pt() > 180.) { 
                        h_jetet_bx0_bxm2_u->Fill((*reco_jets)[i].Pt()); 
                        h_jeteta_bx0_bxm2_u->Fill((*reco_jets)[i].Eta()); 
                        h_jetetaphi_bx0_bxm2_u->Fill((*reco_jets)[i].Eta(), (*reco_jets)[i].Phi()); 
                     }              
                  }
               }
            } // unprefirable
            if(isFirstBunchInTrain == 1){
               for(unsigned int i = 0; i < reco_jets->size(); i++){
                  if(match_l1_bx0->size() > 0){
                     if((*match_l1_bx0)[i].Pt() > 180.) {
                        h_jetet_bx0_f->Fill((*reco_jets)[i].Pt());
                        h_jeteta_bx0_f->Fill((*reco_jets)[i].Eta());
                        h_jetetaphi_bx0_f->Fill((*reco_jets)[i].Eta(), (*reco_jets)[i].Phi());
                     }
                  }
                  if(match_l1_bxm1->size() > 0){
                     if((*match_l1_bxm1)[i].Pt() > 180.) {
                        h_jetet_bxm1_f->Fill((*reco_jets)[i].Pt());
                        h_jeteta_bxm1_f->Fill((*reco_jets)[i].Eta());
                        h_jetetaphi_bxm1_f->Fill((*reco_jets)[i].Eta(), (*reco_jets)[i].Phi());
                     }
                  }
                  if(match_l1_bx0->size() > 0 && match_l1_bxm1->size() > 0){
                     if((*match_l1_bx0)[i].Pt() > 180. || (*match_l1_bxm1)[i].Pt() > 180.) {
                        h_jetet_bx0_bxm1_f->Fill((*reco_jets)[i].Pt());
                        h_jeteta_bx0_bxm1_f->Fill((*reco_jets)[i].Eta());
                        h_jetetaphi_bx0_bxm1_f->Fill((*reco_jets)[i].Eta(), (*reco_jets)[i].Phi()); 
                     }
                  }
                  if(match_l1_bxm2->size() > 0){
                     if((*match_l1_bxm2)[i].Pt() > 180.) { 
                        h_jetet_bxm2_f->Fill((*reco_jets)[i].Pt()); 
                        h_jeteta_bxm2_f->Fill((*reco_jets)[i].Eta()); 
                        h_jetetaphi_bxm2_f->Fill((*reco_jets)[i].Eta(), (*reco_jets)[i].Phi()); 
                     }
                  }  
                  if(match_l1_bx0->size() > 0 && match_l1_bxm2->size() > 0){
                     if((*match_l1_bx0)[i].Pt() > 180. || (*match_l1_bxm2)[i].Pt() > 180.) { 
                        h_jetet_bx0_bxm2_f->Fill((*reco_jets)[i].Pt()); 
                        h_jeteta_bx0_bxm2_f->Fill((*reco_jets)[i].Eta()); 
                        h_jetetaphi_bx0_bxm2_f->Fill((*reco_jets)[i].Eta(), (*reco_jets)[i].Phi()); 
                     }
                  }
               }
            } // firstbunchintrain
      //   } // dijet condition
      //} //njet condition
   }
}

void postan::BookHistos(const char* file2){
   fileName = new TFile(file2, "RECREATE");
   fileName->cd();

   h_jetet_bx0_u = new TH1F("JetEt_bx0_unprefirable","Jet E_{T} at BX = 0", 40, 0, 250);
   h_jetet_bxm1_u = new TH1F("JetEt_bxm1_unprefirable","Jet E_{T} at BX = -1", 40, 0, 250);
   h_jetet_bx0_bxm1_u = new TH1F("JetEt_bx0_bxm1_unprefirable","Jet E_{T} at BX = 0 or -1", 40, 0, 250);
   h_jetet_bxm2_u = new TH1F("JetEt_bxm2_unprefirable","Jet E_{T} at BX = -2", 40, 0, 250);
   h_jetet_bx0_bxm2_u = new TH1F("JetEt_bx0_bxm2_unprefirable","Jet E_{T} at BX = 0 or -2", 40, 0, 250);

   h_jeteta_bx0_u = new TH1F("JetEta_bx0_unprefirable","Jet E_{T} at BX = 0", 40, -5, 5);
   h_jeteta_bxm1_u = new TH1F("JetEta_bxm1_unprefirable","Jet E_{T} at BX = -1", 40, -5, 5);
   h_jeteta_bx0_bxm1_u = new TH1F("JetEta_bx0_bxm1_unprefirable","Jet E_{T} at BX = 0 or -1", 40, -5, 5);
   h_jeteta_bxm2_u = new TH1F("JetEta_bxm2_unprefirable","Jet E_{T} at BX = -2", 40, -5, 5);
   h_jeteta_bx0_bxm2_u = new TH1F("JetEta_bx0_bxm2_unprefirable","Jet E_{T} at BX = 0 or -2", 40, -5, 5);

   h_jetetaphi_bx0_u = new TH2F("JetEtaPhi_bx0_unprefirable","#eta vs #phi of jets with p_T>180 GeV (BX=0)",40, -5, 5, 40, -M_PI, M_PI);
   h_jetetaphi_bxm1_u = new TH2F("JetEtaPhi_bxm1_unprefirable","#eta vs #phi of jets with p_T>180 GeV (BX=-1)",40, -5, 5, 40, -M_PI, M_PI);
   h_jetetaphi_bx0_bxm1_u = new TH2F("JetEtaPhi_bx0_bxm1_unprefirable","#eta vs #phi of jets with p_T>180 GeV (BX=0 or -1)",40, -5, 5, 40, -M_PI, M_PI);
   h_jetetaphi_bxm2_u = new TH2F("JetEtaPhi_bxm2_unprefirable","#eta vs #phi of jets with p_T>180 GeV (BX=-2)",40, -5, 5, 40, -M_PI, M_PI);
   h_jetetaphi_bx0_bxm2_u = new TH2F("JetEtaPhi_bx0_bxm2_unprefirable","#eta vs #phi of jets with p_T>180 GeV (BX=0 or -2)",40, -5, 5, 40, -M_PI, M_PI); 

   h_jetet_bx0_f = new TH1F("JetEt_bx0_firstbunch","Jet E_{T} at BX = 0", 40, 0, 250);
   h_jetet_bxm1_f = new TH1F("JetEt_bxm1_firstbunch","Jet E_{T} at BX = -1", 40, 0, 250);
   h_jetet_bx0_bxm1_f = new TH1F("JetEt_bx0_bxm1_firstbunch","Jet E_{T} at BX = 0 or -1", 40, 0, 250);
   h_jetet_bxm2_f = new TH1F("JetEt_bxm2_firstbunch","Jet E_{T} at BX = -2", 40, 0, 250);
   h_jetet_bx0_bxm2_f = new TH1F("JetEt_bx0_bxm2_firstbunch","Jet E_{T} at BX = 0 or -2", 40, 0, 250);

   h_jeteta_bx0_f = new TH1F("JetEta_bx0_firstbunch","Jet E_{T} at BX = 0", 40, -5, 5);
   h_jeteta_bxm1_f = new TH1F("JetEta_bxm1_firstbunch","Jet E_{T} at BX = -1", 40, -5, 5);
   h_jeteta_bx0_bxm1_f = new TH1F("JetEta_bx0_bxm1_firstbunch","Jet E_{T} at BX = 0 or -1", 40, -5, 5);
   h_jeteta_bxm2_f = new TH1F("JetEta_bxm2_firstbunch","Jet E_{T} at BX = -2", 40, -5, 5);
   h_jeteta_bx0_bxm2_f = new TH1F("JetEta_bx0_bxm2_firstbunch","Jet E_{T} at BX = 0 or -2", 40, -5, 5);

   h_jetetaphi_bx0_f= new TH2F("JetEtaPhi_bx0_firstbunch","#eta vs #phi of jets with p_T>180 GeV (BX=0)",40, -5, 5, 40, -M_PI, M_PI);
   h_jetetaphi_bxm1_f = new TH2F("JetEtaPhi_bxm1_firstbunch","#eta vs #phi of jets with p_T>180 GeV (BX=-1)",40, -5, 5, 40, -M_PI, M_PI);
   h_jetetaphi_bx0_bxm1_f = new TH2F("JetEtaPhi_bx0_bxm1_firstbunch","#eta vs #phi of jets with p_T>180 GeV (BX=0 or -1)",40, -5, 5, 40, -M_PI, M_PI);
   h_jetetaphi_bxm2_f = new TH2F("JetEtaPhi_bxm2_firstbunch","#eta vs #phi of jets with p_T>180 GeV (BX=-2)",40, -5, 5, 40, -M_PI, M_PI);
   h_jetetaphi_bx0_bxm2_f = new TH2F("JetEtaPhi_bx0_bxm2_firstbunch","#eta vs #phi of jets with p_T>180 GeV (BX=0 or -2)",40, -5, 5, 40, -M_PI, M_PI);

}
