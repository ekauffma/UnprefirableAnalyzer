#ifndef postan_h
#define postan_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <stdio.h>
#include <TString.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TH2F.h>

using namespace std;

#include "vector"

class postan {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TFile *fileName;
   TH1F *h_jetet_bx0_u, *h_jetet_bxm1_u, *h_jetet_bx0_bxm1_u, *h_jetet_bxm2_u, *h_jetet_bx0_bxm2_u;
   TH1F *h_jeteta_bx0_u, *h_jeteta_bxm1_u, *h_jeteta_bx0_bxm1_u, *h_jeteta_bxm2_u, *h_jeteta_bx0_bxm2_u;
   TH2F *h_jetetaphi_bx0_u, *h_jetetaphi_bxm1_u, *h_jetetaphi_bx0_bxm1_u, *h_jetetaphi_bxm2_u, *h_jetetaphi_bx0_bxm2_u;
   TH1F *h_jetet_bx0_f, *h_jetet_bxm1_f, *h_jetet_bx0_bxm1_f, *h_jetet_bxm2_f, *h_jetet_bx0_bxm2_f;
   TH1F *h_jeteta_bx0_f, *h_jeteta_bxm1_f, *h_jeteta_bx0_bxm1_f, *h_jeteta_bxm2_f, *h_jeteta_bx0_bxm2_f;
   TH2F *h_jetetaphi_bx0_f, *h_jetetaphi_bxm1_f, *h_jetetaphi_bx0_bxm1_f, *h_jetetaphi_bxm2_f, *h_jetetaphi_bx0_bxm2_f;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run_num;
   Int_t           lumi;
   Int_t           event_num;
   Int_t           isUnprefirable;
   Int_t           isFirstBunchInTrain;
   Int_t           L1FinalOR_bxm1;
   Int_t           idx_L1_FirstBunchBeforeTrain;
   vector<bool>    *trigger_bits;
   vector<TLorentzVector> *reco_jets;
   vector<bool>    *reco_jetId;
   vector<TLorentzVector> *match_l1_bx0;
   vector<TLorentzVector> *match_l1_bxm1;
   vector<TLorentzVector> *match_l1_bxm2;
   vector<TLorentzVector> *match_l1_bx1;
   vector<TLorentzVector> *match_l1_bx2;

   // List of branches
   TBranch        *b_run_num;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event_num;   //!
   TBranch        *b_isUnprefirable;   //!
   TBranch        *b_isFirstBunchInTrain;   //!
   TBranch        *b_L1FinalOR_bxm1;   //!
   TBranch        *b_idx_L1_FirstBunchBeforeTrain;   //!
   TBranch        *b_trigger_bits;   //!
   TBranch        *b_reco_jets;   //!
   TBranch        *b_reco_jetId;   //!
   TBranch        *b_match_l1_bx0;   //!
   TBranch        *b_match_l1_bxm1;   //!
   TBranch        *b_match_l1_bxm2;   //!
   TBranch        *b_match_l1_bx1;   //!
   TBranch        *b_match_l1_bx2;   //!

   postan(string fileToOpen, const char* file2);
   virtual ~postan();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     BookHistos(const char* file2);
};

#endif

#ifdef postan_cxx
postan::postan(string fileToOpen, const char* file2) 
{
  
   BookHistos(file2);
   TChain *chain = new TChain("demo/eventTree");
   ifstream file;
   file.open(fileToOpen.c_str(), ifstream::in );
   char filename[2000];
   while (true) {
     file >> filename;
     if( file.eof() ) break;
         chain->Add(filename);
         cout<<"Added "<<filename<<endl;
   }//loop over while
   Init(chain);
}

postan::~postan()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   fileName->cd();
   fileName->Write();
   fileName->Close();
}

Int_t postan::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t postan::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void postan::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trigger_bits = 0;
   reco_jets = 0;
   reco_jetId = 0;
   match_l1_bx0 = 0;
   match_l1_bxm1 = 0;
   match_l1_bxm2 = 0;
   match_l1_bx1 = 0;
   match_l1_bx2 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run_num", &run_num, &b_run_num);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event_num", &event_num, &b_event_num);
   fChain->SetBranchAddress("isUnprefirable", &isUnprefirable, &b_isUnprefirable);
   fChain->SetBranchAddress("isFirstBunchInTrain", &isFirstBunchInTrain, &b_isFirstBunchInTrain);
   fChain->SetBranchAddress("L1FinalOR_bxm1", &L1FinalOR_bxm1, &b_L1FinalOR_bxm1);
   fChain->SetBranchAddress("idx_L1_FirstBunchBeforeTrain", &idx_L1_FirstBunchBeforeTrain, &b_idx_L1_FirstBunchBeforeTrain);
   fChain->SetBranchAddress("trigger_bits", &trigger_bits, &b_trigger_bits);
   fChain->SetBranchAddress("reco_jets", &reco_jets, &b_reco_jets);
   fChain->SetBranchAddress("reco_jetId", &reco_jetId, &b_reco_jetId);
   fChain->SetBranchAddress("match_l1_bx0", &match_l1_bx0, &b_match_l1_bx0);
   fChain->SetBranchAddress("match_l1_bxm1", &match_l1_bxm1, &b_match_l1_bxm1);
   fChain->SetBranchAddress("match_l1_bxm2", &match_l1_bxm2, &b_match_l1_bxm2);
   fChain->SetBranchAddress("match_l1_bx1", &match_l1_bx1, &b_match_l1_bx1);
   fChain->SetBranchAddress("match_l1_bx2", &match_l1_bx2, &b_match_l1_bx2);
   Notify();
}

Bool_t postan::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void postan::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t postan::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef postan_cxx
