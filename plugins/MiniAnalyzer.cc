// -*- C++ -*-
//
// Package:    L1JetMini/MiniAnalyzer
// Class:      MiniAnalyzer
//
/**\class MiniAnalyzer MiniAnalyzer.cc L1JetMini/MiniAnalyzer/plugins/MiniAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Elliott Kauffman
//         Created:  Tue, 13 Jun 2023 09:23:09 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Jet.h"

#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"

#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"

#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

#include "TH1.h"
#include "TH2.h"

using namespace l1extra;
using namespace std;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


class MiniAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MiniAnalyzer(const edm::ParameterSet&);
  ~MiniAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  //edm::EDGetTokenT<BXVector<l1t::Jet>> stage2JetToken_;
  edm::EDGetTokenT<GlobalExtBlkBxCollection> UnprefirableEventToken_;
  const edm::EDGetTokenT<l1t::JetBxCollection> jetBXCollectionToken_;
  edm::EDGetTokenT< BXVector<GlobalAlgBlk> > gtAlgBlkToken;
  edm::Handle< BXVector<GlobalAlgBlk> > gtAlgBlkHandle;


  bool Flag_IsUnprefirable;
  bool Flag_FirstBunchInTrain;

  edm::ESGetToken<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd> L1TUtmTriggerMenuEventToken;
  const L1TUtmTriggerMenu* l1GtMenu;
  const std::map<std::string, L1TUtmAlgorithm>* algorithmMap;

  // define histograms
  TH1F *n1_u;
  TH1F *n2_u; 
  TH1F *n3_u;
  TH1F *n4_u;
  TH1F *n5_u;
  TH1F *n6_u;
  TH1F *n7_u;
  TH2F *n8_u;

  TH1F *h1_u;
  TH1F *h2_u;
  TH1F *h3_u;
  TH1F *h4_u;
  TH1F *h5_u;
  TH1F *h6_u;
  TH1F *h7_u;
  TH1F *h8_u;
  TH1F *h9_u;
  TH1F *h10_u;
  TH1F *h11_u;
  TH1F *h12_u;
  TH1F *h13_u;
  TH1F *h14_u;
  TH1F *h15_u;

  TH1F *n1_f;
  TH1F *n2_f;
  TH1F *n3_f;
  TH1F *n4_f;
  TH1F *n5_f;
  TH1F *n6_f;
  TH1F *n7_f;
  TH2F *n8_f;

  TH1F *h1_f;
  TH1F *h2_f;
  TH1F *h3_f;
  TH1F *h4_f;
  TH1F *h5_f;
  TH1F *h6_f;
  TH1F *h7_f;
  TH1F *h8_f;
  TH1F *h9_f;
  TH1F *h10_f;
  TH1F *h11_f;
  TH1F *h12_f;
  TH1F *h13_f;
  TH1F *h14_f;
  TH1F *h15_f;

  int nJets_f, nJets_u;
  float HT0_u, HTM1_u, HT0_f, HTM1_f;
  float jetEtBx0_u, jetEtBxM1_u, jetEtaBx0_u, jetEtaBxM1_u, jetPhiBx0_u, jetPhiBxM1_u, jetEtBxM2_u;
  float jetEtBx0_f, jetEtBxM1_f, jetEtaBx0_f, jetEtaBxM1_f, jetPhiBx0_f, jetPhiBxM1_f, jetEtBxM2_f;


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig):
  //stage2JetToken_(consumes<BXVector<l1t::Jet>>( edm::InputTag("caloStage2Digis","Jet","RECO"))),
  UnprefirableEventToken_(consumes<GlobalExtBlkBxCollection>(edm::InputTag("simGtExtUnprefireable"))),
  jetBXCollectionToken_(consumes<l1t::JetBxCollection>(edm::InputTag("caloStage2Digis","Jet","RECO"))),
  gtAlgBlkToken( consumes< BXVector<GlobalAlgBlk> >(edm::InputTag("gtStage2Digis","","RECO")) )

{
  L1TUtmTriggerMenuEventToken = consumesCollector().esConsumes<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd>();
  edm::Service<TFileService> fs;
  n1_u = fs->make<TH1F>("nJets_unprefirable","Number of jets",15,0,15);
  n2_u = fs->make<TH1F>("nJetTh_unprefirable","Jets passing E_{T} threshold", 5, -0.5, 4.5);
  n3_u = fs->make<TH1F>("nJetPre60_unprefirable","Jets passing Pre-trigger threshold 60 GeV", 5, -0.5, 4.5);
  n4_u = fs->make<TH1F>("nJetPre90_unprefirable","Jets passing Pre-trigger threshold 90 GeV", 5, -0.5, 4.5);
  n5_u = fs->make<TH1F>("nJetPre120_unprefirable","Jets passing Pre-trigger threshold 120 GeV", 5, -0.5, 4.5);
  n6_u = fs->make<TH1F>("nJetPre150_unprefirable","Jets passing Pre-trigger threshold 150 GeV", 5, -0.5, 4.5);
  n7_u = fs->make<TH1F>("nJetPre180_unprefirable","Jets passing Pre-trigger threshold 180 GeV", 5, -0.5, 4.5);
  n8_u = fs->make<TH2F>("EtaPhi_unprefirable","#eta vs #phi of jets passing Pre-trigger threshold", 40, -5, 5, 40, -M_PI, M_PI);
  
  h1_u = fs->make<TH1F>("JetEtbx0_unprefirable","Jet E_{T} at BX = 0", 40, 0, 250);
  h2_u = fs->make<TH1F>("JetEtbxm1_unprefirable","Jet E_{T} at BX = -1", 40, 0, 250);
  h3_u = fs->make<TH1F>("LeadJetEtbx0_unprefirable","Leading Jet E_{T} at BX = 0", 80, 0, 500);
  h4_u = fs->make<TH1F>("LeadJetEtbxm1_unprefirable","Leading Jet E_{T} at BX = -1", 80, 0, 500);
  h5_u = fs->make<TH1F>("LeadJetEtbxm2_unprefirable","Leading Jet E_{T} at BX = -2", 80, 0, 500);
  h6_u = fs->make<TH1F>("LeadJetEtabx0_unprefirable","Leading Jet #eta at BX = 0", 40, -5, 5);
  h7_u = fs->make<TH1F>("LeadJetEtabxm1_unprefirable","Leading Jet #eta at BX = -1", 40, -5, 5);
  h8_u = fs->make<TH1F>("nJetHT_unprefirable","Jets sum passing H_{T} threshold", 5, -0.5, 4.5);
  h9_u = fs->make<TH1F>("nJetHTPre150_unprefirable","Jets sum passing Pre-trigger H_{T} threshold 150 GeV", 5, -0.5, 4.5);
  h10_u = fs->make<TH1F>("nJetHTPre200_unprefirable","Jets sum passing Pre-trigger H_{T} threshold 200 GeV", 5, -0.5, 4.5);
  h11_u = fs->make<TH1F>("nJetHTPre250_unprefirable","Jets sum passing Pre-trigger H_{T} threshold 250 GeV", 5, -0.5, 4.5);
  h12_u = fs->make<TH1F>("nJetHTPre300_unprefirable","Jets sum passing Pre-trigger H_{T} threshold 300 GeV", 5, -0.5, 4.5);
  h13_u = fs->make<TH1F>("nJetHTPre350_unprefirable","Jets sum passing Pre-trigger H_{T} threshold 350 GeV", 5, -0.5, 4.5);
  h14_u = fs->make<TH1F>("HTbx0_unprefirable","H_{T} at at BX = 0", 80, 0, 1000.);
  h15_u = fs->make<TH1F>("HTbxm1_unprefirable","H_{T} at at BX = -1", 80, 0, 1000.);

  n1_f = fs->make<TH1F>("nJets_firstbunch","Number of jets",15,0,15);
  n2_f = fs->make<TH1F>("nJetTh_firstbunch","Jets passing E_{T} threshold", 5, -0.5, 4.5);
  n3_f = fs->make<TH1F>("nJetPre60_firstbunch","Jets passing Pre-trigger threshold 60 GeV", 5, -0.5, 4.5);
  n4_f = fs->make<TH1F>("nJetPre90_firstbunch","Jets passing Pre-trigger threshold 90 GeV", 5, -0.5, 4.5);
  n5_f = fs->make<TH1F>("nJetPre120_firstbunch","Jets passing Pre-trigger threshold 120 GeV", 5, -0.5, 4.5);
  n6_f = fs->make<TH1F>("nJetPre150_firstbunch","Jets passing Pre-trigger threshold 150 GeV", 5, -0.5, 4.5);
  n7_f = fs->make<TH1F>("nJetPre180_firstbunch","Jets passing Pre-trigger threshold 180 GeV", 5, -0.5, 4.5);
  n8_f = fs->make<TH2F>("EtaPhi_firstbunch","#eta vs #phi of jets passing Pre-trigger threshold", 40, -5, 5, 40, -M_PI, M_PI);

  h1_f = fs->make<TH1F>("JetEtbx0_firstbunch","Jet E_{T} at BX = 0", 40, 0, 250);
  h2_f = fs->make<TH1F>("JetEtbxm1_firstbunch","Jet E_{T} at BX = -1", 40, 0, 250);
  h3_f = fs->make<TH1F>("LeadJetEtbx0_firstbunch","Leading Jet E_{T} at BX = 0", 80, 0, 500);
  h4_f = fs->make<TH1F>("LeadJetEtbxm1_firstbunch","Leading Jet E_{T} at BX = -1", 80, 0, 500);
  h5_f = fs->make<TH1F>("LeadJetEtbxm2_firstbunch","Leading Jet E_{T} at BX = -2", 80, 0, 500);
  h6_f = fs->make<TH1F>("LeadJetEtabx0_firstbunch","Leading Jet #eta at BX = 0", 40, -5, 5);
  h7_f = fs->make<TH1F>("LeadJetEtabxm1_firstbunch","Leading Jet #eta at BX = -1", 40, -5, 5);
  h8_f = fs->make<TH1F>("nJetHT_firstbunch","Jets sum passing H_{T} threshold", 5, -0.5, 4.5);
  h9_f = fs->make<TH1F>("nJetHTPre150_firstbunch","Jets sum passing Pre-trigger H_{T} threshold 150 GeV", 5, -0.5, 4.5);
  h10_f = fs->make<TH1F>("nJetHTPre200_firstbunch","Jets sum passing Pre-trigger H_{T} threshold 200 GeV", 5, -0.5, 4.5);
  h11_f = fs->make<TH1F>("nJetHTPre250_firstbunch","Jets sum passing Pre-trigger H_{T} threshold 250 GeV", 5, -0.5, 4.5);
  h12_f = fs->make<TH1F>("nJetHTPre300_firstbunch","Jets sum passing Pre-trigger H_{T} threshold 300 GeV", 5, -0.5, 4.5);
  h13_f = fs->make<TH1F>("nJetHTPre350_firstbunch","Jets sum passing Pre-trigger H_{T} threshold 350 GeV", 5, -0.5, 4.5);
  h14_f = fs->make<TH1F>("HTbx0_firstbunch","H_{T} at at BX = 0", 80, 0, 1000.);
  h15_f = fs->make<TH1F>("HTbxm1_firstbunch","H_{T} at at BX = -1", 80, 0, 1000.);

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

MiniAnalyzer::~MiniAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty

  h8_u->Sumw2(); h9_u->Sumw2(); h10_u->Sumw2(); h11_u->Sumw2(); h12_u->Sumw2(); h13_u->Sumw2();
  h9_u->Divide(h8_u); h10_u->Divide(h8_u); h11_u->Divide(h8_u); h12_u->Divide(h8_u); h13_u->Divide(h8_u);
   
  h8_f->Sumw2(); h9_f->Sumw2(); h10_f->Sumw2(); h11_f->Sumw2(); h12_f->Sumw2(); h13_f->Sumw2();
  h9_f->Divide(h8_f); h10_f->Divide(h8_f); h11_f->Divide(h8_f); h12_f->Divide(h8_f); h13_f->Divide(h8_f);

  n2_u->Sumw2(); n3_u->Sumw2(); n4_u->Sumw2(); n5_u->Sumw2(); n6_u->Sumw2(); n7_u->Sumw2();
  n3_u->Divide(n2_u); n4_u->Divide(n2_u); n5_u->Divide(n2_u); n6_u->Divide(n2_u); n7_u->Divide(n2_u);

  n2_f->Sumw2(); n3_f->Sumw2(); n4_f->Sumw2(); n5_f->Sumw2(); n6_f->Sumw2(); n7_f->Sumw2();
  n3_f->Divide(n2_f); n4_f->Divide(n2_f); n5_f->Divide(n2_f); n6_f->Divide(n2_f); n7_f->Divide(n2_f);
}

//
// member functions
//

// ------------ method called for each event  ------------
void MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  
  auto menuRcd = iSetup.get<L1TUtmTriggerMenuRcd>();
  l1GtMenu = &menuRcd.get(L1TUtmTriggerMenuEventToken);
  algorithmMap = &(l1GtMenu->getAlgorithmMap());
  
  //FirstBunchInTrain
  Flag_FirstBunchInTrain = false;
  iEvent.getByToken(gtAlgBlkToken, gtAlgBlkHandle);
  if(gtAlgBlkHandle.isValid()){
    std::vector<GlobalAlgBlk>::const_iterator algBlk = gtAlgBlkHandle->begin(0);
    // std::cout<<"BxVector First BX: "<<gtAlgBlkHandle->getFirstBX()<<std::endl;
    // std::cout<<"BxVector Last BX: "<<gtAlgBlkHandle->getLastBX()<<std::endl;
    if(algBlk != gtAlgBlkHandle->end(0)){
      for (std::map<std::string, L1TUtmAlgorithm>::const_iterator itAlgo = algorithmMap->begin(); itAlgo != algorithmMap->end(); itAlgo++) {
  std::string algName = itAlgo->first;
        int algBit = itAlgo->second.getIndex();
        bool initialDecision = algBlk->getAlgoDecisionInitial(algBit);
        if(algBit==473 && initialDecision==1){
          Flag_FirstBunchInTrain = true;
        }
      }
    }
  }

  if(Flag_FirstBunchInTrain){
    
    nJets_f = 0;
    HT0_f = 0.; HTM1_f = 0.;
    jetEtBx0_f=0.; jetEtBxM1_f=0.; jetEtaBx0_f=-99.; jetEtaBxM1_f=-99.; jetPhiBx0_f=-99.; jetPhiBxM1_f=-99.; jetEtBxM2_f=0.;
  
    edm::Handle<l1t::JetBxCollection> jetColl;
    iEvent.getByToken(jetBXCollectionToken_, jetColl);
    l1t::JetBxCollection jets;
    jets = (*jetColl.product());

    // jets in bx=0
    for (auto it = jets.begin(0); it!=jets.end(0); it++){
      
      nJets_f = nJets_f + 1;
      
      h1_f->Fill(it->pt());
      HT0_f=HT0_f+it->pt();
      
      if(it->pt() > jetEtBx0_f){
        jetEtBx0_f = it->pt();
        jetEtaBx0_f = it->eta();
        jetPhiBx0_f = it->phi();
      }

    }
    h14_f->Fill(HT0_f);

    // jets in bx=-1
    for (auto it = jets.begin(-1); it!=jets.end(-1); it++){
      
      nJets_f = nJets_f + 1;

      h2_f->Fill(it->pt());
      HTM1_f=HTM1_f+it->pt();
    
      if(it->pt() > jetEtBxM1_f){
        jetEtBxM1_f = it->pt();
        jetEtaBxM1_f = it->eta();
        jetPhiBxM1_f = it->phi();
      }

    }
    h15_f->Fill(HTM1_f); 

    // jets in bx=-2
    for (auto it = jets.begin(-2); it!=jets.end(-2); it++){
      
      nJets_f = nJets_f + 1;

      if(it->pt() > jetEtBxM2_f){
        jetEtBxM2_f = it->pt(); 
      }

    }

    // jets in bx=1
    for (auto it = jets.begin(1); it!=jets.end(1); it++){
      nJets_f = nJets_f + 1;
    }
    // jets in bx=2
    for (auto it = jets.begin(2); it!=jets.end(2); it++){
      nJets_f = nJets_f + 1;
    }

    //HT
    for(int ibin = 0; ibin < 5; ibin++){   //loop over # of thresholds
      float j = 300. + ibin*50.;    //HT thresholds at BX=0
      std::string s = std::to_string(j);
      char const *pchar = s.c_str();
      h8_f->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h9_f->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h10_f->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h11_f->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h12_f->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h13_f->GetXaxis()->SetBinLabel(ibin+1, pchar);
      if(HT0_f > j){
        h8_f->Fill(ibin);   //fill the BX=0 denominator
        float k = 150. + ibin*50.;    //HT thresholds at BX=-1
        if(HTM1_f > k){   //fill the BX=-1 numerators
          h9_f->Fill(ibin);
          h10_f->Fill(ibin);
          h11_f->Fill(ibin);
          h12_f->Fill(ibin);
          h13_f->Fill(ibin);
        }
      }
    }

    if(jetEtBx0_f > 0.){
      if(jetEtBx0_f > 120. && jetEtBxM1_f > 60.) n8_f->Fill(jetEtaBxM1_f, jetPhiBxM1_f);
      for(int ibin = 0; ibin < 5; ibin++){   //loop over # of thresholds
        float j = 120. + ibin*30.;    //jet pt thresholds at BX=0
        std::string s = std::to_string(j);
        char const *pchar = s.c_str();
        n2_f->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n3_f->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n4_f->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n5_f->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n6_f->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n7_f->GetXaxis()->SetBinLabel(ibin+1, pchar);

        if(jetEtBx0_f > j){
          n2_f->Fill(ibin);   //fill the BX=0 denominator
          float k = 60. + ibin*30.;    //jet pt thresholds at BX=-1
          if(jetEtBxM1_f > k){   //fill the BX=-1 numerators
            n8_f->Fill(jetEtaBxM1_f, jetPhiBxM1_f);
            if(jetEtBxM1_f > k) n3_f->Fill(ibin);
            if(jetEtBxM1_f > k) n4_f->Fill(ibin);
            if(jetEtBxM1_f > k) n5_f->Fill(ibin);
            if(jetEtBxM1_f > k) n6_f->Fill(ibin);
            if(jetEtBxM1_f > k) n7_f->Fill(ibin);
          }
        }
      }
    }

    if(jetEtBxM1_f > 60. ) std::cout<<"jetEtBx0 (FirstBunchInTrain): "<<jetEtBx0_f<<"\t"<<"jetEtBxM1 (FirstBunchInTrain): "<<jetEtBxM1_f<<std::endl;
    n1_f->Fill(nJets_f);
    h3_f->Fill(jetEtBx0_f);
    h4_f->Fill(jetEtBxM1_f);
    h5_f->Fill(jetEtBxM2_f);
    h6_f->Fill(jetEtaBx0_f);
    h7_f->Fill(jetEtaBxM1_f);

  }

  //Unprefirable
  Flag_IsUnprefirable = false;
  edm::Handle<GlobalExtBlkBxCollection> handleUnprefEventResults;
  iEvent.getByToken(UnprefirableEventToken_, handleUnprefEventResults);
  if(handleUnprefEventResults.isValid()){
    if (handleUnprefEventResults->size() != 0) {
      Flag_IsUnprefirable = handleUnprefEventResults->at(0, 0).getExternalDecision(GlobalExtBlk::maxExternalConditions-1);
    }
  }
  
  if(Flag_IsUnprefirable){

    nJets_u = 0;
    HT0_u = 0.; HTM1_u = 0.;
    jetEtBx0_u=0.; jetEtBxM1_u=0.; jetEtaBx0_u=-99.; jetEtaBxM1_u=-99.; jetPhiBx0_u=-99.; jetPhiBxM1_u=-99.; jetEtBxM2_u=0.;
    
    edm::Handle<l1t::JetBxCollection> jetColl;
    iEvent.getByToken(jetBXCollectionToken_, jetColl);
    l1t::JetBxCollection jets;
    jets = (*jetColl.product()); 

    // jets in bx=0
    for (auto it = jets.begin(0); it!=jets.end(0); it++){
      
      nJets_u = nJets_u + 1;
      
      h1_u->Fill(it->pt());
      HT0_u=HT0_u+it->pt();
      
      if(it->pt() > jetEtBx0_u){
        jetEtBx0_u = it->pt();
        jetEtaBx0_u = it->eta();
        jetPhiBx0_u = it->phi();
      }

    }
    h14_u->Fill(HT0_u);

    // jets in bx=-1
    for (auto it = jets.begin(-1); it!=jets.end(-1); it++){
      
      nJets_u = nJets_u + 1;

      h2_u->Fill(it->pt());
      HTM1_u=HTM1_u+it->pt();
    
      if(it->pt() > jetEtBxM1_u){
        jetEtBxM1_u = it->pt();
        jetEtaBxM1_u = it->eta();
        jetPhiBxM1_u = it->phi();
      }

    }
    h15_u->Fill(HTM1_u); 

    // jets in bx=-2
    for (auto it = jets.begin(-2); it!=jets.end(-2); it++){
      
      nJets_u = nJets_u + 1;

      if(it->pt() > jetEtBxM2_u){
        jetEtBxM2_u = it->pt(); 
      }

    }

    // jets in bx=1
    for (auto it = jets.begin(1); it!=jets.end(1); it++){
      nJets_u = nJets_u + 1;
    }
    // jets in bx=2
    for (auto it = jets.begin(2); it!=jets.end(2); it++){
      nJets_u = nJets_u + 1;
    }
   
    //HT
    for(int ibin = 0; ibin < 5; ibin++){   //loop over # of thresholds
      float j = 300. + ibin*50.;    //HT thresholds at BX=0
      std::string s = std::to_string(j);
      char const *pchar = s.c_str();
      h8_u->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h9_u->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h10_u->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h11_u->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h12_u->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h13_u->GetXaxis()->SetBinLabel(ibin+1, pchar);
      if(HT0_u > j){
        h8_u->Fill(ibin);   //fill the BX=0 denominator
        float k = 150. + ibin*50.;    //HT thresholds at BX=-1
        if(HTM1_u > k){   //fill the BX=-1 numerators
          h9_u->Fill(ibin);
          h10_u->Fill(ibin);
          h11_u->Fill(ibin);
          h12_u->Fill(ibin);
          h13_u->Fill(ibin);
        }
      }
    }
  
    if(jetEtBx0_u > 0.){
      if(jetEtBx0_u > 120. && jetEtBxM1_u > 60.) n8_u->Fill(jetEtaBxM1_u, jetPhiBxM1_u);
      for(int ibin = 0; ibin < 5; ibin++){   //loop over # of thresholds
        float j = 120. + ibin*30.;    //jet pt thresholds at BX=0
        std::string s = std::to_string(j);
        char const *pchar = s.c_str();
        n2_u->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n3_u->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n4_u->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n5_u->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n6_u->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n7_u->GetXaxis()->SetBinLabel(ibin+1, pchar);

        if(jetEtBx0_u > j){
          n2_u->Fill(ibin);   //fill the BX=0 denominator
          float k = 60. + ibin*30.;    //jet pt thresholds at BX=-1
          if(jetEtBxM1_u > k){   //fill the BX=-1 numerators
            n8_u->Fill(jetEtaBxM1_u, jetPhiBxM1_u);
            if(jetEtBxM1_u > k) n3_u->Fill(ibin);
            if(jetEtBxM1_u > k) n4_u->Fill(ibin);
            if(jetEtBxM1_u > k) n5_u->Fill(ibin);
            if(jetEtBxM1_u > k) n6_u->Fill(ibin);
            if(jetEtBxM1_u > k) n7_u->Fill(ibin);
          }
        }
      }
    }
   
    if(jetEtBxM1_u > 60. ) std::cout<<"jetEtBx0 (Unprefirable): "<<jetEtBx0_u<<"\t"<<"jetEtBxM1 (Unprefirable): "<<jetEtBxM1_u<<std::endl;
    n1_u->Fill(nJets_u);
    h3_u->Fill(jetEtBx0_u);
    h4_u->Fill(jetEtBxM1_u);
    h5_u->Fill(jetEtBxM2_u);
    h6_u->Fill(jetEtaBx0_u);
    h7_u->Fill(jetEtaBxM1_u);
  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void MiniAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void MiniAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MiniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);