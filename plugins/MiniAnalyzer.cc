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

  edm::ESGetToken<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd> L1TUtmTriggerMenuEventToken;
  const L1TUtmTriggerMenu* l1GtMenu;
  const std::map<std::string, L1TUtmAlgorithm>* algorithmMap;

  // define histograms
  TH1F *n1;
  TH1F *n2; 
  TH1F *n3;
  TH1F *n4;
  TH1F *n5;
  TH1F *n6;
  TH1F *n7;
  TH2F *n8;

  TH1F *h1;
  TH1F *h2;
  TH1F *h3;
  TH1F *h4;
  TH1F *h5;
  TH1F *h6;
  TH1F *h7;
  TH1F *h8;
  TH1F *h9;
  TH1F *h10;
  TH1F *h11;
  TH1F *h12;
  TH1F *h13;
  TH1F *h14;
  TH1F *h15;

  int nJets;
  float HT0, HTM1;
  float jetEtBx0, jetEtBxM1, jetEtaBx0, jetEtaBxM1, jetPhiBx0, jetPhiBxM1, jetEtBxM2;


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
  n1 = fs->make<TH1F>("nJets","Number of jets",15,0,15);
  n2 = fs->make<TH1F>("nJetTh","Jets passing E_{T} threshold", 5, -0.5, 4.5);
  n3 = fs->make<TH1F>("nJetPre60","Jets passing Pre-trigger threshold 60 GeV", 5, -0.5, 4.5);
  n4 = fs->make<TH1F>("nJetPre90","Jets passing Pre-trigger threshold 90 GeV", 5, -0.5, 4.5);
  n5 = fs->make<TH1F>("nJetPre120","Jets passing Pre-trigger threshold 120 GeV", 5, -0.5, 4.5);
  n6 = fs->make<TH1F>("nJetPre150","Jets passing Pre-trigger threshold 150 GeV", 5, -0.5, 4.5);
  n7 = fs->make<TH1F>("nJetPre180","Jets passing Pre-trigger threshold 180 GeV", 5, -0.5, 4.5);
  n8 = fs->make<TH2F>("EtaPhi","#eta vs #phi of jets passing Pre-trigger threshold", 40, -5, 5, 40, -M_PI, M_PI);
  
  h1 = fs->make<TH1F>("JetEtbx0","Jet E_{T} at BX = 0", 40, 0, 250);
  h2 = fs->make<TH1F>("JetEtbxm1","Jet E_{T} at BX = -1", 40, 0, 250);
  h3 = fs->make<TH1F>("LeadJetEtbx0","Leading Jet E_{T} at BX = 0", 40, 0, 500);
  h4 = fs->make<TH1F>("LeadJetEtbxm1","Leading Jet E_{T} at BX = -1", 40, 0, 200);
  h5 = fs->make<TH1F>("LeadJetEtbxm2","Leading Jet E_{T} at BX = -2", 40, 0, 200);
  h6 = fs->make<TH1F>("LeadJetEtabx0","Leading Jet #eta at BX = 0", 40, -5, 5);
  h7 = fs->make<TH1F>("LeadJetEtabxm1","Leading Jet #eta at BX = -1", 40, -5, 5);
  h8 = fs->make<TH1F>("nJetHT","Jets sum passing H_{T} threshold", 5, -0.5, 4.5);
  h9 = fs->make<TH1F>("nJetHTPre150","Jets sum passing Pre-trigger H_{T} threshold 150 GeV", 5, -0.5, 4.5);
  h10 = fs->make<TH1F>("nJetHTPre200","Jets sum passing Pre-trigger H_{T} threshold 200 GeV", 5, -0.5, 4.5);
  h11 = fs->make<TH1F>("nJetHTPre250","Jets sum passing Pre-trigger H_{T} threshold 250 GeV", 5, -0.5, 4.5);
  h12 = fs->make<TH1F>("nJetHTPre300","Jets sum passing Pre-trigger H_{T} threshold 300 GeV", 5, -0.5, 4.5);
  h13 = fs->make<TH1F>("nJetHTPre350","Jets sum passing Pre-trigger H_{T} threshold 350 GeV", 5, -0.5, 4.5);
  h14 = fs->make<TH1F>("HTbx0","H_{T} at at BX = 0", 40, 0, 1000.);
  h15 = fs->make<TH1F>("HTbxm1","H_{T} at at BX = -1", 40, 0, 500.);

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

  h8->Sumw2(); h9->Sumw2(); h10->Sumw2(); h11->Sumw2(); h12->Sumw2(); h13->Sumw2();
  h9->Divide(h8); h10->Divide(h8); h11->Divide(h8); h12->Divide(h8); h13->Divide(h8);
   
  n2->Sumw2(); n3->Sumw2(); n4->Sumw2(); n5->Sumw2(); n6->Sumw2(); n7->Sumw2();
  n3->Divide(n2); n4->Divide(n2); n5->Divide(n2); n6->Divide(n2); n7->Divide(n2);
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
  
  //cout<<iEvent.id().event()<<endl;

  iEvent.getByToken(gtAlgBlkToken, gtAlgBlkHandle);
  if(gtAlgBlkHandle.isValid() && iEvent.id().event()==2874541){
    std::vector<GlobalAlgBlk>::const_iterator algBlk = gtAlgBlkHandle->begin(0);
    std::cout<<"BxVector First BX: "<<gtAlgBlkHandle->getFirstBX()<<std::endl;
    std::cout<<"BxVector Last BX: "<<gtAlgBlkHandle->getLastBX()<<std::endl;
    if(algBlk != gtAlgBlkHandle->end(0)){
      for (std::map<std::string, L1TUtmAlgorithm>::const_iterator itAlgo = algorithmMap->begin(); itAlgo != algorithmMap->end(); itAlgo++) {
  std::string algName = itAlgo->first;
        int algBit = itAlgo->second.getIndex();
        //int prescaleColumn = algBlk->getPreScColumn();
        bool initialDecision = algBlk->getAlgoDecisionInitial(algBit);
        bool intermDecision = algBlk->getAlgoDecisionInterm(algBit);
        bool decisionFinal = algBlk->getAlgoDecisionFinal(algBit);
        if(algBit==473 && initialDecision==1){
          std::cout<<"L1 Path: "<<algName<<" Bit: "<<algBit<<" Initial Decision: "<<initialDecision<<" Interm Decision: "<<intermDecision<<" Final Decision: "<<decisionFinal<<std::endl;
        }
      }
    }
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

    nJets = 0;
    HT0 = 0.; HTM1 = 0.;
    jetEtBx0=0.; jetEtBxM1=0.; jetEtaBx0=-99.; jetEtaBxM1=-99.; jetPhiBx0=-99.; jetPhiBxM1=-99.; jetEtBxM2=0.;
    
    edm::Handle<l1t::JetBxCollection> jetColl;
    iEvent.getByToken(jetBXCollectionToken_, jetColl);
    l1t::JetBxCollection jets;
    jets = (*jetColl.product()); 

    for (auto it = jets.begin(0); it!=jets.end(0); it++){
      
      nJets = nJets + 1;
      
      h1->Fill(it->pt());
      HT0=HT0+it->pt();
      
      if(it->pt() > jetEtBx0){
        jetEtBx0 = it->pt();
        jetEtaBx0 = it->eta();
        jetPhiBx0 = it->phi();
      }

    }
    h14->Fill(HT0);

    for (auto it = jets.begin(-1); it!=jets.end(-1); it++){
      
      nJets = nJets + 1;

      h2->Fill(it->pt());
      HTM1=HTM1+it->pt();
    
      if(it->pt() > jetEtBxM1){
        jetEtBxM1 = it->pt();
        jetEtaBxM1 = it->eta();
        jetPhiBxM1 = it->phi();
      }

    }
    h15->Fill(HTM1); 

    for (auto it = jets.begin(-2); it!=jets.end(-2); it++){
      
      nJets = nJets + 1;

      if(it->pt() > jetEtBxM2){
        jetEtBxM2 = it->pt(); 
      }

    }

    for (auto it = jets.begin(1); it!=jets.end(1); it++){
      nJets = nJets + 1;
    }
    for (auto it = jets.begin(2); it!=jets.end(2); it++){
      nJets = nJets + 1;
    }
   
    //HT
    for(int ibin = 0; ibin < 5; ibin++){   //loop over # of thresholds
      float j = 300. + ibin*50.;    //HT thresholds at BX=0
      std::string s = std::to_string(j);
      char const *pchar = s.c_str();
      h8->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h9->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h10->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h11->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h12->GetXaxis()->SetBinLabel(ibin+1, pchar);
      h13->GetXaxis()->SetBinLabel(ibin+1, pchar);
      if(HT0 > j){
        h8->Fill(ibin);   //fill the BX=0 denominator
        float k = 150. + ibin*50.;    //HT thresholds at BX=-1
        if(HTM1 > k){   //fill the BX=-1 numerators
          h9->Fill(ibin);
          h10->Fill(ibin);
          h11->Fill(ibin);
          h12->Fill(ibin);
          h13->Fill(ibin);
        }
      }
    }
  
    if(jetEtBx0 > 0.){
      if(jetEtBx0 > 120. && jetEtBxM1 > 60.) n8->Fill(jetEtaBxM1, jetPhiBxM1);
      for(int ibin = 0; ibin < 5; ibin++){   //loop over # of thresholds
        float j = 120. + ibin*30.;    //jet pt thresholds at BX=0
        std::string s = std::to_string(j);
        char const *pchar = s.c_str();
        n2->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n3->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n4->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n5->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n6->GetXaxis()->SetBinLabel(ibin+1, pchar);
        n7->GetXaxis()->SetBinLabel(ibin+1, pchar);

        if(jetEtBx0 > j){
          n2->Fill(ibin);   //fill the BX=0 denominator
          float k = 60. + ibin*30.;    //jet pt thresholds at BX=-1
          if(jetEtBxM1 > k){   //fill the BX=-1 numerators
            n8->Fill(jetEtaBxM1, jetPhiBxM1);
            if(jetEtBxM1 > k) n3->Fill(ibin);
            if(jetEtBxM1 > k) n4->Fill(ibin);
            if(jetEtBxM1 > k) n5->Fill(ibin);
            if(jetEtBxM1 > k) n6->Fill(ibin);
            if(jetEtBxM1 > k) n7->Fill(ibin);
          }
        }
      }
    }
   
    if(jetEtBxM1 > 60. ) std::cout<<"jetEtBx0: "<<jetEtBx0<<"\t"<<"jetEtBxM1: "<<jetEtBxM1<<std::endl;
    n1->Fill(nJets);
    h3->Fill(jetEtBx0);
    h4->Fill(jetEtBxM1);
    h5->Fill(jetEtBxM2);
    h6->Fill(jetEtaBx0);
    h7->Fill(jetEtaBxM1);
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