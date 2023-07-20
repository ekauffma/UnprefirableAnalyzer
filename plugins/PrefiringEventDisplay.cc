// -*- C++ -*-
//
// Package:    L1Trigger/UnprefirableAnalyzer
// Class:      PrefiringEventDisplay
//
/**\class PrefiringEventDisplay PrefiringEventDisplay.cc L1Trigger/UnprefirableAnalyzer/plugins/PrefiringEventDisplay.cc

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
#include <iostream>
#include <fstream>

// user include files
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
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


class PrefiringEventDisplay : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PrefiringEventDisplay(const edm::ParameterSet&);
  ~PrefiringEventDisplay() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  //edm::EDGetTokenT<BXVector<l1t::Jet>> stage2JetToken_;
  edm::EDGetTokenT<GlobalExtBlkBxCollection> UnprefirableEventToken_;
  const edm::EDGetTokenT<l1t::JetBxCollection> jetBXCollectionToken_; // l1 jets
  edm::EDGetTokenT<vector<pat::Jet>> slimmedJetsToken_; // reco jets to match to l1 jets
  edm::EDGetTokenT<vector<pat::Muon>> slimmedMuonsToken_; // needed for HLT_IsoMu20
  edm::EDGetTokenT<edm::TriggerResults> trgresultsORIGToken_; // need to require HLT_IsoMu20 in order to match reco jets

  unsigned int Run;
  unsigned int Event;
  unsigned int lumi;

  bool Flag_IsUnprefirable;
  bool Flag_FirstBunchInTrain;

  edm::ESGetToken<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd> L1TUtmTriggerMenuEventToken;
  const L1TUtmTriggerMenu* l1GtMenu;
  const std::map<std::string, L1TUtmAlgorithm>* algorithmMap;

  // define histograms
  TH2F *offlineJetsDisplay_bx0;
  TH2F *offlineJetsDisplay_bxm1;
  TH2F *offlineJetsDisplay_bx0_bxm1;
  TH2F *onlineJetsDisplay_bx0;
  TH2F *onlineJetsDisplay_bxm1;
  TH2F *onlineJetsDisplay_bx0_bxm1;

  TTree *searchForEvent;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constructors and destructor
//
PrefiringEventDisplay::PrefiringEventDisplay(const edm::ParameterSet& iConfig):
  UnprefirableEventToken_(consumes<GlobalExtBlkBxCollection>(edm::InputTag("simGtExtUnprefireable"))),
  jetBXCollectionToken_(consumes<l1t::JetBxCollection>(edm::InputTag("caloStage2Digis","Jet","RECO")))
{
  
  // tokens
  L1TUtmTriggerMenuEventToken = consumesCollector().esConsumes<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd>();
  slimmedJetsToken_ = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJetsPuppi") );
  slimmedMuonsToken_ = consumes< std::vector<pat::Muon> >(edm::InputTag("slimmedMuons") );

  trgresultsORIGToken_ = consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT") );

  edm::Service<TFileService> fs;
  
  // histograms for UnprefirableEvent
  offlineJetsDisplay_bx0 = fs->make<TH2F>("offlineJetsDisplay_bx0","Event Display",(34*5), -1.4841, 1.4841, (72*5), -3.142, 3.142);
  offlineJetsDisplay_bxm1 = fs->make<TH2F>("offlineJetsDisplay_bxm1","Event Display",(34*5), -1.4841, 1.4841, (72*5), -3.142, 3.142);
  offlineJetsDisplay_bx0_bxm1 = fs->make<TH2F>("offlineJetsDisplay_bx0_bxm1","Event Display",(34*5), -1.4841, 1.4841, (72*5), -3.142, 3.142);
  onlineJetsDisplay_bx0 = fs->make<TH2F>("onlineJetsDisplay_bx0","Event Display",(34*5), -1.4841, 1.4841, (72*5), -3.142, 3.142);
  onlineJetsDisplay_bxm1 = fs->make<TH2F>("onlineJetsDisplay_bxm1","Event Display",(34*5), -1.4841, 1.4841, (72*5), -3.142, 3.142);
  onlineJetsDisplay_bx0_bxm1 = fs->make<TH2F>("onlineJetsDisplay_bx0_bxm1","Event Display",(34*5), -1.4841, 1.4841, (72*5), -3.142, 3.142);

  searchForEvent = fs->make<TTree>("searchForEvent","searchForEvent");

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

PrefiringEventDisplay::~PrefiringEventDisplay() {
}

//
// member functions
//

// deltaR method (written using https://github.com/cms-sw/cmssw/blob/master/DataFormats/Math/interface/deltaR.h)
template <typename L1Jet, typename RecoJet>
constexpr auto deltaR(const L1Jet& l1jet, const RecoJet& recojet) -> decltype(recojet.eta()) {
  typedef decltype(recojet.eta()) Float;

  Float p1 = l1jet->phi();
  Float p2 = recojet.phi();
  Float e1 = l1jet->eta();
  Float e2 = recojet.eta();
  auto dp = std::abs(p1 - p2);
  if (dp > Float(M_PI)) dp -= Float(2 * M_PI);
  return sqrt((e1 - e2) * (e1 - e2) + dp * dp);
}

template <typename RecoJet, typename L1JetCollection>
vector<float> checkMatchBX(const RecoJet& recojet, const L1JetCollection& l1jetcollection, int bx) {

  float minDeltaR = 999.0;
  float matchedPt = 0.0;
  float matchedEta = 0.0;
  float matchedPhi = 0.0;

  vector<float> returnVec;

  for (auto it = l1jetcollection.begin(bx); it!=l1jetcollection.end(bx); it++){
    // check if deltaR value is smallest so far (finds best match if there is one)
    float currentDeltaR = deltaR(it, recojet);
    if(currentDeltaR<minDeltaR && recojet.pt()>30) {
      minDeltaR = currentDeltaR;
      matchedPt = it->pt();
      matchedEta = it->eta();
      matchedPhi = it->phi();
    }
  }

  if(minDeltaR < 0.4){
    returnVec.push_back(1.0); // match found
  }
  else{
    returnVec.push_back(0.0); // match not found
  }

  returnVec.push_back(matchedPt);
  returnVec.push_back(matchedEta);
  returnVec.push_back(matchedPhi);

  return returnVec;
}

// ------------ method called for each event  ------------
void PrefiringEventDisplay::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
 
  // get reco jets and muons
  edm::Handle< std::vector<pat::Jet> > slimmedJets;
  iEvent.getByToken(slimmedJetsToken_,slimmedJets ); 
  edm::Handle< std::vector<pat::Muon> > slimmedMuons;
  iEvent.getByToken(slimmedMuonsToken_,slimmedMuons );


  //get HLT_IsoMu20 result
  bool passHLT_IsoMu20(false); 
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsORIGToken_, trigResults);
  if( !trigResults.failedToGet() ) {
    int N_Triggers = trigResults->size();
    const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      if (trigResults.product()->accept(i_Trig)) {
        TString TrigPath =trigName.triggerName(i_Trig);
        if(TrigPath.Index("HLT_IsoMu20_v") >=0) passHLT_IsoMu20=true; 
      }
    }
  }

  if(passHLT_IsoMu20) cout<<"Event #"<<iEvent.id().event()<<", passHLT_IsoMu20 = "<<passHLT_IsoMu20<<endl;
 
  //Unprefirable
  Flag_IsUnprefirable = false;
  edm::Handle<GlobalExtBlkBxCollection> handleUnprefEventResults;
  iEvent.getByToken(UnprefirableEventToken_, handleUnprefEventResults);
  if(handleUnprefEventResults.isValid()){
    if (handleUnprefEventResults->size() != 0) {
      Flag_IsUnprefirable = handleUnprefEventResults->at(0, 0).getExternalDecision(GlobalExtBlk::maxExternalConditions-1);
    }
  }


  if(Flag_IsUnprefirable && passHLT_IsoMu20){

    cout<<"    isUnprefirable = "<<Flag_IsUnprefirable<<endl;

    edm::Handle<l1t::JetBxCollection> jetColl;
    iEvent.getByToken(jetBXCollectionToken_, jetColl);
    l1t::JetBxCollection jets;
    jets = (*jetColl.product());

    cout<<"    slimmedJets size = "<<(*slimmedJets).size()<<endl;

    Event = iEvent.id().event();

    // iterate through reco jets
    for(long unsigned int i = 0; i<(*slimmedJets).size(); i++){

      if (Event!=1687892580) continue;

      // check for match in bx0
      vector<float> matchInfo_bx0 = checkMatchBX((*slimmedJets)[i], jets, 0);
      // check for match in bxm1
      vector<float> matchInfo_bxm1 = checkMatchBX((*slimmedJets)[i], jets, -1);

      if(matchInfo_bx0.at(0)==1.0 && matchInfo_bxm1.at(0)==1.0){
        cout<<"    matched to 0 and -1"<<endl;
        cout<<"      offline eta = "<<(*slimmedJets)[i].eta()<<", phi = "<<(*slimmedJets)[i].phi()<<", pt = "<<(*slimmedJets)[i].pt()<<endl;
        cout<<"      online eta = "<<matchInfo_bx0.at(2)<<", phi = "<<matchInfo_bx0.at(3)<<", pt = "<<matchInfo_bx0.at(1)<<endl;
        offlineJetsDisplay_bx0_bxm1->Fill((*slimmedJets)[i].eta(), (*slimmedJets)[i].phi(), (*slimmedJets)[i].pt());
        onlineJetsDisplay_bx0_bxm1->Fill(matchInfo_bx0.at(2), matchInfo_bx0.at(3), matchInfo_bx0.at(1));
      }
      if(matchInfo_bx0.at(0)==1.0 && matchInfo_bxm1.at(0)==0.0){
        cout<<"    matched to 0, not -1"<<endl;
        cout<<"      offline eta = "<<(*slimmedJets)[i].eta()<<", phi = "<<(*slimmedJets)[i].phi()<<", pt = "<<(*slimmedJets)[i].pt()<<endl;
        cout<<"      online eta = "<<matchInfo_bx0.at(2)<<", phi = "<<matchInfo_bx0.at(3)<<", pt = "<<matchInfo_bx0.at(1)<<endl;
        offlineJetsDisplay_bx0->Fill((*slimmedJets)[i].eta(), (*slimmedJets)[i].phi(), (*slimmedJets)[i].pt());
        onlineJetsDisplay_bx0->Fill(matchInfo_bx0.at(2), matchInfo_bx0.at(3), matchInfo_bx0.at(1));
      }
      if(matchInfo_bx0.at(0)==0.0 && matchInfo_bxm1.at(0)==1.0){
        cout<<"    matched to -1, not 0"<<endl;
        cout<<"      offline eta = "<<(*slimmedJets)[i].eta()<<", phi = "<<(*slimmedJets)[i].phi()<<", pt = "<<(*slimmedJets)[i].pt()<<endl;
        cout<<"      online eta = "<<matchInfo_bxm1.at(2)<<", phi = "<<matchInfo_bxm1.at(3)<<", pt = "<<matchInfo_bxm1.at(1)<<endl;
        offlineJetsDisplay_bxm1->Fill((*slimmedJets)[i].eta(), (*slimmedJets)[i].phi(), (*slimmedJets)[i].pt());
        onlineJetsDisplay_bxm1->Fill(matchInfo_bxm1.at(2), matchInfo_bxm1.at(3), matchInfo_bxm1.at(1));
        Run   = iEvent.id().run();
        Event = iEvent.id().event();
        lumi  = iEvent.luminosityBlock();
        searchForEvent->Fill();
      }
    }
  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void PrefiringEventDisplay::beginJob() {
  // please remove this method if not needed
  searchForEvent->Branch("event_runNo",  &Run,   "event_runNo/I");
  searchForEvent->Branch("event_evtNo",  &Event, "event_evtNo/I");
  searchForEvent->Branch("event_lumi",   &lumi,  "event_lumi/I");
}

// ------------ method called once each job just after ending the event loop  ------------
void PrefiringEventDisplay::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PrefiringEventDisplay::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(PrefiringEventDisplay);
