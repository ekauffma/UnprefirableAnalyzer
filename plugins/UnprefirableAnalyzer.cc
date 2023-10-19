// -*- C++ -*-
//
// Package:    L1Trigger/UnprefirableAnalyzer
// Class:      UnprefirableAnalyzer
//
/**\class UnprefirableAnalyzer UnprefirableAnalyzer.cc L1Trigger/UnprefirableAnalyzer/plugins/UnprefirableAnalyzer.cc

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
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
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
#include "DataFormats/JetReco/interface/JetID.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "JetMETCorrections/Type1MET/interface/PFJetMETcorrInputProducerT.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "PhysicsTools/PatUtils/interface/PATJetCorrExtractor.h"

#include "TLorentzVector.h"
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


class UnprefirableAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit UnprefirableAnalyzer(const edm::ParameterSet&);
  ~UnprefirableAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  //edm::EDGetTokenT<BXVector<l1t::Jet>> stage2JetToken_;
  edm::EDGetTokenT<GlobalExtBlkBxCollection> UnprefirableEventToken_;
  const edm::EDGetTokenT<l1t::JetBxCollection> jetBXCollectionToken_; // l1 jets
  edm::EDGetTokenT< BXVector<GlobalAlgBlk> > gtAlgBlkToken;
  edm::Handle< BXVector<GlobalAlgBlk> > gtAlgBlkHandle;
  edm::EDGetTokenT<vector<pat::Jet>> slimmedJetsToken_; // reco jets to match to l1 jets
  edm::EDGetTokenT<vector<pat::Muon>> slimmedMuonsToken_; // needed for HLT_IsoMu24
  edm::EDGetTokenT<edm::TriggerResults> trgresultsORIGToken_; // need to require HLT_IsoMu24 in order to match reco jets
  const edm::ESGetToken<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd> l1GtMenuToken_;

  edm::ESGetToken<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd> L1TUtmTriggerMenuEventToken;
  const L1TUtmTriggerMenu* l1GtMenu;
  const std::map<std::string, L1TUtmAlgorithm>* algorithmMap;

  TTree *eventTree;
  int run_num;
  int lumi;
  int event_num;
  bool isUnprefirable;
  bool isFirstBunchInTrain;
  bool L1FinalOR_bxm1;
  vector<bool> trigger_bits;
  vector<TLorentzVector> reco_jets;
  vector<bool> reco_jetId;
  vector<TLorentzVector> match_l1_bx0;
  vector<TLorentzVector> match_l1_bxm1;
  vector<TLorentzVector> match_l1_bxm2;
  vector<TLorentzVector> match_l1_bx1;
  vector<TLorentzVector> match_l1_bx2; 
  int idx_L1_FirstBunchBeforeTrain;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constructors and destructor
//
UnprefirableAnalyzer::UnprefirableAnalyzer(const edm::ParameterSet& iConfig):
  UnprefirableEventToken_(consumes<GlobalExtBlkBxCollection>(edm::InputTag("simGtExtUnprefireable"))),
  jetBXCollectionToken_(consumes<l1t::JetBxCollection>(edm::InputTag("caloStage2Digis","Jet","RECO"))),
  gtAlgBlkToken( consumes< BXVector<GlobalAlgBlk> >(edm::InputTag("gtStage2Digis","","RECO")) ),
  l1GtMenuToken_(esConsumes<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd>())
{
  
  // tokens
  L1TUtmTriggerMenuEventToken = consumesCollector().esConsumes<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd>();
  slimmedJetsToken_ = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJetsPuppi") );
  slimmedMuonsToken_ = consumes< std::vector<pat::Muon> >(edm::InputTag("slimmedMuons") );

  trgresultsORIGToken_ = consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT") );

  edm::Service<TFileService> fs;
 
  // TTree
  eventTree = fs->make<TTree>("eventTree", "Events");

  eventTree->Branch("run_num",    &run_num,     "run_num/I");
  eventTree->Branch("lumi",   &lumi,    "lumi/I");
  eventTree->Branch("event_num",  &event_num,   "event_num/I");
  eventTree->Branch("isUnprefirable",  &isUnprefirable,   "isUnprefirable/I");
  eventTree->Branch("isFirstBunchInTrain",  &isFirstBunchInTrain,   "isFirstBunchInTrain/I");
  eventTree->Branch("L1FinalOR_bxm1", &L1FinalOR_bxm1, "L1FinalOR_bxm1/I");
  eventTree->Branch("idx_L1_FirstBunchBeforeTrain", &idx_L1_FirstBunchBeforeTrain, "idx_L1_FirstBunchBeforeTrain/I"); 
 
  eventTree->Branch("trigger_bits", "vector<bool>", &trigger_bits, 32000, 0);

  eventTree->Branch("reco_jets", "vector<TLorentzVector>", &reco_jets, 32000, 0);
  eventTree->Branch("reco_jetId", "vector<bool>", &reco_jetId, 32000, 0);

  eventTree->Branch("match_l1_bx0", "vector<TLorentzVector>", &match_l1_bx0, 32000, 0);
  eventTree->Branch("match_l1_bxm1", "vector<TLorentzVector>", &match_l1_bxm1, 32000, 0);
  eventTree->Branch("match_l1_bxm2", "vector<TLorentzVector>", &match_l1_bxm2, 32000, 0);
  eventTree->Branch("match_l1_bx1", "vector<TLorentzVector>", &match_l1_bx1, 32000, 0);
  eventTree->Branch("match_l1_bx2", "vector<TLorentzVector>", &match_l1_bx2, 32000, 0);
  
 
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

UnprefirableAnalyzer::~UnprefirableAnalyzer() {
}

//
// member functions
//

// AKpuppi jetID
template <typename RecoJet>
bool puppiJetID(const RecoJet& jet) {
  bool tmp = true;
  if (std::abs(jet.eta())<2.6) {
    tmp &= jet.neutralHadronEnergyFraction() < 0.9;
    tmp &= jet.muonEnergyFraction() < 0.8;
    tmp &= jet.chargedEmEnergyFraction() < 0.8;
    tmp &= jet.chargedMultiplicity() > 0;
    tmp &= jet.chargedHadronEnergyFraction() > 0.01;
    tmp &= (jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1;
    tmp &= jet.neutralEmEnergyFraction() < 0.9;
  }
  if (std::abs(jet.eta())>2.6 && std::abs(jet.eta())<=2.7){
    tmp &= jet.chargedEmEnergyFraction() < 0.8;
    tmp &= jet.neutralEmEnergyFraction() < 0.99;
    tmp &= jet.muonEnergyFraction() < 0.8;
    tmp &= jet.neutralHadronEnergyFraction() < 0.9;
  }
  if (std::abs(jet.eta())>2.7 && std::abs(jet.eta())<= 3.0){
    tmp &= jet.neutralEmEnergyFraction() < 0.99;
  }
  if (std::abs(jet.eta())>3.0){
    tmp &= jet.neutralEmEnergyFraction() < 0.9;
    tmp &= jet.neutralMultiplicity() > 2;
  }
  return tmp;
}


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
bool checkMatchBX(const RecoJet& recojet, const L1JetCollection& l1jetcollection, int bx, TLorentzVector* l1_jet ) {

  bool match = false;
  float temp = 99.0;
  float matchedPt = 0.0;
  float matchedEta = 0.0;
  float matchedPhi = 0.0;
  float matchedE = 0.0;

  for (auto it = l1jetcollection.begin(bx); it!=l1jetcollection.end(bx); it++){
    // check if deltaR is lower than temp
    if(deltaR(it, recojet)<temp) {
      temp = deltaR(it, recojet);
      matchedPt = it->pt();
      matchedEta = it->eta();
      matchedPhi = it->phi();
      matchedE = it->energy();
    }
  } 

  if(temp < 0.4) {
    match = true;
    l1_jet->SetPtEtaPhiE(matchedPt, matchedEta, matchedPhi, matchedE);
  }
   
  return match;
}

float phiAbs(float phi1, float phi2){
  float dp = abs(phi1-phi2);
  if (dp > M_PI){
    dp -= 2*M_PI;
  }
  return abs(dp);
}

// ------------ method called for each event  ------------
void UnprefirableAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  run_num = iEvent.id().run();
  lumi = iEvent.id().luminosityBlock();
  event_num = iEvent.id().event();

  trigger_bits.clear();
  reco_jets.clear();
  reco_jetId.clear();
  match_l1_bx0.clear();
  match_l1_bxm1.clear();
  match_l1_bxm1.clear();
  match_l1_bx1.clear();
  match_l1_bx2.clear();

  // get firstbunchbeforetrain
  edm::ESHandle<L1TUtmTriggerMenu> menu;
  menu = iSetup.getHandle(l1GtMenuToken_);
  int idx_L1_FirstBunchBeforeTrain = -1;
  for (auto const &keyval : menu->getAlgorithmMap()) {
    std::string const &name = keyval.second.getName();
    unsigned int indx = keyval.second.getIndex();
    if(name.find("L1_FirstBunchBeforeTrain")!=string::npos) idx_L1_FirstBunchBeforeTrain = indx;
  }

  // get reco jets and muons
  edm::Handle< std::vector<pat::Jet> > slimmedJets;
  iEvent.getByToken(slimmedJetsToken_,slimmedJets ); 
  edm::Handle< std::vector<pat::Muon> > slimmedMuons;
  iEvent.getByToken(slimmedMuonsToken_,slimmedMuons );

  iEvent.getByToken(gtAlgBlkToken, gtAlgBlkHandle);

  //check whether any seed fired in bx=-1
  L1FinalOR_bxm1 = false;
  for(int i =0; i <512; i++){
    trigger_bits.push_back(gtAlgBlkHandle->begin(-1)->getAlgoDecisionFinal(i));
    if (i!=idx_L1_FirstBunchBeforeTrain){
      L1FinalOR_bxm1 = L1FinalOR_bxm1 || gtAlgBlkHandle->begin(-1)->getAlgoDecisionFinal(i);      
    }
  }

  //get HLT_AK8PFJet500 result
  bool passHLT_AK8PFJet500(false); 
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsORIGToken_, trigResults);
  if( !trigResults.failedToGet() ) {
    int N_Triggers = trigResults->size();
    const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      if (trigResults.product()->accept(i_Trig)) {
        TString TrigPath =trigName.triggerName(i_Trig);
        if(TrigPath.Index("HLT_AK8PFJet500") >=0) passHLT_AK8PFJet500=true; 
      }
    }
  }

  auto menuRcd = iSetup.get<L1TUtmTriggerMenuRcd>();
  l1GtMenu = &menuRcd.get(L1TUtmTriggerMenuEventToken);
  algorithmMap = &(l1GtMenu->getAlgorithmMap());
  
  //FirstBunchInTrain
  isFirstBunchInTrain = false;
  if(gtAlgBlkHandle.isValid()){
    std::vector<GlobalAlgBlk>::const_iterator algBlk = gtAlgBlkHandle->begin(0);
    if(algBlk != gtAlgBlkHandle->end(0)){
      for (std::map<std::string, L1TUtmAlgorithm>::const_iterator itAlgo = algorithmMap->begin(); itAlgo != algorithmMap->end(); itAlgo++) {
        std::string algName = itAlgo->first;
        int algBit = itAlgo->second.getIndex();
        bool initialDecision = algBlk->getAlgoDecisionInitial(algBit);
        if(algBit==473 && initialDecision==1){
          isFirstBunchInTrain = true;
        }
      }
    }
  }

  //Unprefirable
  isUnprefirable = false;
  edm::Handle<GlobalExtBlkBxCollection> handleUnprefEventResults;
  iEvent.getByToken(UnprefirableEventToken_, handleUnprefEventResults);
  if(handleUnprefEventResults.isValid()){
    if (handleUnprefEventResults->size() != 0) {
      isUnprefirable = handleUnprefEventResults->at(0, 0).getExternalDecision(GlobalExtBlk::maxExternalConditions-1);
    }
  }

  if(passHLT_AK8PFJet500){
    
    edm::Handle<l1t::JetBxCollection> jetColl;
    iEvent.getByToken(jetBXCollectionToken_, jetColl);
    l1t::JetBxCollection jets;
    jets = (*jetColl.product());

    // iterate through reco jets
    for(long unsigned int i = 0; i<(*slimmedJets).size(); i++){
     
      reco_jetId.push_back(puppiJetID((*slimmedJets)[i]));
      
      TLorentzVector reco_jet;
      reco_jet.SetPtEtaPhiE((*slimmedJets)[i].pt(), (*slimmedJets)[i].eta(), (*slimmedJets)[i].phi(), (*slimmedJets)[i].energy());
      reco_jets.push_back(reco_jet); 
 
      TLorentzVector l1_jet_bx0; // keep track of matched l1 pt
      TLorentzVector l1_jet_bxm1;
      TLorentzVector l1_jet_bxm2;
      TLorentzVector l1_jet_bx1;
      TLorentzVector l1_jet_bx2;

      // match jets
      bool match_bx0 = checkMatchBX((*slimmedJets)[i], jets, 0, &l1_jet_bx0);
      bool match_bxm1 = checkMatchBX((*slimmedJets)[i], jets, -1, &l1_jet_bxm1);
      bool match_bxm2 = checkMatchBX((*slimmedJets)[i], jets, -2, &l1_jet_bxm2);
      bool match_bx1 = checkMatchBX((*slimmedJets)[i], jets, 1, &l1_jet_bx1);
      bool match_bx2 = checkMatchBX((*slimmedJets)[i], jets, 2, &l1_jet_bx2);
      
      if(match_bx0){
        match_l1_bx0.push_back(l1_jet_bx0);
      }
      if(match_bxm1){
        match_l1_bxm1.push_back(l1_jet_bxm1);
      }
      if(match_bxm2){
        match_l1_bxm2.push_back(l1_jet_bxm2);
      }
      if(match_bx1){
        match_l1_bx1.push_back(l1_jet_bx1);
      }
      if(match_bx2){
        match_l1_bx2.push_back(l1_jet_bx2);
      }
    }
    eventTree->Fill();
  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void UnprefirableAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void UnprefirableAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void UnprefirableAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(UnprefirableAnalyzer);
