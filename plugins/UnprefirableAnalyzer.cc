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

  bool Flag_IsUnprefirable;
  bool Flag_FirstBunchInTrain;

  edm::ESGetToken<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd> L1TUtmTriggerMenuEventToken;
  const L1TUtmTriggerMenu* l1GtMenu;
  const std::map<std::string, L1TUtmAlgorithm>* algorithmMap;

  // define histograms
  TH1F *h_jetet_bx0;
  TH1F *h_jetet_bxm1;
  TH1F *h_jetet_bx0_bxm1;
  TH1F *h_jetet_bxm2;
  TH1F *h_jetet_bx1;
  TH1F *h_jetet_bx2;

  TH1F *h_jeteta_bx0;
  TH1F *h_jeteta_bxm1;
  TH1F *h_jeteta_bx0_bxm1;
  TH1F *h_jeteta_bxm2;
  TH1F *h_jeteta_bx1;
  TH1F *h_jeteta_bx2;

  TH1F *h_jeteres_bx0;
  TH1F *h_jeteres_bxm1;
  TH1F *h_jeteres_bxm2;
  TH1F *h_jeteres_bx1;
  TH1F *h_jeteres_bx2;

  TH1F *h_lowrespt_bx0;
  TH1F *h_lowrespt_bxm1;
  TH1F *h_lowrespt_bxm2;
  TH1F *h_lowrespt_bx1;
  TH1F *h_lowrespt_bx2;

  TH2F *h_jetetaphi_bx0;
  TH2F *h_jetetaphi_bxm1;
  TH2F *h_jetetaphi_bxm2;
  TH2F *h_jetetaphi_bx1;
  TH2F *h_jetetaphi_bx2;

  TH2F *h_jetetaphi_bx0_on;
  TH2F *h_jetetaphi_bxm1_on;
  TH2F *h_jetetaphi_bxm2_on;
  TH2F *h_jetetaphi_bx1_on;
  TH2F *h_jetetaphi_bx2_on;

  TH2F *h_jetpteta_bx0;
  TH2F *h_jetpteta_bxm1;
  TH2F *h_jetpteta_bxm2;
  TH2F *h_jetpteta_bx1;
  TH2F *h_jetpteta_bx2;
  TH2F *h_jetpteta_bxm1_nobx0;

  TH1F *h_jeteta_lowpt_bx0;
  TH1F *h_jeteta_medpt_bx0;
  TH1F *h_jeteta_highpt_bx0;
  TH1F *h_jeteta_lowpt_bxm1;
  TH1F *h_jeteta_medpt_bxm1;
  TH1F *h_jeteta_highpt_bxm1;
  TH1F *h_jeteta_lowpt_bx0_bxm1;
  TH1F *h_jeteta_medpt_bx0_bxm1;
  TH1F *h_jeteta_highpt_bx0_bxm1;

  TH1I *nbx;
  TH1I *h_nevt;

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
  gtAlgBlkToken( consumes< BXVector<GlobalAlgBlk> >(edm::InputTag("gtStage2Digis","","RECO")) )
{
  
  // tokens
  L1TUtmTriggerMenuEventToken = consumesCollector().esConsumes<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd>();
  slimmedJetsToken_ = consumes< std::vector<pat::Jet> >(edm::InputTag("slimmedJetsPuppi") );
  slimmedMuonsToken_ = consumes< std::vector<pat::Muon> >(edm::InputTag("slimmedMuons") );

  trgresultsORIGToken_ = consumes<edm::TriggerResults>( edm::InputTag("TriggerResults::HLT") );

  edm::Service<TFileService> fs;
  
  // histograms
  h_jetet_bx0 = fs->make<TH1F>("JetEt_bx0","Jet E_{T} at BX = 0", 40, 0, 250);
  h_jetet_bxm1 = fs->make<TH1F>("JetEt_bxm1","Jet E_{T} at BX = -1", 40, 0, 250);
  h_jetet_bx0_bxm1 = fs->make<TH1F>("JetEt_bx0_bxm1","Jet E_{T} at BX = 0 or -1", 40, 0, 250);  
  h_jetet_bxm2 = fs->make<TH1F>("JetEt_bxm2","Jet E_{T} at BX = -2", 40, 0, 250);
  h_jetet_bx1 = fs->make<TH1F>("JetEt_bx1","Jet E_{T} at BX = 1", 40, 0, 250);
  h_jetet_bx2 = fs->make<TH1F>("JetEt_bx2","Jet E_{T} at BX = 2", 40, 0, 250);

  h_jeteta_bx0 = fs->make<TH1F>("JetEta_bx0","Jet #eta at BX = 0", 40, -5, 5);
  h_jeteta_bxm1 = fs->make<TH1F>("JetEta_bxm1","Jet #eta at BX = -1", 40, -5, 5);
  h_jeteta_bx0_bxm1 = fs->make<TH1F>("JetEta_bx0_bxm1","Jet #eta at BX = 0 or -1", 40, -5, 5);
  h_jeteta_bxm2 = fs->make<TH1F>("JetEta_bxm2","Jet #eta at BX = -2", 40, -5, 5);
  h_jeteta_bx1 = fs->make<TH1F>("JetEta_bx1","Jet #eta at BX = 1", 40, -5, 5);
  h_jeteta_bx2 = fs->make<TH1F>("JetEta_bx2","Jet #eta at BX = 2", 40, -5, 5);

  h_jeteres_bx0 = fs->make<TH1F>("JetERes_bx0","Jet Energy Resolution at BX = 0", 100, 0, 10);
  h_jeteres_bxm1 = fs->make<TH1F>("JetERes_bxm1","Jet Energy Resolution at BX = -1", 100, 0, 10);
  h_jeteres_bxm2 = fs->make<TH1F>("JetERes_bxm2","Jet Energy Resolution at BX = -2", 100, 0, 10);
  h_jeteres_bx1 = fs->make<TH1F>("JetERes_bx1","Jet Energy Resolution at BX = 1", 100, 0, 10);
  h_jeteres_bx2 = fs->make<TH1F>("JetERes_bx2","Jet Energy Resolution at BX = 2", 100, 0, 10);

  h_lowrespt_bx0 = fs->make<TH1F>("JetPt_lowres_bx0", "Online Jet p_T for Low Resolution (<0.5) Jets", 40, 30, 250);
  h_lowrespt_bxm1 = fs->make<TH1F>("JetPt_lowres_bxm1", "Online Jet p_T for Low Resolution (<0.5) Jets", 40, 30, 250);
  h_lowrespt_bxm2 = fs->make<TH1F>("JetPt_lowres_bxm2", "Online Jet p_T for Low Resolution (<0.5) Jets", 40, 30, 250);
  h_lowrespt_bx1 = fs->make<TH1F>("JetPt_lowres_bx1", "Online Jet p_T for Low Resolution (<0.5) Jets", 40, 30, 250);
  h_lowrespt_bx2 = fs->make<TH1F>("JetPt_lowres_bx2", "Online Jet p_T for Low Resolution (<0.5) Jets", 40, 30, 250);

  h_jetetaphi_bx0 = fs->make<TH2F>("JetEtaPhi_bx0","#eta vs #phi of jets with p_T>30 GeV (BX=0)",40, -5, 5, 40, -M_PI, M_PI);
  h_jetetaphi_bxm1 = fs->make<TH2F>("JetEtaPhi_bxm1","#eta vs #phi of jets with p_T>30 GeV (BX=-1)",40, -5, 5, 40, -M_PI, M_PI);
  h_jetetaphi_bxm2 = fs->make<TH2F>("JetEtaPhi_bxm2","#eta vs #phi of jets with p_T>30 GeV (BX=-2)",40, -5, 5, 40, -M_PI, M_PI);
  h_jetetaphi_bx1 = fs->make<TH2F>("JetEtaPhi_bx1","#eta vs #phi of jets with p_T>30 GeV (BX=1)",40, -5, 5, 40, -M_PI, M_PI);
  h_jetetaphi_bx2 = fs->make<TH2F>("JetEtaPhi_bx2","#eta vs #phi of jets with p_T>30 GeV (BX=2)",40, -5, 5, 40, -M_PI, M_PI);
  
  h_jetetaphi_bx0_on = fs->make<TH2F>("JetEtaPhi_bx0_online","#eta vs #phi of offline jets with p_T>30 GeV (BX=0)",40, -5, 5, 40, -M_PI, M_PI);
  h_jetetaphi_bxm1_on = fs->make<TH2F>("JetEtaPhi_bxm1_online","#eta vs #phi of offline jets with p_T>30 GeV (BX=-1)",40, -5, 5, 40, -M_PI, M_PI);
  h_jetetaphi_bxm2_on = fs->make<TH2F>("JetEtaPhi_bxm2_online","#eta vs #phi of offline jets with p_T>30 GeV (BX=-2)",40, -5, 5, 40, -M_PI, M_PI);
  h_jetetaphi_bx1_on = fs->make<TH2F>("JetEtaPhi_bx1_online","#eta vs #phi of offline jets with p_T>30 GeV (BX=1)",40, -5, 5, 40, -M_PI, M_PI);
  h_jetetaphi_bx2_on = fs->make<TH2F>("JetEtaPhi_bx2_online","#eta vs #phi of offline jets with p_T>30 GeV (BX=2)",40, -5, 5, 40, -M_PI, M_PI);

  h_jetpteta_bx0 = fs->make<TH2F>("JetPtEta_bx0", "Jet p_T vs #eta", 40, 30, 250, 40, -5, 5);
  h_jetpteta_bxm1 = fs->make<TH2F>("JetPtEta_bxm1", "Jet p_T vs #eta", 40, 30, 250, 40, -5, 5);
  h_jetpteta_bxm2 = fs->make<TH2F>("JetPtEta_bxm2", "Jet p_T vs #eta", 40, 30, 250, 40, -5, 5);
  h_jetpteta_bx1 = fs->make<TH2F>("JetPtEta_bx1", "Jet p_T vs #eta", 40, 30, 250, 40, -5, 5);
  h_jetpteta_bx2 = fs->make<TH2F>("JetPtEta_bx2", "Jet p_T vs #eta", 40, 30, 250, 40, -5, 5);
  h_jetpteta_bxm1_nobx0 = fs->make<TH2F>("JetPtEta_bxm1_nobx0", "Jet p_T vs #eta", 40, 30, 250, 40, -5, 5);

  h_jeteta_lowpt_bx0 = fs->make<TH1F>("JetEta_lowpt_bx0", "Jet #eta for L1 pT<=15 GeV at BX=0", 40, -5, 5);
  h_jeteta_medpt_bx0 = fs->make<TH1F>("JetEta_medpt_bx0", "Jet #eta for L1 15 GeV<pT<=30 GeV at BX=0", 40, -5, 5);
  h_jeteta_highpt_bx0 = fs->make<TH1F>("JetEta_high_bx0", "Jet #eta for L1 pT>30 GeV at BX=0", 40, -5, 5);
  h_jeteta_lowpt_bxm1 = fs->make<TH1F>("JetEta_lowpt_bxm1", "Jet #eta for L1 pT<=15 GeV at BX=-1", 40, -5, 5);
  h_jeteta_medpt_bxm1 = fs->make<TH1F>("JetEta_medpt_bxm1", "Jet #eta for L1 15 GeV<pT<=30 GeV at BX=-1", 40, -5, 5);
  h_jeteta_highpt_bxm1 = fs->make<TH1F>("JetEta_high_bxm1", "Jet #eta for L1 pT>30 GeV at BX=-1", 40, -5, 5);
  h_jeteta_lowpt_bx0_bxm1 = fs->make<TH1F>("JetEta_lowpt_bx0_bxm1", "Jet #eta for L1 pT<=15 GeV at BX=0 or BX=-1", 40, -5, 5);
  h_jeteta_medpt_bx0_bxm1 = fs->make<TH1F>("JetEta_medpt_bx0_bxm1", "Jet #eta for L1 15 GeV<pT<=30 GeV at BX=0 or BX=-1", 40, -5, 5);
  h_jeteta_highpt_bx0_bxm1 = fs->make<TH1F>("JetEta_high_bx0_bxm1", "Jet #eta for L1 pT>30 GeV at BX=0 or BX=-1", 40, -5, 5);

  nbx = fs->make<TH1I>("nJets_bx", "Number of jets per bx",5,-2.0,2.0);

  h_nevt = fs->make<TH1I>("nEvt_category", "Number of events passing each flag",10,0.0,10.0);

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
bool checkMatchBX(const RecoJet& recojet, const L1JetCollection& l1jetcollection, int bx, float* l1_pt, 
                  TH1I* nbx, TH1F* h_jetet, TH1F* h_jeteta, TH2F* h_jetetaphi, TH2F* h_jetetaphi_on, TH2F* h_jetpteta, TH1F* h_jeteres, TH1F* h_lowrespt) {

  bool match = false;
  float temp = 99.0;
  float matchedPt = 0.0;
  float matchedEta = 0.0;
  float matchedPhi = 0.0;

  for (auto it = l1jetcollection.begin(bx); it!=l1jetcollection.end(bx); it++){
    // check if deltaR is lower than temp
    if(deltaR(it, recojet)<temp) {
      temp = deltaR(it, recojet);
      matchedPt = it->pt();
      matchedEta = it->eta();
      matchedPhi = it->phi();
    }
  } 

  if(temp < 0.4) {
    match = true;
    if(recojet.pt()>30){
      nbx->Fill(bx);
      h_jetet->Fill(recojet.pt());
      h_jeteta->Fill(recojet.eta());
      h_jetetaphi->Fill(recojet.eta(), recojet.phi());
      h_jetetaphi_on->Fill(matchedEta,matchedPhi);
      h_jetpteta->Fill(recojet.pt(), recojet.eta());
      h_jeteres->Fill(matchedPt/recojet.pt());
      if(matchedPt/recojet.pt()<0.5){
        h_lowrespt->Fill(recojet.pt());
      }
    }
    
    // assign l1 jet pt
    *l1_pt = matchedPt;
  }
  return match;
}

// ------------ method called for each event  ------------
void UnprefirableAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
 
  // get reco jets and muons
  edm::Handle< std::vector<pat::Jet> > slimmedJets;
  iEvent.getByToken(slimmedJetsToken_,slimmedJets ); 
  edm::Handle< std::vector<pat::Muon> > slimmedMuons;
  iEvent.getByToken(slimmedMuonsToken_,slimmedMuons );


  //get HLT_IsoMu24 result
  bool passHLT_IsoMu24(false); 
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsORIGToken_, trigResults);
  if( !trigResults.failedToGet() ) {
    int N_Triggers = trigResults->size();
    const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      if (trigResults.product()->accept(i_Trig)) {
        TString TrigPath =trigName.triggerName(i_Trig);
        if(TrigPath.Index("HLT_IsoMu24_v") >=0) passHLT_IsoMu24=true; 
      }
    }
  }

  if(passHLT_IsoMu24) cout<<"Event #"<<iEvent.id().event()<<", passHLT_IsoMu24 = "<<passHLT_IsoMu24<<endl;
 
  auto menuRcd = iSetup.get<L1TUtmTriggerMenuRcd>();
  l1GtMenu = &menuRcd.get(L1TUtmTriggerMenuEventToken);
  algorithmMap = &(l1GtMenu->getAlgorithmMap());
  
  //FirstBunchInTrain
  Flag_FirstBunchInTrain = false;
  iEvent.getByToken(gtAlgBlkToken, gtAlgBlkHandle);
  if(gtAlgBlkHandle.isValid()){
    std::vector<GlobalAlgBlk>::const_iterator algBlk = gtAlgBlkHandle->begin(0);
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

  //Unprefirable
  Flag_IsUnprefirable = false;
  edm::Handle<GlobalExtBlkBxCollection> handleUnprefEventResults;
  iEvent.getByToken(UnprefirableEventToken_, handleUnprefEventResults);
  if(handleUnprefEventResults.isValid()){
    if (handleUnprefEventResults->size() != 0) {
      Flag_IsUnprefirable = handleUnprefEventResults->at(0, 0).getExternalDecision(GlobalExtBlk::maxExternalConditions-1);
    }
  }

  // for checking which events pass which flag 
  if (!Flag_IsUnprefirable && !Flag_FirstBunchInTrain && !passHLT_IsoMu24 ) h_nevt->Fill(0);
  if (Flag_IsUnprefirable && !Flag_FirstBunchInTrain && !passHLT_IsoMu24 ) h_nevt->Fill(1);
  if (!Flag_IsUnprefirable && Flag_FirstBunchInTrain && !passHLT_IsoMu24 ) h_nevt->Fill(2); 
  if (!Flag_IsUnprefirable && !Flag_FirstBunchInTrain && passHLT_IsoMu24 ) h_nevt->Fill(3);
  if (!Flag_IsUnprefirable && Flag_FirstBunchInTrain && passHLT_IsoMu24 ) h_nevt->Fill(4);
  if (Flag_IsUnprefirable && !Flag_FirstBunchInTrain && passHLT_IsoMu24 ) h_nevt->Fill(5);
  if (Flag_IsUnprefirable && Flag_FirstBunchInTrain && !passHLT_IsoMu24 ) h_nevt->Fill(6);
  if (Flag_IsUnprefirable && Flag_FirstBunchInTrain && passHLT_IsoMu24 ) h_nevt->Fill(7);


  if((Flag_IsUnprefirable || Flag_FirstBunchInTrain) && passHLT_IsoMu24){

    cout<<"    isUnprefirable = "<<Flag_IsUnprefirable<<endl;
    cout<<"    FirstBunchInTrain = "<<Flag_FirstBunchInTrain<<endl;

    edm::Handle<l1t::JetBxCollection> jetColl;
    iEvent.getByToken(jetBXCollectionToken_, jetColl);
    l1t::JetBxCollection jets;
    jets = (*jetColl.product());

    cout<<"    slimmedJets size = "<<(*slimmedJets).size()<<endl;

    // iterate through reco jets
    for(long unsigned int i = 0; i<(*slimmedJets).size(); i++){

      float l1_pt;

      // match jets and fill histograms
      bool match_bx0 = checkMatchBX((*slimmedJets)[i], 
                                    jets, 
                                    0, 
                                    &l1_pt,
                                    nbx, 
                                    h_jetet_bx0, 
                                    h_jeteta_bx0, 
                                    h_jetetaphi_bx0, 
                                    h_jetetaphi_bx0_on, 
                                    h_jetpteta_bx0, 
                                    h_jeteres_bx0, 
                                    h_lowrespt_bx0);
      bool match_bxm1 = checkMatchBX((*slimmedJets)[i],
                                     jets,
                                     -1,
                                     &l1_pt,
                                     nbx,
                                     h_jetet_bxm1,
                                     h_jeteta_bxm1,
                                     h_jetetaphi_bxm1,
                                     h_jetetaphi_bxm1_on,
                                     h_jetpteta_bxm1,
                                     h_jeteres_bxm1,
                                     h_lowrespt_bxm1);
      bool match_bxm2 = checkMatchBX((*slimmedJets)[i],
                                     jets,
                                     -2,
                                     &l1_pt,
                                     nbx,
                                     h_jetet_bxm2,
                                     h_jeteta_bxm2,
                                     h_jetetaphi_bxm2,
                                     h_jetetaphi_bxm2_on,
                                     h_jetpteta_bxm2,
                                     h_jeteres_bxm2,
                                     h_lowrespt_bxm2);
      bool match_bx1 = checkMatchBX((*slimmedJets)[i],
                                    jets,
                                    1,
                                    &l1_pt,
                                    nbx,
                                    h_jetet_bx1,
                                    h_jeteta_bx1,
                                    h_jetetaphi_bx1,
                                    h_jetetaphi_bx1_on,
                                    h_jetpteta_bx1,
                                    h_jeteres_bx1,
                                    h_lowrespt_bx1);
      bool match_bx2 = checkMatchBX((*slimmedJets)[i],
                                    jets,
                                    2,
                                    &l1_pt,
                                    nbx,
                                    h_jetet_bx2,
                                    h_jeteta_bx2,
                                    h_jetetaphi_bx2,
                                    h_jetetaphi_bx2_on,
                                    h_jetpteta_bx2,
                                    h_jeteres_bx2,
                                    h_lowrespt_bx2);

      if((match_bx0 || match_bxm1) && (*slimmedJets)[i].pt()>30){
        h_jetet_bx0_bxm1->Fill((*slimmedJets)[i].pt());
        h_jeteta_bx0_bxm1->Fill((*slimmedJets)[i].eta());
      }

      if((match_bxm1 && !match_bx0) && (*slimmedJets)[i].pt()>30){
        h_jetpteta_bxm1_nobx0->Fill((*slimmedJets)[i].pt(), (*slimmedJets)[i].eta());
      }

      // pT separated hists for bx=0
      if(match_bx0){
        if(l1_pt<=15){
          h_jeteta_lowpt_bx0->Fill((*slimmedJets)[i].eta());
        }
        if((l1_pt>15) && (l1_pt<=30)){
          h_jeteta_medpt_bx0->Fill((*slimmedJets)[i].eta());
        }
        if(l1_pt>30){
          h_jeteta_highpt_bx0->Fill((*slimmedJets)[i].eta());
        }
      }

      // pT separated hists for bx=-1
      if(match_bxm1){
        if(l1_pt<=15){
          h_jeteta_lowpt_bxm1->Fill((*slimmedJets)[i].eta());
        }
        if((l1_pt>15) && (l1_pt<=30)){
          h_jeteta_medpt_bxm1->Fill((*slimmedJets)[i].eta());
        }
        if(l1_pt>30){
          h_jeteta_highpt_bxm1->Fill((*slimmedJets)[i].eta());
        }
      }

      // pT separated hists for bx=0 or bx=-1
      if(match_bx0 || match_bxm1){
        if(l1_pt<=15){
          h_jeteta_lowpt_bx0_bxm1->Fill((*slimmedJets)[i].eta());
        }
        if((l1_pt>15) && (l1_pt<=30)){
          h_jeteta_medpt_bx0_bxm1->Fill((*slimmedJets)[i].eta());
        }
        if(l1_pt>30){
          h_jeteta_highpt_bx0_bxm1->Fill((*slimmedJets)[i].eta());
        }
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
