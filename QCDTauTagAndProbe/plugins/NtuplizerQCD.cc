// -*- C++ -*-
//
// Package:    QCDTagAndProbe/QCDTagAndProbe
// Class:      NtuplizerQCD
// 
/**\class NtuplizerQCD NtuplizerQCD.cc QCDTagAndProbe/QCDTagAndProbe/plugins/NtuplizerQCD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maximilian Burkart
//         Created:  Tue, 03 Jul 2018 12:50:23 GMT
//
//

#ifndef NTUPLIZERQCD_H
#define NTUPLIZERQCD_H
// system include files
#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>


// user include files
#include <TNtuple.h>
#include <TString.h>

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>

#include <FWCore/Framework/interface/Event.h>

#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <FWCore/ServiceRegistry/interface/Service.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <FWCore/Common/interface/TriggerNames.h> 
#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>
#include <HLTrigger/HLTcore/interface/HLTConfigProvider.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include <DataFormats/PatCandidates/interface/GenericParticle.h>

#include "tParameterSet.h"

#include <CommonTools/UtilAlgos/interface/TFileService.h>



//Set this variable to decide the number of triggers that you want to check simultaneously
#define NUMBER_OF_MAXIMUM_TRIGGERS 64
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class NtuplizerQCD : public edm::EDAnalyzer {
   public:
      explicit NtuplizerQCD(const edm::ParameterSet&);
      virtual ~NtuplizerQCD();

   private:
      virtual void beginJob();
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob();
      virtual void endRun(edm::Run const&, edm::EventSetup const&);

      void Initialize();
      bool hasFilters(const pat::TriggerObjectStandAlone&, const std::vector<std::string>&);
      bool passesJetID(const pat::JetRef&, const int);
      int GenIndex(const pat::TauRef&, const edm::View<pat::GenericParticle>*);
      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::TauRefVector> tauTag_;
      edm::EDGetTokenT<pat::JetRefVector> jetTag_;
      edm::EDGetTokenT<edm::View<pat::GenericParticle>> genParticlesTag_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<std::vector<reco::Vertex>> VtxTag_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puTag_;
      edm::EDGetTokenT<pat::METCollection> metTag_;
      edm::InputTag processName_;
      std::string treeName_;
      bool useGenMatch_;
      bool useHLTMatch_;
      bool isMC_;
      std::string filterPath_;
      int year_;

      TTree* tree_;
      TTree* triggerNamesTree_;
      
      // ---------------------------------------
      // variables to be filled in output tree
      ULong64_t indexevents_;
      Int_t runNumber_;
      Int_t lumi_;
      Int_t PS_column_;

      float MC_weight_;

      unsigned long tauTriggerBits_;
      unsigned long jetTagTriggerBits_;
      float tauPt_;
      float tauEta_;
      float tauPhi_;
      int tauCharge_;
      int tauDM_;
      int tau_genindex_;
      float tauTrkPt_;

      float jetTagPt_;
      float jetTagEta_;
      float jetTagPhi_;
      float jetProbePt_;
      float jetProbeEta_;
      float jetProbePhi_;
      int jetProbePartonFlavour_;
      float jetProbeNeutHadFrac_;
      float jetProbeNeutEMFrac_;
      int jetProbeNumConsts_;
      float jetProbeMuonFrac_;
      float jetProbeChargedHadFrac_;
      int jetProbeChargedMult_;
      float jetProbeChargedEMFrac_;
      int jetProbeN60_;
      int jetProbeN90_;
      float jetProbeQgLikelihood_;
      int nJets_;

      float metPt_;
      float metPhi_;

      bool tauDecayModeFinding_;
      bool tauDecayModeFindingNewDMs_;

      bool tauByLooseCombinedIsolationDeltaBetaCorr3Hits_;
      bool tauByMediumCombinedIsolationDeltaBetaCorr3Hits_;
      bool tauByTightCombinedIsolationDeltaBetaCorr3Hits_;

      float tauByIsolationMVArun2017v2DBoldDMwLTraw2017_;
      bool tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_;
      bool tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_;
      bool tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_;
      bool tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_;
      bool tauByTightIsolationMVArun2017v2DBoldDMwLT2017_;
      bool tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_;
      bool tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_;

      float tauByIsolationMVArun2017v2DBnewDMwLTraw2017_;
      bool tauByVVLooseIsolationMVArun2017v2DBnewDMwLT2017_;
      bool tauByVLooseIsolationMVArun2017v2DBnewDMwLT2017_;
      bool tauByLooseIsolationMVArun2017v2DBnewDMwLT2017_;
      bool tauByMediumIsolationMVArun2017v2DBnewDMwLT2017_;
      bool tauByTightIsolationMVArun2017v2DBnewDMwLT2017_;
      bool tauByVTightIsolationMVArun2017v2DBnewDMwLT2017_;
      bool tauByVVTightIsolationMVArun2017v2DBnewDMwLT2017_;

      float tauByDeepTau2017v2p1VSjetraw_;
      bool tauByVVVLooseDeepTau2017v2p1VSjet_;
      bool tauByVVLooseDeepTau2017v2p1VSjet_;
      bool tauByVLooseDeepTau2017v2p1VSjet_;
      bool tauByLooseDeepTau2017v2p1VSjet_;
      bool tauByMediumDeepTau2017v2p1VSjet_;
      bool tauByTightDeepTau2017v2p1VSjet_;
      bool tauByVTightDeepTau2017v2p1VSjet_;
      bool tauByVVTightDeepTau2017v2p1VSjet_;

      float tauByDeepTau2017v2p1VSeraw_;
      bool tauByVVVLooseDeepTau2017v2p1VSe_;
      bool tauByVVLooseDeepTau2017v2p1VSe_;
      bool tauByVLooseDeepTau2017v2p1VSe_;
      bool tauByLooseDeepTau2017v2p1VSe_;
      bool tauByMediumDeepTau2017v2p1VSe_;
      bool tauByTightDeepTau2017v2p1VSe_;
      bool tauByVTightDeepTau2017v2p1VSe_;
      bool tauByVVTightDeepTau2017v2p1VSe_;
        
      float tauByDeepTau2017v2p1VSmuraw_;
      bool tauByVLooseDeepTau2017v2p1VSmu_;
      bool tauByLooseDeepTau2017v2p1VSmu_;
      bool tauByMediumDeepTau2017v2p1VSmu_;
      bool tauByTightDeepTau2017v2p1VSmu_;

      bool tauAgainstMuonLoose3_;
      bool tauAgainstMuonTight3_;
      bool tauAgainstElectronVLooseMVA6_;
      bool tauAgainstElectronLooseMVA6_;
      bool tauAgainstElectronMediumMVA6_;
      bool tauAgainstElectronTightMVA6_;
      bool tauAgainstElectronVTightMVA6_;

      int Nvtx_;
      float nTruePU_;

      bool isTagHLTmatched_;
      bool isProbeHLTmatched_;


      //!Contains the parameters
      tVParameterSet parameters_Tag_;
      tVParameterSet parameters_;
      
      //! Maximum
      std::bitset<NUMBER_OF_MAXIMUM_TRIGGERS> tauTriggerBitSet_;
      std::bitset<NUMBER_OF_MAXIMUM_TRIGGERS> jetTagTriggerBitSet_;

      HLTConfigProvider hltConfig_;

      UInt_t lastFilter_;
      std::vector <std::string> triggerModules_;
      TString filterLabel_;
      TTree* filterLabelsTree_;
      unsigned int lastFilterInd_;

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
NtuplizerQCD::NtuplizerQCD(const edm::ParameterSet& iConfig) :
    tauTag_ (consumes<pat::TauRefVector>  (iConfig.getParameter<edm::InputTag>("taus"))),
    jetTag_ (consumes<pat::JetRefVector>  (iConfig.getParameter<edm::InputTag>("jets"))),
    genParticlesTag_ (consumes<edm::View<pat::GenericParticle>>  (iConfig.getParameter<edm::InputTag>("genParticles"))),
    triggerObjects_ (consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("triggerSet"))),
    triggerBits_ (consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggerResultsLabel"))),
    VtxTag_ (consumes<std::vector<reco::Vertex>> (iConfig.getParameter<edm::InputTag>("Vertices"))),
    puTag_ (consumes<std::vector<PileupSummaryInfo>> (iConfig.getParameter<edm::InputTag>("puInfo"))),
    metTag_ (consumes<pat::METCollection> (iConfig.getParameter<edm::InputTag>("met"))),
    processName_ (iConfig.getParameter<edm::InputTag>("triggerResultsLabel")),
    treeName_ (iConfig.getParameter<std::string>("treeName")),
    useHLTMatch_ (iConfig.getParameter<bool>("useHLTMatch")),
    isMC_ (iConfig.getParameter<bool>("isMC")),
    filterPath_ (iConfig.getParameter<std::string>("filterPath")),
    year_ (iConfig.getParameter<int>("year"))
{
   TString triggerName;
   edm::Service<TFileService> fs;
   triggerNamesTree_ = fs->make<TTree>("triggerNames", "triggerNames");
   triggerNamesTree_->Branch("triggerNames", &triggerName);

   filterLabelsTree_ = fs->make<TTree>("filterLabels", "filterLabels");
   filterLabelsTree_->Branch("filterLabels", &filterLabel_);

   //Building the trigger arrays
   const std::vector<edm::ParameterSet>& HLTList = iConfig.getParameter <std::vector<edm::ParameterSet> > ("triggerListProbe");
   for (const edm::ParameterSet& parameterSet : HLTList) {
       tParameterSet pSet;
       pSet.hltPath = parameterSet.getParameter<std::string>("HLT");
       triggerName = pSet.hltPath;
       pSet.hltFilters1 = parameterSet.getParameter<std::vector<std::string> >("path1");
       pSet.hltFilters2 = parameterSet.getParameter<std::vector<std::string> >("path2");
       pSet.leg1 = parameterSet.getParameter<int>("leg1");
       pSet.leg2 = parameterSet.getParameter<int>("leg2");
       parameters_.push_back(pSet);

       triggerNamesTree_->Fill();
   }

   const std::vector<edm::ParameterSet>& HLTList_Tag = iConfig.getParameter <std::vector<edm::ParameterSet> > ("triggerListTag");
   for (const edm::ParameterSet& parameterSet : HLTList_Tag) {
       tParameterSet pSet;
       pSet.hltPath = parameterSet.getParameter<std::string>("HLT");
       pSet.hltFilters1 = parameterSet.getParameter<std::vector<std::string>>("path1");
       pSet.hltFilters2 = parameterSet.getParameter<std::vector<std::string>>("path2");
       pSet.leg1 = parameterSet.getParameter<int>("leg1");
       pSet.leg2 = parameterSet.getParameter<int>("leg2");
       parameters_Tag_.push_back(pSet);
   }
   Initialize();
   return;

}


NtuplizerQCD::~NtuplizerQCD()
{}


//
// member functions
//

// ------------ method called for each event  ------------
void NtuplizerQCD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    Initialize();

    indexevents_ = iEvent.id().event();
    runNumber_ = iEvent.id().run();
    lumi_ = iEvent.luminosityBlock();

    edm::Handle<pat::TauRefVector> taus;
    edm::Handle<pat::JetRefVector> jets;
    edm::Handle<edm::View<pat::GenericParticle> > genParticles;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<std::vector<reco::Vertex> >  vertices;
    edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
    edm::Handle<pat::METCollection> metHandle;


    iEvent.getByToken(tauTag_, taus);
    iEvent.getByToken(jetTag_, jets);
    if(useHLTMatch_)
      iEvent.getByToken(triggerObjects_, triggerObjects);
    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(VtxTag_, vertices);
    iEvent.getByToken(puTag_, puInfo);
    iEvent.getByToken(metTag_, metHandle);

    if(isMC_)
      iEvent.getByToken(genParticlesTag_, genParticles);

    std::vector<unsigned int> probesWritten;
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    for (unsigned int i = 0; i < jets->size(); i+=1)
    {
        const pat::JetRef jetTag = jets->at(i);
        bool tagId = passesJetID(jetTag, year_);
        // int tagId = jetTag->userInt("tightJetIdLepVeto");
        if (tagId)
        {
            for (unsigned int j = 0; j < jets->size(); j+=1)
            {
                if (j != i)
                {
                    bool isTagHLTmatched = false;
                    isProbeHLTmatched_ = false;
                    const pat::JetRef jetProbe = jets->at(j);
                    tauTriggerBitSet_.reset();
                    jetTagTriggerBitSet_.reset();
                    if (useHLTMatch_)
                    {
                        for (pat::TriggerObjectStandAlone obj: *triggerObjects)
                        {
                            const float dR_tag = deltaR(*jetTag, obj);
                            if (dR_tag < 0.3)
                            {
                                const edm::TriggerNames::Strings& triggerNames = names.triggerNames();
                                // Looking for the path index.
                                unsigned int x = 0;
                                bool foundTrigger = false;
                                for (const tParameterSet& parameter : parameters_Tag_)
                                {
                                    if ((parameter.hltPathIndex >= 0) && obj.hasPathName(triggerNames[parameter.hltPathIndex], true, false))
                                    {
                                        foundTrigger = true;
                                        const std::vector<std::string>& filters = parameter.hltFilters1;
                                        if (hasFilters(obj, filters))
                                        {
                                            jetTagTriggerBitSet_[x] = true;
                                        }
                                    }
                                    x += 1;
                                }
                                if (foundTrigger)
                                {
                                    isTagHLTmatched = true;
                                }
                            }
                            
                            const float dR_probe = deltaR(*jetProbe, obj);
                            if (dR_probe < 0.3)
                            {
                                const edm::TriggerNames::Strings& triggerNames = names.triggerNames();
                                // Looking for the path index.
                                unsigned int x = 0;
                                bool foundTrigger = false;
                                for (const tParameterSet& parameter : parameters_)
                                {
                                    if ((parameter.hltPathIndex >= 0) && obj.hasPathName(triggerNames[parameter.hltPathIndex], true, false))
                                    {
                                        // isTagHLTmatched ? std::cout << "Tag is HLT matched." : std::cout << "Tag is NOT HLT matched.";
                                        // std::cout << std::endl;
                                        // std::cout << "Object is matched to probe jet: " << std::endl;
                                        // std::cout << "jet kinematics: " << "pT: " << jetProbe->pt() << " eta: " << jetProbe->eta() << " phi: " << jetProbe->phi() << std::endl;
                                        // std::cout << "Kinematics: " << "pT: " << obj.pt() << " eta: " << obj.eta() << " phi: " << obj.phi() << std::endl;
                                        foundTrigger = true;
                                        const std::vector<std::string>& filters = parameter.hltFilters1;
                                        if (hasFilters(obj, filters))
                                        {
                                            tauTriggerBitSet_[x] = true;
                                        }
                                    }
                                    x += 1;
                                }
                                if (foundTrigger)
                                {
                                    isProbeHLTmatched_ = true;
                                }
                            }

                            if (obj.hasPathName(filterPath_ + "*", false, false))
                            {
                                for (std::vector<std::string>::reverse_iterator filterName = triggerModules_.rbegin(); filterName != triggerModules_.rend(); filterName+=1)
                                {
                                    if (obj.hasFilterLabel(*filterName))
                                    {
                                        if (triggerModules_.rend() - filterName > lastFilter_)
                                        {
                                            lastFilter_ = triggerModules_.rend() - filterName;
                                            // std::cout << "FilterName: " << *filterName  << " is accepted? " << obj.hasPathLastFilterAccepted() << std::endl;
                                            // std::cout << "lastFilter updated to: " << lastFilter_ << std::endl;
                                            // std::cout << "Kinematics: " << "pT: " << obj.pt() << " eta: " << obj.eta() << " phi: " << obj.phi() << std::endl;
                                            // bool isBoth = obj.hasPathName( filterPath_ + "*", true, true );
                                            // bool isL3   = obj.hasPathName( filterPath_ + "*", false, true );
                                            // bool isLF   = obj.hasPathName( filterPath_ + "*", true, false );
                                            // bool isNone = obj.hasPathName( filterPath_ + "*", false, false );
                                            // std::cout << "   " << filterPath_ + "*";
                                            // if (isBoth) std::cout << "(L,3)";
                                            // if (isL3 && !isBoth) std::cout << "(*,3)";
                                            // if (isLF && !isBoth) std::cout << "(L,*)";
                                            // if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
                                            //std::cout << std::endl;
                                        }
                                    }
                                }
                            }
                        }
                    }   
                    else
                    {
                        isTagHLTmatched = true;
                    }

                    tauTriggerBits_ = tauTriggerBitSet_.to_ulong();
                    jetTagTriggerBits_ = jetTagTriggerBitSet_.to_ulong();

                    //std::cout << "Tag pT:" << jetTag->pt() << std::endl;
                    //std::cout << "Probe pT:" << jetProbe->pt() << std::endl;
                    bool foundTau = false;
                    for (pat::TauRef tau : *taus)
                    {
                        //std::cout << "Tau pT:" << tau->pt() << std::endl;
                        const float dR_tau = deltaR(*jetProbe, *tau);
                        if (dR_tau < 0.3)
                        {
                            foundTau = true;
                            tauPt_ = tau->pt();
                            tauEta_ = tau->eta();
                            tauPhi_ = tau->phi();
                            tauCharge_ = tau->charge();
                            tauDM_ = tau->decayMode();
                            tauTrkPt_ = tau->leadChargedHadrCand()->pt();

                            tauDecayModeFinding_ = tau->tauID("decayModeFinding");
                            tauDecayModeFindingNewDMs_ = tau->tauID("decayModeFindingNewDMs");

                            tauByLooseCombinedIsolationDeltaBetaCorr3Hits_ = tau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
                            tauByMediumCombinedIsolationDeltaBetaCorr3Hits_ = tau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
                            tauByTightCombinedIsolationDeltaBetaCorr3Hits_ = tau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");

                            tauByIsolationMVArun2017v2DBoldDMwLTraw2017_ = tau->tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017");
                            tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_ = tau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017");
                            tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_ = tau->tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017");
                            tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_ = tau->tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017");
                            tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_ = tau->tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017");
                            tauByTightIsolationMVArun2017v2DBoldDMwLT2017_ = tau->tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017");
                            tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_ = tau->tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017");
                            tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_ = tau->tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017");

                            tauByIsolationMVArun2017v2DBnewDMwLTraw2017_ = tau->tauID("byIsolationMVArun2017v2DBnewDMwLTraw2017");
                            tauByVVLooseIsolationMVArun2017v2DBnewDMwLT2017_ = tau->tauID("byVVLooseIsolationMVArun2017v2DBnewDMwLT2017");
                            tauByVLooseIsolationMVArun2017v2DBnewDMwLT2017_ = tau->tauID("byVLooseIsolationMVArun2017v2DBnewDMwLT2017");
                            tauByLooseIsolationMVArun2017v2DBnewDMwLT2017_ = tau->tauID("byLooseIsolationMVArun2017v2DBnewDMwLT2017");
                            tauByMediumIsolationMVArun2017v2DBnewDMwLT2017_ = tau->tauID("byMediumIsolationMVArun2017v2DBnewDMwLT2017");
                            tauByTightIsolationMVArun2017v2DBnewDMwLT2017_ = tau->tauID("byTightIsolationMVArun2017v2DBnewDMwLT2017");
                            tauByVTightIsolationMVArun2017v2DBnewDMwLT2017_ = tau->tauID("byVTightIsolationMVArun2017v2DBnewDMwLT2017");
                            tauByVVTightIsolationMVArun2017v2DBnewDMwLT2017_ = tau->tauID("byVVTightIsolationMVArun2017v2DBnewDMwLT2017");

                            tauByDeepTau2017v2p1VSjetraw_ = tau->tauID("byDeepTau2017v2p1VSjetraw");
                            tauByVVVLooseDeepTau2017v2p1VSjet_ = tau->tauID("byVVVLooseDeepTau2017v2p1VSjet");
                            tauByVVLooseDeepTau2017v2p1VSjet_ = tau->tauID("byVVLooseDeepTau2017v2p1VSjet");
                            tauByVLooseDeepTau2017v2p1VSjet_ = tau->tauID("byVLooseDeepTau2017v2p1VSjet");
                            tauByLooseDeepTau2017v2p1VSjet_ = tau->tauID("byLooseDeepTau2017v2p1VSjet");
                            tauByMediumDeepTau2017v2p1VSjet_ = tau->tauID("byMediumDeepTau2017v2p1VSjet");
                            tauByTightDeepTau2017v2p1VSjet_ = tau->tauID("byTightDeepTau2017v2p1VSjet");
                            tauByVTightDeepTau2017v2p1VSjet_ = tau->tauID("byVTightDeepTau2017v2p1VSjet");
                            tauByVVTightDeepTau2017v2p1VSjet_ = tau->tauID("byVVTightDeepTau2017v2p1VSjet");

                            tauByDeepTau2017v2p1VSeraw_ = tau->tauID("byDeepTau2017v2p1VSeraw");
                            tauByVVVLooseDeepTau2017v2p1VSe_ = tau->tauID("byVVVLooseDeepTau2017v2p1VSe");
                            tauByVVLooseDeepTau2017v2p1VSe_ = tau->tauID("byVVLooseDeepTau2017v2p1VSe");
                            tauByVLooseDeepTau2017v2p1VSe_ = tau->tauID("byVLooseDeepTau2017v2p1VSe");
                            tauByLooseDeepTau2017v2p1VSe_ = tau->tauID("byLooseDeepTau2017v2p1VSe");
                            tauByMediumDeepTau2017v2p1VSe_ = tau->tauID("byMediumDeepTau2017v2p1VSe");
                            tauByTightDeepTau2017v2p1VSe_ = tau->tauID("byTightDeepTau2017v2p1VSe");
                            tauByVTightDeepTau2017v2p1VSe_ = tau->tauID("byVTightDeepTau2017v2p1VSe");
                            tauByVVTightDeepTau2017v2p1VSe_ = tau->tauID("byVVTightDeepTau2017v2p1VSe");

                            tauByDeepTau2017v2p1VSmuraw_ = tau->tauID("byDeepTau2017v2p1VSmuraw");
                            tauByVLooseDeepTau2017v2p1VSmu_ = tau->tauID("byVLooseDeepTau2017v2p1VSmu");
                            tauByLooseDeepTau2017v2p1VSmu_ = tau->tauID("byLooseDeepTau2017v2p1VSmu");
                            tauByMediumDeepTau2017v2p1VSmu_ = tau->tauID("byMediumDeepTau2017v2p1VSmu");
                            tauByTightDeepTau2017v2p1VSmu_ = tau->tauID("byTightDeepTau2017v2p1VSmu");
                            
                            tauAgainstMuonLoose3_ = tau->tauID("againstMuonLoose3");
                            tauAgainstMuonTight3_ = tau->tauID("againstMuonTight3");
                            tauAgainstElectronVLooseMVA6_ = tau->tauID("againstElectronVLooseMVA6");
                            tauAgainstElectronLooseMVA6_ = tau->tauID("againstElectronLooseMVA6");
                            tauAgainstElectronMediumMVA6_ = tau->tauID("againstElectronMediumMVA6");
                            tauAgainstElectronTightMVA6_ = tau->tauID("againstElectronTightMVA6");
                            tauAgainstElectronVTightMVA6_ = tau->tauID("againstElectronVTightMVA6");
                            if (isMC_)
                            {
                                const edm::View<pat::GenericParticle>* genparts = genParticles.product();
                                tau_genindex_ = GenIndex(tau, genparts);
                            }
                            break;
                        }
                    }
                    
                    jetTagPt_ = jetTag->pt();
                    jetTagEta_ = jetTag->eta();
                    jetTagPhi_ = jetTag->phi();

                    jetProbePt_ = jetProbe->pt();
                    jetProbeEta_ = jetProbe->eta();
                    jetProbePhi_ = jetProbe->phi();

                    jetProbePartonFlavour_ = jetProbe->partonFlavour();
                    jetProbeNeutHadFrac_ = jetProbe->neutralHadronEnergyFraction();
                    jetProbeNeutEMFrac_ = jetProbe->neutralEmEnergyFraction();
                    jetProbeNumConsts_ = jetProbe->nConstituents();
                    jetProbeMuonFrac_ = jetProbe->muonEnergyFraction();
                    jetProbeChargedHadFrac_ = jetProbe->chargedHadronEnergyFraction();
                    jetProbeChargedMult_ = jetProbe->chargedMultiplicity();
                    jetProbeChargedEMFrac_ = jetProbe->chargedEmEnergyFraction();
                    jetProbeN60_ = jetProbe->n60();
                    jetProbeN90_ = jetProbe->n90();
                    jetProbeQgLikelihood_ = jetProbe->userFloat("QGTagger:qgLikelihood");

                    nJets_ = jets->size();

                    metPt_ = (*metHandle)[0].pt();
                    metPhi_ = (*metHandle)[0].phi();

                    Nvtx_ = vertices->size();

                    nTruePU_ = -99;
                    if (isMC_)
                    {
                        std::vector<PileupSummaryInfo>::const_iterator PVI;
                        for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI)
                        {
                            if(PVI->getBunchCrossing() == 0)
                            {
                                float nTrueInt = PVI->getTrueNumInteractions();
                                nTruePU_ = nTrueInt;
                                break;
                            }
                        }
                    }
                    if (isTagHLTmatched && foundTau && (std::find(probesWritten.begin(), probesWritten.end(), j) == probesWritten.end()))
                    {
                        probesWritten.push_back(j);
                        //std::cout << "Fill event with: EventNumber " << indexevents_ << " RunNumber " << runNumber_ << " and LumiSection: " << lumi_ << std::endl;
                        tree_->Fill();
                    }
                    // else
                    // {
                    //     std::cout << "Did not fill event with: EventNumber " << indexevents_ << " RunNumber " << runNumber_ << " and LumiSection: " << lumi_ << std::endl;
                    //     std::cout << "reason: " << isTagHLTmatched ? std::cout << " no tau found" : std::cout << " tag not HLT matched";
                    //     std::cout << std::endl;
                    // }
                }
            }
        }
    }
}


// ------------ method called once each job just before starting event loop  ------------
void NtuplizerQCD::beginJob()
{
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>(treeName_.c_str(), treeName_.c_str());

    tree_->Branch("EventNumber", &indexevents_, "EventNumber/l");
    tree_->Branch("RunNumber", &runNumber_, "RunNumber/I");
    tree_->Branch("lumi", &lumi_, "lumi/I");

    tree_->Branch("tauPt", &tauPt_, "tauPt/F");
    tree_->Branch("tauEta", &tauEta_, "tauEta/F");
    tree_->Branch("tauPhi", &tauPhi_, "tauPhi/F");
    tree_->Branch("tauCharge", &tauCharge_, "tauCharge/I");
    tree_->Branch("tauDM", &tauDM_, "tauDM/I");
    tree_->Branch("tauTrkPt", &tauTrkPt_, "tauTrkPt/F");
    tree_->Branch("tau_genindex", &tau_genindex_, "tau_genindex/I");

    tree_->Branch("decayModeFinding", &tauDecayModeFinding_, "decayModeFinding/O");
    tree_->Branch("decayModeFindingNewDMs", &tauDecayModeFindingNewDMs_, "decayModeFindingsNewDMs/O");

    tree_->Branch("byLooseCombinedIsolationDeltaBetaCorr3Hits", &tauByLooseCombinedIsolationDeltaBetaCorr3Hits_, "byLooseCombinedIsolationDeltaBetaCorr3Hits/O");
    tree_->Branch("byMediumCombinedIsolationDeltaBetaCorr3Hits", &tauByMediumCombinedIsolationDeltaBetaCorr3Hits_, "byMediumCombinedIsolationDeltaBetaCorr3Hits/O");
    tree_->Branch("byTightCombinedIsolationDeltaBetaCorr3Hits", &tauByTightCombinedIsolationDeltaBetaCorr3Hits_, "byTightCombinedIsolationDeltaBetaCorr3Hits/O");
    
    tree_->Branch("byIsolationMVArun2017v2DBoldDMwLTraw2017", &tauByIsolationMVArun2017v2DBoldDMwLTraw2017_, "byIsolationMVArun2017v2DBoldDMwLTraw2017/F");
    tree_->Branch("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017", &tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_, "byVVLooseIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("byVLooseIsolationMVArun2017v2DBoldDMwLT2017", &tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_, "byVLooseIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("byLooseIsolationMVArun2017v2DBoldDMwLT2017", &tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_, "byLooseIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("byMediumIsolationMVArun2017v2DBoldDMwLT2017", &tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_, "byMediumIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("byTightIsolationMVArun2017v2DBoldDMwLT2017", &tauByTightIsolationMVArun2017v2DBoldDMwLT2017_, "byTightIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("byVTightIsolationMVArun2017v2DBoldDMwLT2017", &tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_, "byVTightIsolationMVArun2017v2DBoldDMwLT2017/O");
    tree_->Branch("byVVTightIsolationMVArun2017v2DBoldDMwLT2017", &tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_, "byVVTightIsolationMVArun2017v2DBoldDMwLT2017/O");

    tree_->Branch("byIsolationMVArun2017v2DBnewDMwLTraw2017", &tauByIsolationMVArun2017v2DBnewDMwLTraw2017_, "byIsolationMVArun2017v2DBnewDMwLTraw2017/F");
    tree_->Branch("byVVLooseIsolationMVArun2017v2DBnewDMwLT2017", &tauByVVLooseIsolationMVArun2017v2DBnewDMwLT2017_, "byVVLooseIsolationMVArun2017v2DBnewDMwLT2017/O");
    tree_->Branch("byVLooseIsolationMVArun2017v2DBnewDMwLT2017", &tauByVLooseIsolationMVArun2017v2DBnewDMwLT2017_, "byVLooseIsolationMVArun2017v2DBnewDMwLT2017/O");
    tree_->Branch("byLooseIsolationMVArun2017v2DBnewDMwLT2017", &tauByLooseIsolationMVArun2017v2DBnewDMwLT2017_, "byLooseIsolationMVArun2017v2DBnewDMwLT2017/O");
    tree_->Branch("byMediumIsolationMVArun2017v2DBnewDMwLT2017", &tauByMediumIsolationMVArun2017v2DBnewDMwLT2017_, "byMediumIsolationMVArun2017v2DBnewDMwLT2017/O");
    tree_->Branch("byTightIsolationMVArun2017v2DBnewDMwLT2017", &tauByTightIsolationMVArun2017v2DBnewDMwLT2017_, "byTightIsolationMVArun2017v2DBnewDMwLT2017/O");
    tree_->Branch("byVTightIsolationMVArun2017v2DBnewDMwLT2017", &tauByVTightIsolationMVArun2017v2DBnewDMwLT2017_, "byVTightIsolationMVArun2017v2DBnewDMwLT2017/O");
    tree_->Branch("byVVTightIsolationMVArun2017v2DBnewDMwLT2017", &tauByVVTightIsolationMVArun2017v2DBnewDMwLT2017_, "byVVTightIsolationMVArun2017v2DBnewDMwLT2017/O");

    tree_->Branch("byDeepTau2017v2p1VSjetraw", &tauByDeepTau2017v2p1VSjetraw_, "byDeepTau2017v2p1VSjetraw/F");
    tree_->Branch("byVVVLooseDeepTau2017v2p1VSjet", &tauByVVVLooseDeepTau2017v2p1VSjet_, "byVVVLooseDeepTau2017v2p1VSjet/O");
    tree_->Branch("byVVLooseDeepTau2017v2p1VSjet", &tauByVVLooseDeepTau2017v2p1VSjet_, "byVVLooseDeepTau2017v2p1VSjet/O");
    tree_->Branch("byVLooseDeepTau2017v2p1VSjet", &tauByVLooseDeepTau2017v2p1VSjet_, "byVLooseDeepTau2017v2p1VSjet/O");
    tree_->Branch("byLooseDeepTau2017v2p1VSjet", &tauByLooseDeepTau2017v2p1VSjet_, "byLooseDeepTau2017v2p1VSjet/O");
    tree_->Branch("byMediumDeepTau2017v2p1VSjet", &tauByMediumDeepTau2017v2p1VSjet_, "byMediumDeepTau2017v2p1VSjet/O");
    tree_->Branch("byTightDeepTau2017v2p1VSjet", &tauByTightDeepTau2017v2p1VSjet_, "byTightDeepTau2017v2p1VSjet/O");
    tree_->Branch("byVTightDeepTau2017v2p1VSjet", &tauByVTightDeepTau2017v2p1VSjet_, "byVTightDeepTau2017v2p1VSjet/O");
    tree_->Branch("byVVTightDeepTau2017v2p1VSjet", &tauByVVTightDeepTau2017v2p1VSjet_, "byVVTightDeepTau2017v2p1VSjet/O");

    tree_->Branch("byDeepTau2017v2p1VSeraw", &tauByDeepTau2017v2p1VSeraw_, "byDeepTau2017v2p1VSeraw/F");
    tree_->Branch("byVVVLooseDeepTau2017v2p1VSe", &tauByVVVLooseDeepTau2017v2p1VSe_, "byVVVLooseDeepTau2017v2p1VSe/O");
    tree_->Branch("byVVLooseDeepTau2017v2p1VSe", &tauByVVLooseDeepTau2017v2p1VSe_, "byVVLooseDeepTau2017v2p1VSe/O");
    tree_->Branch("byVLooseDeepTau2017v2p1VSe", &tauByVLooseDeepTau2017v2p1VSe_, "byVLooseDeepTau2017v2p1VSe/O");
    tree_->Branch("byLooseDeepTau2017v2p1VSe", &tauByLooseDeepTau2017v2p1VSe_, "byLooseDeepTau2017v2p1VSe/O");
    tree_->Branch("byMediumDeepTau2017v2p1VSe", &tauByMediumDeepTau2017v2p1VSe_, "byMediumDeepTau2017v2p1VSe/O");
    tree_->Branch("byTightDeepTau2017v2p1VSe", &tauByTightDeepTau2017v2p1VSe_, "byTightDeepTau2017v2p1VSe/O");
    tree_->Branch("byVTightDeepTau2017v2p1VSe", &tauByVTightDeepTau2017v2p1VSe_, "byVTightDeepTau2017v2p1VSe/O");
    tree_->Branch("byVVTightDeepTau2017v2p1VSe", &tauByVVTightDeepTau2017v2p1VSe_, "byVVTightDeepTau2017v2p1VSe/O");
    
    tree_->Branch("byDeepTau2017v2p1VSmuraw", &tauByDeepTau2017v2p1VSmuraw_, "byDeepTau2017v2p1VSmuraw/F");
    tree_->Branch("byVLooseDeepTau2017v2p1VSmu", &tauByVLooseDeepTau2017v2p1VSmu_, "byVLooseDeepTau2017v2p1VSmu/O");
    tree_->Branch("byLooseDeepTau2017v2p1VSmu", &tauByLooseDeepTau2017v2p1VSmu_, "byLooseDeepTau2017v2p1VSmu/O");
    tree_->Branch("byMediumDeepTau2017v2p1VSmu", &tauByMediumDeepTau2017v2p1VSmu_, "byMediumDeepTau2017v2p1VSmu/O");
    tree_->Branch("byTightDeepTau2017v2p1VSmu", &tauByTightDeepTau2017v2p1VSmu_, "byTightDeepTau2017v2p1VSmu/O");
          
    tree_->Branch("againstMuonLoose3", &tauAgainstMuonLoose3_,"againstMuonLoose3/O");
    tree_->Branch("againstMuonTight3", &tauAgainstMuonTight3_, "againstMuonTight3/O");
    tree_->Branch("againstElectronVLooseMVA6", &tauAgainstElectronVLooseMVA6_, "againstElectronVLooseMVA6/O");
    tree_->Branch("againstElectronLooseMVA6", &tauAgainstElectronLooseMVA6_, "againstElectronLooseMVA6/O");
    tree_->Branch("againstElectronMediumMVA6", &tauAgainstElectronMediumMVA6_, "againstElectronMediumMVA6/O");
    tree_->Branch("againstElectronTightMVA6", &tauAgainstElectronTightMVA6_, "againstElectronTightMVA6/O");
    tree_->Branch("againstElectronVTightMVA6", &tauAgainstElectronVTightMVA6_, "againstElectronVTightMVA6/O");
    tree_->Branch("tauTriggerBits", &tauTriggerBits_, "tauTriggerBits/l");
    tree_->Branch("jetTagTriggerBits", &jetTagTriggerBits_, "jetTagTriggerBits/l");

    tree_->Branch("jetTagPt", &jetTagPt_, "jetTagPt/F");
    tree_->Branch("jetTagEta", &jetTagEta_, "jetTagEta/F");
    tree_->Branch("jetTagPhi", &jetTagPhi_, "jetTagPhi/F");

    tree_->Branch("jetProbePt", &jetProbePt_, "jetProbePt/F");
    tree_->Branch("jetProbeEta", &jetProbeEta_, "jetProbeEta/F");
    tree_->Branch("jetProbePhi", &jetProbePhi_, "jetProbePhi/F");
    tree_->Branch("jetProbePartonFlavour", &jetProbePartonFlavour_, "jetProbePartonFlavour/I");
    tree_->Branch("jetProbeNeutHadFrac", &jetProbeNeutHadFrac_, "jetProbeNeutHadFrac/F");
    tree_->Branch("jetProbeNeutEMFrac", &jetProbeNeutEMFrac_, "jetProbeNeutEMFrac/F");
    tree_->Branch("jetProbeNumConsts", &jetProbeNumConsts_, "jetProbeNumConsts/I");
    tree_->Branch("jetProbeMuonFrac", &jetProbeMuonFrac_, "jetProbeMuonFrac/F");
    tree_->Branch("jetProbeChargedHadFrac", &jetProbeChargedHadFrac_, "jetProbeChargedHadFrac/F");
    tree_->Branch("jetProbeChargedMult", &jetProbeChargedMult_, "jetProbeChargedMult/I");
    tree_->Branch("jetProbeChargedEMFrac", &jetProbeChargedEMFrac_, "jetProbeChargedEMFrac/F");
    tree_->Branch("jetProbeN60", &jetProbeN60_, "jetProbeN60/I");
    tree_->Branch("jetProbeN90", &jetProbeN90_, "jetProbeN90/I");
    tree_->Branch("jetProbeQgLikelihood", &jetProbeQgLikelihood_, "jetProbeQgLikelihood/F");
    tree_->Branch("nJets", &nJets_, "nJets/I");

    tree_->Branch("metPt", &metPt_, "metPt/F");
    tree_->Branch("metPhi", &metPhi_, "metPhi/F");

    tree_->Branch("isMatched", &isProbeHLTmatched_, "isMatched/O");

    tree_->Branch("Nvtx", &Nvtx_, "Nvtx/I");
    tree_->Branch("nTruePU", &nTruePU_, "nTruePU/F");

    tree_->Branch("lastFilter", &lastFilter_, "lastFilter/I");

    return;
}

// ------------ method called once each job just after ending the event loop  ------------
void NtuplizerQCD::endJob() 
{
    return;
}


void NtuplizerQCD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    Bool_t changedConfig = false;

    if (!hltConfig_.init(iRun, iSetup, processName_.process(), changedConfig))
    {
        edm::LogError("HLTMatchingFilter") << "Initialization of HLTConfigProvider failed!!";
        return;
    }

    const edm::TriggerNames::Strings& triggerNames = hltConfig_.triggerNames();
    std::cout << " ===== LOOKING FOR THE PATH INDEXES =====" << std::endl;
    for (tParameterSet& parameter : parameters_Tag_){
        const std::string& hltPath = parameter.hltPath;
        bool found = false; 
        for(unsigned int j=0; j < triggerNames.size(); j++)
        {
            //std::cout << triggerNames[j] << std::endl;
            if (triggerNames[j].find(hltPath) != std::string::npos) {
                found = true;
                parameter.hltPathIndex = j;
                                           
               std::cout << "### FOUND AT INDEX #" << j << " --> " << triggerNames[j] << std::endl;
            }
        }
        if (!found) parameter.hltPathIndex = -1;
    }

    for (tParameterSet& parameter : parameters_){
        const std::string& hltPath = parameter.hltPath;
        bool found = false;
        for(unsigned int j=0; j < triggerNames.size(); j++)
        {
            if (triggerNames[j].find(hltPath) != std::string::npos) {
                found = true;
                parameter.hltPathIndex = j;

                std::cout << "### FOUND AT INDEX #" << j << " --> " << triggerNames[j] << std::endl;
                // Look for the trigger filters running in this configuration.
                if (hltPath==filterPath_)
                {
                    lastFilterInd_ = j;
                }
            }
        }
        if (!found) parameter.hltPathIndex = -1;
    }

    // Get trigger modules which ran with saveTags option, e.g. important EDFilters
    triggerModules_= hltConfig_.saveTagsModules(lastFilterInd_);
    for (const std::string triggerModule: triggerModules_)
    {
        filterLabel_ = triggerModule;
        filterLabelsTree_->Fill();
    }

    return;
}

void NtuplizerQCD::Initialize()
{
    indexevents_ = 0;
    runNumber_ = 0;
    lumi_ = 0;

    tauPt_ = -1.;
    tauEta_ = -999;
    tauPhi_ = -999;
    tauCharge_ = 0;
    tauDM_ = -1;
    tauTrkPt_ = -1.;
    tau_genindex_ = -1;
    tauDecayModeFinding_ = 0;
    tauDecayModeFindingNewDMs_ = 0;

    tauByLooseCombinedIsolationDeltaBetaCorr3Hits_ = 0;
    tauByMediumCombinedIsolationDeltaBetaCorr3Hits_ = 0;
    tauByTightCombinedIsolationDeltaBetaCorr3Hits_ = 0;

    tauByIsolationMVArun2017v2DBoldDMwLTraw2017_ = -1.;
    tauByVVLooseIsolationMVArun2017v2DBoldDMwLT2017_ = 0;
    tauByVLooseIsolationMVArun2017v2DBoldDMwLT2017_ = 0;
    tauByLooseIsolationMVArun2017v2DBoldDMwLT2017_ = 0;
    tauByMediumIsolationMVArun2017v2DBoldDMwLT2017_ = 0;
    tauByTightIsolationMVArun2017v2DBoldDMwLT2017_ = 0;
    tauByVTightIsolationMVArun2017v2DBoldDMwLT2017_ = 0;
    tauByVVTightIsolationMVArun2017v2DBoldDMwLT2017_ = 0;

    tauByIsolationMVArun2017v2DBnewDMwLTraw2017_ = -1.;
    tauByVVLooseIsolationMVArun2017v2DBnewDMwLT2017_ = 0;
    tauByVLooseIsolationMVArun2017v2DBnewDMwLT2017_ = 0;
    tauByLooseIsolationMVArun2017v2DBnewDMwLT2017_ = 0;
    tauByMediumIsolationMVArun2017v2DBnewDMwLT2017_ = 0;
    tauByTightIsolationMVArun2017v2DBnewDMwLT2017_ = 0;
    tauByVTightIsolationMVArun2017v2DBnewDMwLT2017_ = 0;
    tauByVVTightIsolationMVArun2017v2DBnewDMwLT2017_ = 0;

    tauByDeepTau2017v2p1VSjetraw_ = -1.;
    tauByVVVLooseDeepTau2017v2p1VSjet_ = 0;
    tauByVVLooseDeepTau2017v2p1VSjet_ = 0;
    tauByVLooseDeepTau2017v2p1VSjet_ = 0;
    tauByLooseDeepTau2017v2p1VSjet_ = 0;
    tauByMediumDeepTau2017v2p1VSjet_ = 0;
    tauByTightDeepTau2017v2p1VSjet_ = 0;
    tauByVTightDeepTau2017v2p1VSjet_ = 0;
    tauByVVTightDeepTau2017v2p1VSjet_ = 0;

    tauByDeepTau2017v2p1VSeraw_ = -1.;
    tauByVVVLooseDeepTau2017v2p1VSe_ = 0;
    tauByVVLooseDeepTau2017v2p1VSe_ = 0;
    tauByVLooseDeepTau2017v2p1VSe_ = 0;
    tauByLooseDeepTau2017v2p1VSe_ = 0;
    tauByMediumDeepTau2017v2p1VSe_ = 0;
    tauByTightDeepTau2017v2p1VSe_ = 0;
    tauByVTightDeepTau2017v2p1VSe_ = 0;
    tauByVVTightDeepTau2017v2p1VSe_ = 0;

    tauByDeepTau2017v2p1VSmuraw_ = -1.;
    tauByVLooseDeepTau2017v2p1VSmu_ = 0;
    tauByLooseDeepTau2017v2p1VSmu_ = 0;
    tauByMediumDeepTau2017v2p1VSmu_ = 0;
    tauByTightDeepTau2017v2p1VSmu_ = 0;

    tauAgainstMuonLoose3_ = 0;
    tauAgainstMuonTight3_ = 0;
    tauAgainstElectronVLooseMVA6_ = 0;
    tauAgainstElectronLooseMVA6_ = 0;
    tauAgainstElectronMediumMVA6_ = 0;
    tauAgainstElectronTightMVA6_ = 0;
    tauAgainstElectronVTightMVA6_ = 0;

    jetTagPt_ = -1.;
    jetTagEta_ = -999;
    jetTagPhi_ = -999;

    jetProbePt_ = -1.;
    jetProbeEta_ = -999;
    jetProbePhi_ = -999;
    jetProbePartonFlavour_ = -99;
    jetProbeNeutHadFrac_ = -1.;
    jetProbeNeutEMFrac_ = -1.;
    jetProbeNumConsts_ = -1;
    jetProbeMuonFrac_ = -1.;
    jetProbeChargedHadFrac_ = -1.;
    jetProbeChargedMult_ = -1;
    jetProbeChargedEMFrac_ = -1.;
    jetProbeN60_ = -1;
    jetProbeN90_ = -1;
    jetProbeQgLikelihood_ = -1.;

    nJets_ = -1;

    metPt_ = -1.;
    metPhi_ = -999.;

    isTagHLTmatched_ = false;
    isProbeHLTmatched_ = false;

    Nvtx_ = 0;
    nTruePU_ = 0.;

    lastFilter_ = 0;
    tauTriggerBits_ = 999;
    jetTagTriggerBits_ = 999;

    return; 
}

void NtuplizerQCD::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    return;
}

bool NtuplizerQCD::hasFilters(const pat::TriggerObjectStandAlone&  obj , const std::vector<std::string>& filtersToLookFor)
{
    const std::vector<std::string>& eventLabels = obj.filterLabels();
    for (const std::string& filter : filtersToLookFor)
    {
        //Looking for matching filters
        bool found = false;
        for (const std::string& label : eventLabels)
        {
            //if (label == std::string("hltOverlapFilterIsoMu17MediumIsoPFTau40Reg"))
            if (label == filter)
            {
                //std::cout << "#### FOUND FILTER " << label << " == " << filter << " ####" << std::endl;
                found = true;
            }
        }
        if(!found)
        {
            return false;
        }
    }
    return true;
}

int NtuplizerQCD::GenIndex(const pat::TauRef& tau, const edm::View<pat::GenericParticle>* genparts)
{
    float dRmin = 1.0;
    int genMatchInd = -1;

    for(edm::View<pat::GenericParticle>::const_iterator genpart = genparts->begin(); genpart!=genparts->end();++genpart)
    {

        int flags = genpart->userInt ("generalGenFlags");
        int apdg = abs(genpart->pdgId());
        float pT = genpart->p4().pt();

        if( !( apdg==11 || apdg==13 || apdg==66615) ) continue;

        if( apdg==11 || apdg==13)
        {
            if( !(pT>8 && (flags&1 || (flags>>5)&1)) ) continue;
        }
        else if(apdg==66615)
        {
            int tauMothInd = genpart->userInt("TauMothIndex");
            pat::GenericParticle mother = (*genparts)[tauMothInd];
            int flags_tau = mother.userInt ("generalGenFlags");
            if( !(pT>15 && flags_tau&1) ) continue;
        }

        float dR = deltaR(*tau,*genpart);
        if(dR<0.2 && dR<dRmin)
        {
            dRmin = dR;
            if(apdg==11)
            {
                if(flags&1) genMatchInd = 1;
                else if((flags>>5)&1) genMatchInd = 3;
            }
            else if(apdg==13)
            {
                if(flags&1) genMatchInd = 2;
                else if((flags>>5)&1) genMatchInd = 4;
            }
            else if(apdg==66615)
                genMatchInd = 5;
        }

    }

    return genMatchInd;
}

bool NtuplizerQCD::passesJetID(const pat::JetRef& jet, const int year)
{
    bool passesID = false;
    double NHF = jet->neutralHadronEnergyFraction();
    double NEMF = jet->neutralEmEnergyFraction();
    double CHF  = jet->chargedHadronEnergyFraction();
    double MUF  = jet->muonEnergyFraction();
    double CEMF = jet->chargedEmEnergyFraction();
    int NumConst = jet->chargedMultiplicity()+jet->neutralMultiplicity();
    int NumNeutralParticles = jet->neutralMultiplicity();
    int CHM      = jet->chargedMultiplicity(); 
    if (year == 2016)
    {
        // Use tightLepVetoJetID
        if (abs(jet->eta()) <= 2.7)
        {
            passesID = passesID || (NHF < 0.90 && NEMF < 0.90 && NumConst > 1 && MUF < 0.8 && ((abs(jet->eta()) <= 2.4 && CHF > 0 && CHM > 0 && CEMF < 0.90) || abs(jet->eta()) > 2.4 ) && abs(jet->eta()) <= 2.7);
        }
        // use tightID where tightLepVetJetID not applicable.
        else if ((abs(jet->eta()) > 2.7) && (abs(jet->eta()) <= 3.0))
        {
            passesID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 && abs(jet->eta())>2.7 && abs(jet->eta())<=3.0);
        }
        else if (abs(jet->eta()) > 3.0)
        {
            passesID = (NEMF<0.90 && NumNeutralParticles>10 && abs(jet->eta())>3.0 );
        }
    }
    else if (year == 2017)
    {
        // Use tightLepVetoJetID
        if (abs(jet->eta()) <= 2.7)
        {
            passesID = passesID || (NHF < 0.90 && NEMF < 0.90 && NumConst > 1 && MUF < 0.8 && ((abs(jet->eta()) <= 2.4 && CHF > 0 && CHM > 0 && CEMF < 0.80) || abs(jet->eta()) > 2.4 ) && abs(jet->eta()) <= 2.7);
        }
        // use tightID where tightLepVetJetID not applicable.
        else if ((abs(jet->eta()) > 2.7) && (abs(jet->eta()) <= 3.0))
        {
            passesID = passesID || (NEMF>0.02 && NEMF < 0.99 && NumNeutralParticles>2 && abs(jet->eta())>2.7 && abs(jet->eta())<=3.0);
        }
        else if (abs(jet->eta()) > 3.0)
        {
            passesID = passesID || (NEMF<0.90 && NumNeutralParticles>10 && NHF > 0.02 && abs(jet->eta())>3.0 );
        }
    }
    else if (year == 2018)
    {
        // Use tightLepVetoJetID
        if (abs(jet->eta()) <= 2.6)
        {
            passesID = passesID || (abs(jet->eta())<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 );
        }
        else if ((abs(jet->eta()) > 2.6) && (abs(jet->eta()) <= 2.7))
        {
            passesID = passesID || ( abs(jet->eta())>2.6 && abs(jet->eta())<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 );
        }
        // use tightID where tightLepVetJetID not applicable.
        else if ((abs(jet->eta()) > 2.7) && (abs(jet->eta()) <= 3.0))
        {
            passesID = passesID || ( NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2 && abs(jet->eta())>2.7 && abs(jet->eta())<=3.0 );
        }
        else if (abs(jet->eta()) > 3.0)
        {
            passesID = passesID || (NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 && abs(jet->eta())>3.0 );
        }
    }
    else
    {
        std::cout << "No JetID defined for this year." << std::endl;
    }
    return passesID;
}

//define this as a plug-in
#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(NtuplizerQCD);
#endif
