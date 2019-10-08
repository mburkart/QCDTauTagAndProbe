import FWCore.ParameterSet.Config as cms
# filter HLT paths for T&P
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
# flake8: noqa

print "Running on MC"


HLTLIST_TAG = cms.VPSet(
    # TODO: Find some QCD Tag Trigger.
    cms.PSet(
        HLT=cms.string("HLT_PFJet40_v"),
        path1=cms.vstring("hltSinglePFJet40"),
        path2=cms.vstring(""),
        leg1=cms.int32(999),
        leg2=cms.int32(999)
    ),
    cms.PSet(
        HLT=cms.string("HLT_PFJet60_v"),
        path1=cms.vstring("hltSinglePFJet60"),
        path2=cms.vstring(""),
        leg1=cms.int32(999),
        leg2=cms.int32(999)
    ),
    cms.PSet(
        HLT=cms.string("HLT_PFJet80_v"),
        path1=cms.vstring("hltSinglePFJet80"),
        path2=cms.vstring(""),
        leg1=cms.int32(999),
        leg2=cms.int32(999)
    ),
    cms.PSet(
        HLT=cms.string("HLT_PFJet140_v"),
        path1=cms.vstring("hltSinglePFJet140"),
        path2=cms.vstring(""),
        leg1=cms.int32(999),
        leg2=cms.int32(999)
    ),
    cms.PSet(
        HLT=cms.string("HLT_PFJet200_v"),
        path1=cms.vstring("hltSinglePFJet200"),
        path2=cms.vstring(""),
        leg1=cms.int32(999),
        leg2=cms.int32(999)
    ),
    cms.PSet(
        HLT=cms.string("HLT_PFJet260_v"),
        path1=cms.vstring("hltSinglePFJet260"),
        path2=cms.vstring(""),
        leg1=cms.int32(999),
        leg2=cms.int32(999)
    ),
    cms.PSet(
        HLT=cms.string("HLT_PFJet320_v"),
        path1=cms.vstring("hltSinglePFJet320"),
        path2=cms.vstring(""),
        leg1=cms.int32(999),
        leg2=cms.int32(999)
    ),
    cms.PSet(
        HLT=cms.string("HLT_PFJet400_v"),
        path1=cms.vstring("hltSinglePFJet400"),
        path2=cms.vstring(""),
        leg1=cms.int32(999),
        leg2=cms.int32(999)
    ),
    cms.PSet(
        HLT=cms.string("HLT_PFJet450_v"),
        path1=cms.vstring("hltSinglePFJet450"),
        path2=cms.vstring(""),
        leg1=cms.int32(999),
        leg2=cms.int32(999)
    ),
    cms.PSet(
        HLT=cms.string("HLT_PFJet500_v"),
        path1=cms.vstring("hltSinglePFJet500"),
        path2=cms.vstring(""),
        leg1=cms.int32(999),
        leg2=cms.int32(999)
    ),
    cms.PSet(
        HLT=cms.string("HLT_PFJet550_v"),
        path1=cms.vstring("hltSinglePFJet550"),
        path2=cms.vstring(""),
        leg1=cms.int32(999),
        leg2=cms.int32(999)
    ),
)


HLTLIST = cms.VPSet(
    # SingleTau
    cms.PSet(
        HLT=cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1=cms.vstring("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched"),
        path2=cms.vstring(""),
        leg1=cms.int32(15),
        leg2=cms.int32(999)
    ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v"),
        path1 = cms.vstring ("hltSelectedPFTau180MediumChargedIsolationL1HLTMatched1Prong"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau200MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),
    cms.PSet (
        HLT = cms.string("HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_v"),
        path1 = cms.vstring ("hltSelectedPFTau220MediumChargedIsolationL1HLTMatched"),
        path2 = cms.vstring (""),
        leg1 = cms.int32(15),
        leg2 = cms.int32(999)
    ),
)

hltFilter = hlt.hltHighLevel.clone(
    TriggerResultsTag=cms.InputTag("TriggerResults", "", "HLT"),
    HLTPaths=["HLT_PFJet40_v*", "HLT_PFJet60_v*", "HLT_PFJet80_v*",
              "HLT_PFJet140_v*", "HLT_PFJet200_v*", "HLT_PFJet260_v*",
              "HLT_PFJet320_v*", "HLT_PFJet400_v*", "HLT_PFJet450_v*",
              "HLT_PFJet500_v*", "HLT_PFJet550_v*"],
    andOr=cms.bool(True),  # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw=cms.bool(True)  # if True: throws exception if a trigger path is invalid)
)

### ----------------------------------------------------------------------
### gen info, only from MC
### ----------------------------------------------------------------------
genInfo = cms.EDProducer("GenFiller",
        src=cms.InputTag("prunedGenParticles"),
        storeLightFlavAndGlu=cms.bool(True)  # if True, store also udcs and gluons (first copy)
)

goodJets = cms.EDFilter("PATJetRefSelector",
        src = cms.InputTag("slimmedJetsWithUserData"),
        cut = cms.string(
                'pt > 20 && abs(eta) < 2.7' # kinematics
                '&& muonMultiplicity == 0' # quality requirements
        ),
        filter = cms.bool(True)
)

tightJetIdLepVeto = cms.EDProducer("PatJetIDValueMapProducer",
        filterParams=cms.PSet(
                version = cms.string('WINTER17'),
                quality = cms.string('TIGHTLEPVETO'),
        ),
        src = cms.InputTag("slimmedJets")
)

slimmedJetsWithUserData = cms.EDProducer("PATJetUserDataEmbedder",
        src = cms.InputTag("slimmedJets"),
        userFloats = cms.PSet(
            qgLikelihood = cms.InputTag("QGTagger:qgLikelihood"),
        ),
)

goodTaus = cms.EDFilter("PATTauRefSelector",
        src = cms.InputTag("NewTauIDsEmbedded"),
        cut = cms.string(
                'pt > 20 && abs(eta) < 2.1 ' #kinematics
                '&& abs(charge) > 0 && abs(charge) < 2 ' #sometimes 2 prongs have charge != 1
                # '&& tauID("decayModeFinding") > 0.5 ' # tau ID
                '&& (tauID("decayModeFinding") > 0.5 || tauID("decayModeFindingNewDMs") > 0.5)'
                '&& (tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") > 0.5 || tauID("byVVVLooseDeepTau2017v2p1VSjet") > 0.5)' # tau iso - NOTE: can as well use boolean discriminators with WP
                '&& (tauID("againstMuonLoose3") > 0.5 || tauID("byVLooseDeepTau2017v2p1VSmu") > 0.5)' # anti Muon tight
                '&& (tauID("againstElectronVLooseMVA6") > 0.5 || tauID("byVVVLooseDeepTau2017v2p1VSe") > 0.5)' # anti-Ele loose
        ),
        filter = cms.bool(True)
)

patTriggerUnpacker = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
                                    patTriggerObjectsStandAlone = cms.InputTag("slimmedPatTrigger"),
                                    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                                    unpackFilterLabels = cms.bool(True)
)

muonsForVeto = cms.EDFilter("PATMuonRefSelector",
        src = cms.InputTag("slimmedMuons"),
        cut = cms.string(
            "pt > 10 && abs(eta) < 2.4 " # kinematics
            "&& ( (pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5 * pfIsolationR04().sumPUPt, 0.0)) / pt() ) < 0.3 " #isolation
            "&& isLooseMuon()" # quality requirement
        ),
)

bJetsForVeto = cms.EDFilter("PATJetRefSelector",
        src = cms.InputTag("selectedUpdatedPatJetsNewDFTraining"),
        cut = cms.string(
                'pt > 20 && abs(eta) < 2.4 ' #kinematics
                '&& (bDiscriminator("pfDeepFlavourJetTags:probb") + bDiscriminator("pfDeepFlavourJetTags:probbb") + bDiscriminator("pfDeepFlavourJetTags:problepb")) > 0.2770' # b tag with medium WP
        ),
        #filter = cms.bool(True)
)

bkgVeto = cms.EDFilter("Background_filter_QCD",
        muons = cms.InputTag("muonsForVeto"),
        electrons = cms.InputTag("slimmedElectrons"),
        electronId = cms.string("mvaEleID-Fall17-iso-V2-wpLoose"),
        bjets = cms.InputTag("bJetsForVeto")
)

Ntuplizer = cms.EDAnalyzer("NtuplizerQCD",
        treeName = cms.string("TagAndProbe"),
        taus = cms.InputTag("goodTaus"),
        jets = cms.InputTag("goodJets"),
        genParticles = cms.InputTag("genInfo"),
        triggerSet = cms.InputTag("patTriggerUnpacker"),
        triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
        Vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        puInfo = cms.InputTag("slimmedAddPileupInfo"),
        met = cms.InputTag("slimmedMETs"),
        triggerListTag = HLTLIST_TAG,
        triggerListProbe = HLTLIST,
        useHLTMatch = cms.bool(True),
        isMC = cms.bool(True),
        filterPath = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"),
        year = cms.int32(2018)
)

TAndPSeq = cms.Sequence(
    # tightJetIdLepVeto +
    slimmedJetsWithUserData +
    goodJets +
    genInfo
)

VetoSeq = cms.Sequence(
    muonsForVeto +
    bJetsForVeto +
    bkgVeto
)

NtupleSeq = cms.Sequence(
    patTriggerUnpacker +
    Ntuplizer
)
