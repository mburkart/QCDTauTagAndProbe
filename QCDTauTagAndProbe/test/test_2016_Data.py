import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import FWCore.PythonUtilities.LumiList as LumiList

process = cms.Process("TagAndProbe")

isMC = False

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

#### handling of cms line options for tier3 submission
#### the following are dummy defaults, so that one can normally use the config changing file list by hand etc.

options = VarParsing.VarParsing ('analysis')
options.register ('skipEvents',
                  -1, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Number of events to skip")
options.register ('JSONfile',
                  "", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "JSON file (empty for no JSON)")
if isMC:
    options.outputFile = 'NTuple_MC.root'
else:
    options.outputFile = 'NTuple_Data.root'
options.inputFiles = []

#options.register('numOfThreads',
#                 10,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.int,
#                 'Number of threads.'
#                )
#
#options.register('numOfStreams',
#                 10,
#                 VarParsing.VarParsing.multiplicity.singleton,
#                 VarParsing.VarParsing.varType.int,
#                 'Number of streams.'
#                )
options.parseArguments()
# START ELECTRON CUT BASED ID SECTION
#
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
#

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       era="2016-Legacy",
                       eleIDModules=[
                           "RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff",

                           "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff",
                           "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff",

                           "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff",
                           "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff",
                       ],
                       phoIDModules=[])

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
        process,
        jetSource = cms.InputTag("slimmedJets"),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag("slimmedSecondaryVertices"),
        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
        btagDiscriminators = [
            'pfDeepFlavourJetTags:probb',
            'pfDeepFlavourJetTags:probbb',
            'pfDeepFlavourJetTags:problepb',
            ],
        postfix='NewDFTraining'
)

process.bTaggingSequence = cms.Sequence(
        process.patJetCorrFactorsNewDFTraining +
        process.updatedPatJetsNewDFTraining +
        process.pfImpactParameterTagInfosNewDFTraining +
        process.pfInclusiveSecondaryVertexFinderTagInfosNewDFTraining +
        process.pfDeepCSVTagInfosNewDFTraining +
        process.pfDeepFlavourTagInfosNewDFTraining +
        process.pfDeepFlavourJetTagsNewDFTraining +
        process.patJetCorrFactorsTransientCorrectedNewDFTraining +
        process.updatedPatJetsTransientCorrectedNewDFTraining +
        process.selectedUpdatedPatJetsNewDFTraining
)

#START RERUNNING OF ID TRAINING
#
# set up the rerunning of the latest tau id trainings
import RecoTauTag.RecoTau.tools.runTauIdMVA as idemb
na = idemb.TauIDEmbedder(process, cms,
        debug=True,
        updatedTauName="NewTauIDsEmbedded",
        toKeep=["2017v2", "newDM2017v2", "deepTau2017v2p1"]
)
na.runTauID()

# Produce the quark gluon likelihood variable
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("RecoJets.JetProducers.QGTagger_cfi")
process.QGTagger.srcJets = cms.InputTag("slimmedJets")
process.QGTagger.jetsLabel = cms.string("QGL_AK4PFchs")


if not isMC: # will use 80X
    from Configuration.AlCa.autoCond import autoCond
    process.GlobalTag.globaltag = '102X_dataRun2_nanoAOD_2016_v1'
    process.load('QCDTauTagAndProbe.QCDTauTagAndProbe.tagAndProbe_2016_cff')
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
            '/store/data/Run2016E/JetHT/MINIAOD/17Jul2018-v1/00000/B25DDF93-0F8D-E811-A8B9-0017A4770470.root'
        ),
    )
else:
    process.GlobalTag.globaltag = '102X_mcRun2_asymptotic_v6' #MC 25 ns miniAODv2
    process.load('QCDTauTagAndProbe.QCDTauTagAndProbe.MCanalysis_2016_cff')
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
             '/store/mc/RunIISummer16MiniAODv3/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/100000/6EA405A4-9AE1-E811-A68D-48FD8E2824B3.root'
        )
    )


if options.JSONfile:
    print "Using JSON: " , options.JSONfile
    process.source.lumisToProcess = LumiList.LumiList(filename = options.JSONfile).getVLuminosityBlockRange()

if options.inputFiles:
    process.source.fileNames = cms.untracked.vstring(options.inputFiles)

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(-1)
            )

options.maxEvents = -1

if options.maxEvents >= -1:
    process.maxEvents.input = cms.untracked.int32(options.maxEvents)
if options.skipEvents >= 0:
    process.source.skipEvents = cms.untracked.uint32(options.skipEvents)



process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
#    numberOfThreads = cms.untracked.uint32(options.numOfThreads),
#    numberOfStreams = cms.untracked.uint32(options.numOfStreams)
)

process.p = cms.Path(
    process.hltFilter +
    process.rerunMvaIsolationSequence +
    process.NewTauIDsEmbedded +
    process.goodTaus +
    process.egammaPostRecoSeq +
    process.QGTagger +
    process.TAndPSeq +
    process.bTaggingSequence +
    process.VetoSeq +
    process.NtupleSeq
)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# Adding ntuplizer
process.TFileService=cms.Service('TFileService',fileName=cms.string(options.outputFile))
