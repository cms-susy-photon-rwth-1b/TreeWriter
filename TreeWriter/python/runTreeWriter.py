import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import re
import getpass

options = VarParsing ('analysis')
options.register ('dataset',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Name of the dataset, used to do further settings")
options.register ('user',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Name the user. If not set by crab, this script will determine it.")
options.register ('fastSim',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.bool,
                  "Whether the sample is simulated with fast-sim. (Default: False)")

# defaults
options.inputFiles = 'root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/0CB41EBB-CA6D-E511-A28E-002618943C3A.root'
options.outputFile = 'photonTree.root'
options.maxEvents = -1
options.fastSim=False
# get and parse the command line arguments
options.parseArguments()

# determine if Data or Simulation
isRealData=True
if options.dataset: # crab sumission
    isRealData=(not options.dataset.endswith("SIM"))
else: # running locally
    isRealData=("SIM" not in options.inputFiles[0])

# the actual TreeWriter module
process = cms.Process("TreeWriter")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.Geometry_cff")

# determine global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
gtName = "auto:run2_data" if isRealData else "auto:run2_mc"

# for further global tags, see here:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
if re.match( "/.*/Run2015.*05Oct2015.*/MINIAOD", options.dataset ) \
    or re.match( "/.*/Run2015D-PromptReco-v4/MINIAOD", options.dataset ): gtName = "74X_dataRun2_v5"
if re.match( "/.*/.*RunIISpring15MiniAODv2.*/MINIAODSIM", options.dataset ): gtName = "74X_mcRun2_asymptotic_v4"
process.GlobalTag = GlobalTag(process.GlobalTag, gtName, '')

hardPUveto=True if options.dataset.startswith("/QCD_HT100to200") else False

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( options.maxEvents ) )
process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(
    options.inputFiles
))


###############################
# Define MET Filters to apply #
###############################
# See https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2?rev=39
applyMetFilters=cms.untracked.vstring(
    "Flag_CSCTightHaloFilter",
    "Flag_eeBadScFilter"
)
# HBHE has to be manually re-run for early data.
# This is not applied as EDFilter, as suggested, but manually
# checked in TreeWriter (otherwise the "initial" event count
# is wrong)
# TODO: remove, when fixed upstream
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")


######################
# Jets               #
######################
pfJetIDSelector = cms.PSet(
    version = cms.string('FIRSTDATA'),
    quality = cms.string('LOOSE')
)

jecLevels = [ 'L1FastJet','L2Relative','L3Absolute' ]
if isRealData: jecLevels.append( 'L2L3Residual' )

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
  src = cms.InputTag("slimmedJets"),
  levels = jecLevels,
  payload = 'AK4PFchs'
)

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = patJetsUpdated.clone(
  jetSource = cms.InputTag("slimmedJets"),
  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
)

################################
# The actual TreeWriter module #
################################
process.TreeWriter = cms.EDAnalyzer('TreeWriter',
                                    # selection configuration
                                    HT_cut=cms.untracked.double(0),
                                    photon_pT_cut=cms.untracked.double(20),
                                    isolatedPhotons=cms.untracked.bool(True),
                                    # physics objects
                                    photons = cms.InputTag("slimmedPhotons"),
                                    jets = cms.InputTag("patJetsReapplyJEC"),
                                    muons = cms.InputTag("slimmedMuons"),
                                    genJets=cms.InputTag("slimmedGenJets"),
                                    electrons = cms.InputTag("slimmedElectrons"),
                                    mets = cms.InputTag("slimmedMETs"),
                                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    prunedGenParticles = cms.InputTag("prunedGenParticles"),
                                    pileUpSummary = cms.InputTag('slimmedAddPileupInfo'),
                                    lheEventProduct = cms.InputTag('externalLHEProducer'),
                                    # electron IDs
                                    electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                    electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                                    electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                                    electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                                    # photon IDs
                                    photonLooseIdMap   = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
                                    photonMediumIdMap  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
                                    photonTightIdMap   = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),
                                    photonMvaValuesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Values"),
                                    # met filters to apply
                                    metFilterNames=applyMetFilters,
                                    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoWorstChargedIsolation"),
                                    pileupHistogramName=cms.untracked.string( "pileupWeight_mix_2015_25ns_Startup_PoissonOOTPU" ),
                                    HBHENoiseFilterResult = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
                                    HBHEIsoNoiseFilterResult = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
                                    hardPUveto=cms.untracked.bool(hardPUveto),
                                    # triggers to be saved
                                    # Warning: To be independent of the version number, the trigger result is saved if the trigger name begins
                                    # with the strings given here. E.g. "HLT" would always be true if any of the triggers fired.
                                    triggerNames=cms.vstring(),
                                    pfJetIDSelector=pfJetIDSelector,
)

# determine user if not set by crab
user=options.user or getpass.getuser()

# user settings
if user=="kiesel":
    process.TreeWriter.HT_cut=500.
    process.TreeWriter.photon_pT_cut=90.
    process.TreeWriter.isolatedPhotons=False
    process.TreeWriter.triggerNames=[
        "HLT_Photon90_CaloIdL_PFHT500_v",
        "HLT_Photon90_v",
        "HLT_PFHT600_v",
        "HLT_Photon175_v",
    ]
elif user=="lange":
    process.TreeWriter.triggerNames=[
        "HLT_Photon90_CaloIdL_PFHT500_v",
        "HLT_Photon90_v",
        "HLT_PFHT600_v",
        'HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40_v',
        "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40_v",
        "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40_v",
        'HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40_v',
        'HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40_v',
        'HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40_v',
        "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_v",
        "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_v",
        "HLT_Photon22_v",
        "HLT_Photon30_v",
        "HLT_Photon36_v",
        "HLT_Photon50_v",
        "HLT_Photon36_R9Id90_HE10_IsoM_v",
        "HLT_Photon165_R9Id90_HE10_IsoM_v",
        "HLT_Photon135_PFMET100_v",
        'HLT_Photon135_PFMET100_NoiseCleaned_v',
        "HLT_Photon175_v",
        "HLT_Photon500_v",
        "HLT_PFMET170_v",
        "HLT_IsoMu18_v",
        "HLT_IsoMu20_v",
        "HLT_Mu20_v",
        "HLT_Mu50_v",
    ]
elif user=="rmeyer":
    process.TreeWriter.triggerNames=[
        "HLT_Photon90_CaloIdL_PFHT500_v",
        "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40_v",
        "HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15_v",
    ]
else:
    print "you shall not pass!"
    print "(unkown user '%s')"%user
    exit()

if options.fastSim:
    process.TreeWriter.metFilterNames = [] # no met filters for fastsim
    process.TreeWriter.pileUpSummary = "addPileupInfo" # for miniaod v1
elif re.match( "/SMS-.*/.*/USER", options.dataset ):
    # signal scan
    process.TreeWriter.metFilterNames = [] # no met filters for fastsim
    process.TreeWriter.pileUpSummary = "addPileupInfo" # for miniaod v1

process.TFileService = cms.Service("TFileService",fileName = cms.string(options.outputFile))


######################
# PHOTONS, ELECTRONS #
######################
from RecoEgamma.PhotonIdentification.PhotonIDValueMapProducer_cfi import *

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD

# turn on VID producer, indicate data format to be DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer  (process, dataFormat)

# define which IDs we want to produce
el_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']
ph_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff',
                 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff']

#add them to the VID producer
for idmod in el_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
for idmod in ph_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

####################
#     RUN          #
####################

process.p = cms.Path(
    process.photonIDValueMapProducer
    *process.egmGsfElectronIDSequence
    *process.egmPhotonIDSequence
    )
process.p += cms.Sequence( process.patJetCorrFactorsReapplyJEC + process.patJetsReapplyJEC )
if not options.fastSim:
    process.p*=process.HBHENoiseFilterResultProducer #produces HBHE bools (applied in TreeWriter manually)
process.p*=process.TreeWriter

