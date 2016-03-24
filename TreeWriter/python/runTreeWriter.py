import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os,re
import getpass

def guessDatasetFromFileName( filename ):
    # This reproduces the dataset roughly if the file is on /store
    # Not reproduced are e.g. the pileup scenario
    # For local files, specify your own rules or run it with the 'dataset' option
    nParts=filename.split("/")
    if len(nParts)>6:
        return "/{}/{}-{}/{}".format(nParts[-5],nParts[-6],nParts[-3],nParts[-4])
    return filename

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

# defaults
options.inputFiles = 'root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/0CB41EBB-CA6D-E511-A28E-002618943C3A.root'
options.outputFile = 'photonTree.root'
options.maxEvents = -1
# get and parse the command line arguments
options.parseArguments()

dataset=options.dataset or guessDatasetFromFileName(options.inputFiles[0])
print "Assumed dataset:", dataset
isRealData=not dataset.endswith("SIM")

# the actual TreeWriter module
process = cms.Process("TreeWriter")
#process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100


# determine global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
gtName="76X_dataRun2_16Dec2015_v0" if isRealData else "76X_mcRun2_asymptotic_RunIIFall15DR76_v1"
if "74X" in dataset and not isRealData: "74X_mcRun2_asymptotic_v5"
process.GlobalTag = GlobalTag(process.GlobalTag, gtName, '')


######################
# PHOTONS, ELECTRONS #
######################
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
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


##########################
# Jet Energy Corrections #
##########################
# where .db files are placed (e.g. for JEC, JER)
# Crab will always be in the $CMSSW_BASE directory, so to run the code locally,
# a symbolic link is added
if not os.path.exists("src"): os.symlink(os.environ["CMSSW_BASE"]+"/src/", "src")
localDataBasePath='sqlite_file:src/TreeWriter/TreeWriter/data/'


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

######################
# Jets               #
######################
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *

process.jec = cms.ESSource("PoolDBESSource",
                           DBParameters = cms.PSet(messageLevel = cms.untracked.int32(0)),
                           timetype = cms.string('runnumber'),
                           toGet = cms.VPSet(
                               cms.PSet(
                                   record = cms.string('JetCorrectionsRecord'),
                                   tag    = (cms.string('JetCorrectorParametersCollection_Summer15_25nsV7_DATA_AK4PFchs') if isRealData
                                             else cms.string('JetCorrectorParametersCollection_Summer15_25nsV7_MC_AK4PFchs')
                                   ),
                                   label  = cms.untracked.string('AK4PFchs')
                               ),
                           ),
                           connect = (cms.string(localDataBasePath+'Summer15_25nsV7_DATA.db') if isRealData
                                      else cms.string(localDataBasePath+'Summer15_25nsV7_MC.db'))
)
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

######################
# MET Significance   #
######################
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMETSignificance
process.load("RecoMET/METProducers.METSignificance_cfi")
process.load("RecoMET/METProducers.METSignificanceParams_cfi")

process.load('Configuration.StandardSequences.Services_cff')
process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jer = cms.ESSource("PoolDBESSource",
                           CondDBSetup,
                           toGet = cms.VPSet(
                               # Pt Resolution
                               cms.PSet(
                                   record = cms.string('JetResolutionRcd'),
                                   tag    = cms.string('JR_MC_PtResolution_Summer15_25nsV6_AK4PFchs'),
                                   label  = cms.untracked.string('AK4PFchs_pt')
                               ),
                               # Phi Resolution
                               cms.PSet(
                                   record = cms.string('JetResolutionRcd'),
                                   tag    = cms.string('JR_MC_PhiResolution_Summer15_25nsV6_AK4PFchs'),
                                   label  = cms.untracked.string('AK4PFchs_phi')
                               ),
                               # Scale factors
                               cms.PSet(
                                   record = cms.string('JetResolutionScaleFactorRcd'),
                                   tag    = cms.string('JR_DATAMCSF_Summer15_25nsV6_AK4PFchs'),
                                   label  = cms.untracked.string('AK4PFchs')
                               ),
                           ),
                           connect = cms.string(localDataBasePath+'Summer15_25nsV6.db')
)
process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')

# rerun metcorrections and uncertainties
#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#runMetCorAndUncFromMiniAOD(process,isData=isRealData,jetColl="slimmedJets")


################################
# Define input and output      #
################################
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outputFile))


################################
# The actual TreeWriter module #
################################
process.TreeWriter = cms.EDAnalyzer('TreeWriter',
                                    # selection configuration
                                    HT_cut=cms.untracked.double(0),
                                    photon_pT_cut=cms.untracked.double(20), # for leading photon
                                    jet_pT_cut=cms.untracked.double(30), # for all jets
                                    isolatedPhotons=cms.untracked.bool(True), # for all photons in the collection
                                    minNumberPhotons_cut=cms.untracked.uint32(1),
                                    minNumberElectrons_cut=cms.untracked.uint32(0),
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
                                    metSig=cms.InputTag("METSignificance","METSignificance"),
                                    packedCandidates=cms.InputTag("packedPFCandidates"),
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
                                    metFilterNames=cms.untracked.vstring(
                                        "Flag_HBHENoiseFilter",
                                        "Flag_HBHENoiseIsoFilter",
                                        "Flag_CSCTightHaloFilter",
                                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                                        "Flag_goodVertices",
                                        "Flag_eeBadScFilter",
                                    ),
                                    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoWorstChargedIsolation"),
                                    pileupHistogramName=cms.untracked.string("pileupWeight_mix_2015_25ns_FallMC_matchData_PoissonOOTPU"),
                                    hardPUveto=cms.untracked.bool(False),
                                    # triggers to be saved
                                    # Warning: To be independent of the version number, the trigger result is saved if the trigger name begins
                                    # with the strings given here. E.g. "HLT" would always be true if any of the triggers fired.
                                    triggerNames=cms.vstring(),
                                    pfJetIDSelector=cms.PSet(version=cms.string('FIRSTDATA'), quality=cms.string('LOOSE'))
)

################################
# Modify the TreeWriter module #
################################

process.TreeWriter.hardPUveto=dataset.startswith("/QCD_HT100to200")
if "RunIISpring15MiniAODv2" in dataset: process.TreeWriter.pileupHistogramName="pileupWeight_mix_2015_25ns_Startup_PoissonOOTPU"

if "Fast" in dataset:
    process.TreeWriter.metFilterNames = [] # no met filters for fastsim
    process.TreeWriter.lheEventProduct = "source"
if dataset.endswith("3d7be4403ea17498be45eb057fcb0278/USER"): # miniAODv1
    process.TreeWriter.pileUpSummary = "addPileupInfo"


# determine user if not set by crab
user=options.user or getpass.getuser()
# user settings
if user=="kiesel":
    process.TreeWriter.HT_cut=500.
    process.TreeWriter.photon_pT_cut=90.
    process.TreeWriter.minNumberPhotons_cut=0
    process.TreeWriter.jet_pT_cut=40
    process.TreeWriter.triggerNames=[
        "HLT_Photon90_CaloIdL_PFHT500_v",
        "HLT_Photon90_v", #prescale: 90
        "HLT_PFHT600_v",
        "HLT_Photon175_v",
        "HLT_PFHT475_v", #prescale: 60
    ]
    if dataset.startswith("SingleElectron") or dataset.startswith("/DY"):
        process.TreeWriter.HT_cut=0.
        process.TreeWriter.photon_pT_cut=90.
        process.TreeWriter.minNumberPhotons_cut=0
        process.TreeWriter.minNumberElectrons_cut=1
        process.TreeWriter.triggerNames=[ "HLT_Ele23_WPLoose_Gsf_v" ]


elif user=="lange" or user=="jschulz":
    process.TreeWriter.jet_pT_cut=30.
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
        "HLT_Photon165_HE10_v",
        "HLT_Photon135_PFMET100_JetIdCleaned_v", # used in early data taking
        "HLT_Photon135_PFMET100_v",              # used in later data taking
        "HLT_Photon135_PFMET100_NoiseCleaned_v", # used for MC
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


####################
#     RUN          #
####################

process.p = cms.Path(
    process.METSignificance
    *process.photonIDValueMapProducer
    *process.egmGsfElectronIDSequence
    *process.egmPhotonIDSequence
    +cms.Sequence( process.patJetCorrFactorsReapplyJEC + process.patJetsReapplyJEC )
    *process.TreeWriter
)
