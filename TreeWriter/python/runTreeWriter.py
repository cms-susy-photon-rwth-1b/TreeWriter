import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os,re
import getpass

cmssw_src=os.environ['CMSSW_BASE']+'/src/'

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
options.register ('miniAODv',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.int,
                  "The MiniAOD version. (Default: 2)")

# defaults
options.inputFiles = 'root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/0CB41EBB-CA6D-E511-A28E-002618943C3A.root'
options.outputFile = 'photonTree.root'
options.maxEvents = -1
options.fastSim=False
options.miniAODv=2
# get and parse the command line arguments
options.parseArguments()

isCrabSubmission=bool(options.dataset) # only set for crab sumission

# determine if Data or Simulation
isRealData=True
if isCrabSubmission:
    isRealData=(not options.dataset.endswith("SIM"))
else: # running locally
    isRealData=("SIM" not in options.inputFiles[0])

options.fastSim = (options.fastSim
                   or bool(re.match( "/SMS-.*/.*/USER", options.dataset )) # signal scan
                   )

# where .db files are placed (e.g. for JEC, JER)
localDataBasePath=('sqlite_file:src/TreeWriter/TreeWriter/data/' if isCrabSubmission
                   else 'sqlite_file:'+cmssw_src+'/TreeWriter/TreeWriter/data/')

# the actual TreeWriter module
process = cms.Process("TreeWriter")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

# determine global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
gtName = "auto:run2_data" if isRealData else "auto:run2_mc"
# for further global tags, see here:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
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
# See https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#MiniAOD_76X_v2_produced_with_the
applyMetFilters=cms.untracked.vstring(
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_CSCTightHalo2015Filter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_goodVertices",
    "Flag_eeBadScFilter",
)

######################
# Jets               #
######################
pfJetIDSelector = cms.PSet(
    version = cms.string('FIRSTDATA'),
    quality = cms.string('LOOSE')
)


######################
# MET Significance   #
######################
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMETSignificance
process.load("RecoMET/METProducers.METSignificance_cfi")
process.load("RecoMET/METProducers.METSignificanceParams_cfi")

# rerun metcorrections and uncertainties
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,isData=isRealData)

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
                                    jets = cms.InputTag("slimmedJets"),
                                    muons = cms.InputTag("slimmedMuons"),
                                    genJets=cms.InputTag("slimmedGenJets"),
                                    electrons = cms.InputTag("slimmedElectrons"),
                                    mets = cms.InputTag("slimmedMETs","","TreeWriter"),
                                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    prunedGenParticles = cms.InputTag("prunedGenParticles"),
                                    pileUpSummary = cms.InputTag('slimmedAddPileupInfo'),
                                    lheEventProduct = cms.InputTag('externalLHEProducer'),
                                    metSig=cms.InputTag("METSignificance","METSignificance","TreeWriter"),
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
                                    metFilterNames=applyMetFilters,
                                    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoWorstChargedIsolation"),
                                    pileupHistogramName=cms.untracked.string("pileupWeight_mix_2015_25ns_FallMC_matchData_PoissonOOTPU"),
                                    hardPUveto=cms.untracked.bool(hardPUveto),
                                    # triggers to be saved
                                    # Warning: To be independent of the version number, the trigger result is saved if the trigger name begins
                                    # with the strings given here. E.g. "HLT" would always be true if any of the triggers fired.
                                    triggerNames=cms.vstring(),
                                    triggerPrescales=cms.vstring(), # also useful to check whether a trigger was run
                                    # if trigger objects are stored
                                    storeTriggerObjects=cms.untracked.bool(True),
                                    pfJetIDSelector=pfJetIDSelector,
)

# determine user if not set by crab
user=options.user or getpass.getuser()

# check for 74X samples and use right PU histogram
if (isCrabSubmission and "RunIISpring15MiniAODv2" in options.dataset) \
   or (len (options.inputFiles)==1 and '/ggm/' in options.inputFiles[0]):
    process.TreeWriter.pileupHistogramName="pileupWeight_mix_2015_25ns_Startup_PoissonOOTPU"

# 2016 pileup
if (isCrabSubmission and "RunIISpring16MiniAODv2-PUSpring16" in options.dataset):
    process.TreeWriter.pileupHistogramName="pileupWeight_mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU"

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
        "HLT_Ele23_WP75_Gsf_v",
        "HLT_Ele22_eta2p1_WP75_Gsf_v",
        "HLT_Ele22_eta2p1_WPLoose_Gsf_v",
        "HLT_Ele27_eta2p1_WPLoose_Gsf_v",
    ]
    if "SingleElectron" in options.dataset or "DY" in options.dataset:
        process.TreeWriter.HT_cut=0.
        process.TreeWriter.minNumberElectrons_cut=1

elif user=="lange" or user=="jschulz":
    process.TreeWriter.jet_pT_cut=30.
    process.TreeWriter.photon_pT_cut=100
    process.TreeWriter.storeTriggerObjects=False
    process.TreeWriter.triggerNames=[
        "HLT_Photon90_CaloIdL_PFHT500_v",
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
        "HLT_Photon75_v",
        "HLT_Photon90_v",
        "HLT_Photon120_v",
        "HLT_Photon36_R9Id90_HE10_IsoM_v",
        "HLT_Photon165_R9Id90_HE10_IsoM_v",
        "HLT_Photon165_HE10_v",
        "HLT_Photon135_PFMET100_JetIdCleaned_v", # used in early data taking
        "HLT_Photon135_PFMET100_v",              # used in later data taking
        "HLT_Photon135_PFMET100_NoiseCleaned_v", # used for MC
        "HLT_Photon175_v",
        "HLT_Photon500_v",
        "HLT_PFMET170_NoiseCleaned_v",
        "HLT_PFMET170_HBHECleaned_v",
        "HLT_PFMET170_JetIdCleaned_v",
        "HLT_PFMET170_NotCleaned_v",
        "HLT_IsoMu18_v",
        "HLT_IsoMu20_v",
        "HLT_IsoMu22_v",
        "HLT_Mu20_v",
        "HLT_Mu45_eta2p1_v",
        "HLT_Mu50_v",
        "HLT_Mu30_TkMu11_v",
        "HLT_DoubleIsoMu17_eta2p1_v",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
        "HLT_Mu17_Photon22_CaloIdL_L1ISO_v",
        "HLT_Ele27_eta2p1_WPTight_Gsf_v",
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
         # HT
        "HLT_PFHT125_v",
        "HLT_PFHT200_v",
        "HLT_PFHT250_v",
        "HLT_PFHT300_v",
        "HLT_PFHT350_v",
        "HLT_PFHT400_v",
        "HLT_PFHT475_v",
        "HLT_PFHT600_v",
        "HLT_PFHT650_v",
        "HLT_PFHT800_v",
    ]
    process.TreeWriter.triggerPrescales=[
        "HLT_Photon135_PFMET100_JetIdCleaned_v", # used in early data taking
        "HLT_Photon135_PFMET100_v",              # used in later data taking
        "HLT_Photon135_PFMET100_NoiseCleaned_v", # used for MC
        "HLT_Photon22_v",
        "HLT_Photon30_v",
        "HLT_Photon36_v",
        "HLT_Photon50_v",
        "HLT_Photon75_v",
        "HLT_Photon90_v",
        "HLT_Photon120_v",
        "HLT_Photon36_R9Id90_HE10_IsoM_v",
        "HLT_IsoMu18_v",
        "HLT_Mu20_v",
        "HLT_Mu17_Photon22_CaloIdL_L1ISO_v",
        # HT
        "HLT_PFHT125_v",
        "HLT_PFHT200_v",
        "HLT_PFHT250_v",
        "HLT_PFHT300_v",
        "HLT_PFHT350_v",
        "HLT_PFHT400_v",
        "HLT_PFHT475_v",
        "HLT_PFHT600_v",
        "HLT_PFHT650_v",
        "HLT_PFHT800_v",
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

for trig in process.TreeWriter.triggerPrescales:
    assert(trig in process.TreeWriter.triggerNames),"Trigger '"+trig+"' is not used, so prescale cannot be stored!"

if options.fastSim:
    process.TreeWriter.metFilterNames = [] # no met filters for fastsim
    process.TreeWriter.lheEventProduct = "source"
if options.miniAODv==1:
    process.TreeWriter.pileUpSummary = "addPileupInfo"

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
    process.METSignificance
    *process.photonIDValueMapProducer
    *process.egmGsfElectronIDSequence
    *process.egmPhotonIDSequence
    )
process.p*=process.TreeWriter

print "#########################################################"
print "This is what I think I am processing..."
print "  user           ",user
print "  CRAB submission",isCrabSubmission
print "  dataset        ",options.dataset
print "  isRealData     ",isRealData
print "  Global Tag     ",gtName
print "  fastSim        ",options.fastSim
print "  MiniAODv       ",options.miniAODv
print "  hardPUveto     ",hardPUveto
print "  PU             ",process.TreeWriter.pileupHistogramName.pythonValue()
print "#########################################################"
