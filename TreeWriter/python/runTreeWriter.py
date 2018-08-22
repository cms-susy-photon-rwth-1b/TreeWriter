import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os,re
import getpass

def guessDatasetFromFileName(filename):
    # This reproduces the dataset roughly if the file is on /store
    # Not reproduced are e.g. the pileup scenario
    # For local files, specify your own rules or run it with the 'dataset' option
    nParts = filename.split("/")
    if "store" in nParts and len(nParts)>6:
        nParts = nParts[nParts.index("store"):]
        return "/{}/{}-{}/{}".format(nParts[3], nParts[2], nParts[5], nParts[4])
    if "user" in nParts:
        return nParts[-1].replace(".root", "")
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
                  
options.register("electronSmearing",
    "Moriond17_23Jan",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "correction type for electron energy smearing"
)
                  
                  

# defaults
#options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/100000/000786F7-3AD0-E611-A6AE-842B2B765E01.root'
#options.inputFiles = 'root://cms-xrd-global.cern.ch//store/data/Run2016B/MuonEG/MINIAOD/03Feb2017_ver2-v2/100000/008C5624-A1EC-E611-8238-0090FAA56F60.root'
options.inputFiles = 'root://cms-xrd-global.cern.ch//store/data/Run2016B/MET/MINIAOD/03Feb2017_ver1-v1/100000/1CDECD1B-0CEB-E611-A2A9-D4AE526A0455.root'
#options.inputFiles = 'root://cms-xrd-global.cern.ch//store/data/Run2016B/DoubleEG/MINIAOD/03Feb2017_ver2-v2/50000/001EB4EF-D3EA-E611-B94E-0CC47A4C8F26.root'
#options.inputFiles = 'root://cms-xrd-global.cern.ch//store/data/Run2016B/DoubleMuon/MINIAOD/03Feb2017_ver2-v2/100000/D628C213-0CEB-E611-B5E2-3417EBE7051F.root'
#options.inputFiles = 'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16MiniAODv2/SMS-TChiNG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/120000/040E9990-AA08-E711-BAAA-0025905B8574.root'
#options.inputFiles = 'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/00099D43-77ED-E611-8889-5065F381E1A1.root'
#options.inputFiles = '/store/mc/RunIISummer16MiniAODv2/GGM_GravitinoLSP_M1-200to1500_M2-200to1500_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSummer16Fast_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/000C44EA-F8EC-E711-8145-0242AC130002.root'
#options.inputFiles = '/store/mc/RunIISummer16MiniAODv2/GGM_GravitinoLSP_M1-50to1500_M3-1000to2500_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSummer16Fast_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/06B8CB51-31EB-E711-B445-0025905A6136.root'
#options.inputFiles = '/store/mc/RunIISummer16MiniAODv2/SMS-T5bbbbZg_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSummer16Fast_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/10000/00BA5D3D-2592-E711-82AE-0242AC110011.root'
#options.inputFiles = '/store/mc/RunIISummer16MiniAODv2/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0030B9D6-72C1-E611-AE49-02163E00E602.root'
#options.inputFiles = '/store/mc/RunIISummer16MiniAODv2/TTGamma_Dilept_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/80000/02915ED5-FBE5-E611-8C52-001E67586A2F.root'


options.outputFile = 'photonTree.root'
options.maxEvents = -1
#options.maxEvents = 1000
#options.maxEvents = 10
# get and parse the command line arguments
options.parseArguments()

dataset=options.dataset or guessDatasetFromFileName(options.inputFiles[0])
print "Assumed dataset:", dataset
isRealData=not dataset.endswith("SIM")

isSignal=False
useHTTrigger=True



electronCollection = cms.InputTag("slimmedElectrons", "", "PAT")
photonCollection   = cms.InputTag("slimmedPhotons", "", "PAT")


# the actual TreeWriter module
process = cms.Process("TreeWriter")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100



process.load('Geometry.CaloEventSetup.CaloTopology_cfi')
process.load("Configuration.StandardSequences.GeometryDB_cff")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# determine global tag
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
if isRealData:
    if "Run2016H" in dataset:
        process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v16"
    else:
        process.GlobalTag.globaltag = "80X_dataRun2_2016SeptRepro_v7"
else:    
    if "80X_mcRun2_asymptotic_2016_TrancheIV_v6" in dataset:
        process.GlobalTag.globaltag = "80X_mcRun2_asymptotic_2016_TrancheIV_v8"
    elif "80X_mcRun2_asymptotic_2016_miniAODv2" in dataset:
        process.GlobalTag.globaltag = "80X_mcRun2_asymptotic_2016_miniAODv2_v1"
    else:
        print "Could not guess correct global tag for", dataset




seq = cms.Sequence()

######################
# PHOTONS, ELECTRONS #
######################


from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)
process.load("EgammaAnalysis.ElectronTools.regressionApplication_cff")
seq += process.regressionApplication

# set the electron and photon sources
process.slimmedElectrons.src = electronCollection
process.slimmedPhotons.src = photonCollection

# overwrite output collections
electronCollection = cms.InputTag("slimmedElectrons", "", process.name_())
photonCollection = cms.InputTag("slimmedPhotons", "", process.name_())


process.selectedElectrons = cms.EDFilter("PATElectronSelector",
  src = electronCollection,
  cut = cms.string("pt>5 && abs(eta)<2.5")
)
electronCollection = cms.InputTag("selectedElectrons", "", process.name_())

process.selectedPhotons = cms.EDFilter("PATPhotonSelector",
  src = photonCollection,
  cut = cms.string("pt>5 && abs(eta)<2.5")
)
photonCollection = cms.InputTag("selectedPhotons", "", process.name_())


# setup the smearing
process.load("EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi")
from EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi import files
process.calibratedPatElectrons.isMC           = cms.bool(not isRealData)
process.calibratedPatElectrons.correctionFile = cms.string(files[options.electronSmearing])
process.calibratedPatElectrons.electrons      = electronCollection
seq += process.calibratedPatElectrons

process.load("EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi")
from EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi import files
process.calibratedPatPhotons.isMC           = cms.bool(not isRealData)
process.calibratedPatPhotons.correctionFile = cms.string(files[options.electronSmearing])
process.calibratedPatPhotons.photons      = photonCollection
seq += process.calibratedPatPhotons


process.load("Configuration.StandardSequences.Services_cff")
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(81),
        engineName  = cms.untracked.string("TRandom3")
    ),
    calibratedPatPhotons = cms.PSet(
        initialSeed = cms.untracked.uint32(81),
        engineName  = cms.untracked.string("TRandom3")
    )
)

# overwrite output collections
electronCollection = cms.InputTag("calibratedPatElectrons", "", process.name_())
photonCollection = cms.InputTag("calibratedPatPhotons", "", process.name_())

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

eleVIDModules = [
    #"RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff",
    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff",
    "RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff"
]

phoVIDModules = [
    'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff',
    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff'
]

switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
switchOnVIDPhotonIdProducer  (process, DataFormat.MiniAOD)

for mod in eleVIDModules:
    setupAllVIDIdsInModule(process, mod, setupVIDElectronSelection)

for idmod in phoVIDModules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

# update some VID modules to work with potentially changed electron collections
process.egmGsfElectronIDs.physicsObjectSrc = electronCollection
process.electronRegressionValueMapProducer.srcMiniAOD = electronCollection
process.electronMVAValueMapProducer.srcMiniAOD = electronCollection


process.egmPhotonIDs.physicsObjectSrc = photonCollection
process.egmPhotonIsolation.srcToIsolate = photonCollection
process.photonIDValueMapProducer.srcMiniAOD = photonCollection
process.photonRegressionValueMapProducer.srcMiniAOD = photonCollection
process.photonMVAValueMapProducer.srcMiniAOD = photonCollection


##########################
# Jet Energy Corrections #
##########################
# where .db files are placed (e.g. for JEC, JER)
# Crab will always be in the $CMSSW_BASE directory, so to run the code locally,
# a symbolic link is added
if not os.path.exists("src"): os.symlink(os.environ["CMSSW_BASE"]+"/src/", "src")
if "Fast" in dataset:
    dbPath = 'sqlite:'
    isSignal = True

    from CondCore.CondDB.CondDB_cfi import CondDB
    CondDB.__delattr__('connect')

    process.jec = cms.ESSource('PoolDBESSource',
        CondDB,
        connect = cms.string(dbPath+'Spring16_25nsFastSimV1_MC.db'),
        toGet = cms.VPSet(
            cms.PSet(
                record = cms.string('JetCorrectionsRecord'),
                tag    = cms.string('JetCorrectorParametersCollection_Spring16_25nsFastSimMC_V1_AK4PFchs'),
                label  = cms.untracked.string('AK4PFchs')
            )
        )
    )
    # Add an ESPrefer to override JEC that might be available from the global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')


jecLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
if isRealData: jecLevels.append('L2L3Residual')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(jecLevels), 'None')
)

##########################
# Top Pt reweighting     #
##########################
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
#process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
#
#process.decaySubset.fillMode = cms.string("kME")
#process.TreeWriter.ttGenEvent = cms.InputTag('genEvt')

##########################
# MET                    #
##########################
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(
    process,
    isData=isRealData,
)


#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
    #ignoreTotal = cms.untracked.int32(1)
#)



################################
# MET Filter                   #
################################
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")


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
                                    isSignalBoolean=cms.untracked.bool(isSignal),
                                    HT_cut=cms.untracked.double(0),
                                    photon_pT_cut=cms.untracked.double(5), # for leading photon
                                    jet_pT_cut=cms.untracked.double(30), # for all jets
                                    isolatedPhotons=cms.untracked.bool(True), # for all photons in the collection
                                    minNumberPhotons_cut=cms.untracked.uint32(0),
                                    minNumberElectrons_cut=cms.untracked.uint32(0),
                                    minNumberLeptons_cut=cms.untracked.uint32(0),
                                    minNumberBinos_cut=cms.untracked.uint32(0),
                                    # physics objects
                                    photons = photonCollection,
                                    jets = cms.InputTag("updatedPatJetsUpdatedJEC"),
                                    muons = cms.InputTag("slimmedMuons"),
                                    genJets= cms.InputTag("slimmedGenJets"),
                                    electrons = electronCollection,
                                    mets = cms.InputTag("slimmedMETs", "", "TreeWriter"),
                                    metCorr = cms.InputTag(""),
                                    metCorrCal = cms.InputTag(""),
                                    caloMets = cms.InputTag("slimmedMETs"),
                                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                    ebRecHits = cms.InputTag("reducedEgamma","reducedEBRecHits"),
                                    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                    prunedGenParticles = cms.InputTag("prunedGenParticles"),
                                    pileUpSummary = cms.InputTag('slimmedAddPileupInfo'),
                                    lheEventProduct = cms.InputTag('externalLHEProducer'),
                                    packedCandidates=cms.InputTag("packedPFCandidates"),
                                    
                                    #ttGenEvent = cms.InputTag("genEvt"),
                                    
                                    # electron IDs
                                    electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                                    electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                                    electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                                    electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
 
                                    # ID decisions (common to all formats)
                                    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
                                    eleTightIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
                                    # ValueMaps with MVA results
                                    mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
                                    mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
                                    
                                    
                                    beamSpot = cms.InputTag('offlineBeamSpot'),
                                    conversionsMiniAOD = cms.InputTag('reducedEgamma:reducedConversions'),
                                    
                                    # photon IDs
                                    photonLooseIdMap   = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-loose"),
                                    photonMediumIdMap  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium"),
                                    photonTightIdMap   = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-tight"),

                                    photonMediumIdBoolMap_mva = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp90"),
                                    photonMediumIdFullInfoMap_mva = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp90"),
                                    photonMvaValuesMap     = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),
                                    photonMvaCategoriesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring16NonTrigV1Categories"),
                                    
                                    # met filters to apply
                                    metFilterNames=cms.untracked.vstring(
                                        "Flag_HBHENoiseFilter",
                                        "Flag_HBHENoiseIsoFilter",
                                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                                        "Flag_goodVertices",
                                        "Flag_eeBadScFilter",
                                        "Flag_globalTightHalo2016Filter",
                                    ),
                                    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoWorstChargedIsolation"),
                                    pileupHistogramName=cms.untracked.string("pileupWeight_mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU"),
                                    hardPUveto=cms.untracked.bool(False),
                                    reMiniAOD=cms.untracked.bool(False),
                                    # triggers to be saved
                                    # Warning: To be independent of the version number, the trigger result is saved if the trigger name begins
                                    # with the strings given here. E.g. "HLT" would always be true if any of the triggers fired.
                                    triggerNames=cms.vstring(),
                                    triggerObjectNames=cms.vstring(),
                                    pfJetIDSelector=cms.PSet(version=cms.string('FIRSTDATA'), quality=cms.string('LOOSE')),
                                    triggerPrescales=cms.vstring(), # also useful to check whether a trigger was run
                                    BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter"),
                                    BadPFMuonFilter = cms.InputTag("BadPFMuonFilter"),
                                    metCorrected = cms.InputTag("slimmedMETs"),
                                    metCalibrated = cms.InputTag("slimmedMETs")
)

#process.TreeWriter.ttGenEvent = cms.InputTag('genEvt')
#process.TreeWriter.ttGenEvent = cms.InputTag('prunedGenParticles')


################################
# Modify the TreeWriter module #
################################

process.TreeWriter.hardPUveto=dataset.startswith("/QCD_HT100to200")

if "03Feb2017" in dataset:
    process.TreeWriter.reMiniAOD = True
    process.TreeWriter.mets = cms.InputTag("slimmedMETsMuEGClean", "", "PAT")
    process.TreeWriter.metCorrected = cms.InputTag("slimmedMETsMuEGClean", "", "TreeWriter")
    process.TreeWriter.metCalibrated = cms.InputTag("slimmedMETsMuEGCleanCalibrated", "", "TreeWriter")
    process.TreeWriter.metFilterNames.extend(["Flag_chargedHadronTrackResolutionFilter", "Flag_muonBadTrackFilter"])

    # Now you are creating the e/g corrected MET on top of the bad muon corrected MET (on re-miniaod)
    from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG
    corMETFromMuonAndEG(
        process,
        pfCandCollection = "",
        electronCollection = "slimmedElectronsBeforeGSFix",
        photonCollection = "slimmedPhotonsBeforeGSFix",
        corElectronCollection = "slimmedElectrons",
        corPhotonCollection = "slimmedPhotons",
        allMETEGCorrected = True,
        muCorrection = False,
        eGCorrection = True,
        runOnMiniAOD = True,
        postfix = "MuEGClean"
    )
    process.slimmedMETsMuEGClean = process.slimmedMETs.clone()
    process.slimmedMETsMuEGClean.src = cms.InputTag("patPFMetT1MuEGClean")
    process.slimmedMETsMuEGClean.rawVariation =  cms.InputTag("patPFMetRawMuEGClean")
    process.slimmedMETsMuEGClean.t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
    del process.slimmedMETsMuEGClean.caloMET

    corMETFromMuonAndEG(
        process,
        pfCandCollection = "",
        electronCollection = "slimmedElectronsBeforeGSFix",
        photonCollection = "slimmedPhotonsBeforeGSFix",
        #corElectronCollection = "calibratedPatElectrons",
        #corPhotonCollection = "calibratedPatPhotons",
        corElectronCollection = electronCollection.value(),
        corPhotonCollection = photonCollection.value(),
        allMETEGCorrected = True,
        muCorrection = False,
        eGCorrection = True,
        runOnMiniAOD = True,
        postfix = "MuEGCleanCalibrated"
    )
    process.slimmedMETsMuEGCleanCalibrated = process.slimmedMETs.clone()
    process.slimmedMETsMuEGCleanCalibrated.src = cms.InputTag("patPFMetT1MuEGClean")
    process.slimmedMETsMuEGCleanCalibrated.rawVariation =  cms.InputTag("patPFMetRawMuEGClean")
    process.slimmedMETsMuEGCleanCalibrated.t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
    del process.slimmedMETsMuEGCleanCalibrated.caloMET

if not isRealData:
    process.TreeWriter.metFilterNames.remove("Flag_eeBadScFilter")
if "Fast" in dataset:
    process.TreeWriter.metFilterNames.remove("Flag_globalTightHalo2016Filter")
    process.TreeWriter.lheEventProduct = "source"
    if "T5Wg" in dataset or "T6Wg" in dataset:
        process.TreeWriter.minNumberBinos_cut = 1

if "PUMoriond17" in dataset:
    process.TreeWriter.pileupHistogramName=cms.untracked.string("pileupWeight_mix_2016_25ns_Moriond17MC_PoissonOOTPU")
if "PUSummer16" in dataset:
    process.TreeWriter.pileupHistogramName=cms.untracked.string("pileupWeight_mix_2016_25ns_Moriond17MC_PoissonOOTPU")
if "GGM_GravitinoLSP_M1" in dataset:
    process.TreeWriter.pileupHistogramName=cms.untracked.string("pileupWeight_mix_2016_25ns_Moriond17MC_PoissonOOTPU")

# determine user if not set by crab
user=options.user or getpass.getuser()
# user settings
if user=="kiesel":
    process.TreeWriter.HT_cut=500.
    process.TreeWriter.photon_pT_cut=90.
    process.TreeWriter.minNumberPhotons_cut=0
    process.TreeWriter.triggerObjectNames = ["hltEG90CaloIdLHEFilter"]
    process.TreeWriter.triggerNames=[
        "HLT_Photon90_CaloIdL_PFHT600_v",
        "HLT_Photon90_v",
        "HLT_PFHT600_v",
        "HLT_PFHT800_v",
        "HLT_Ele27_eta2p1_WPLoose_Gsf_v",
        "HLT_Ele27_eta2p1_WPTight_Gsf_v",
        "HLT_PFJet450_v",
        "HLT_AK8PFJet450_v",
    ]
    process.TreeWriter.triggerPrescales=process.TreeWriter.triggerNames
    if "SingleElectron" in dataset or "DY" in dataset:
        process.TreeWriter.triggerObjectNames = ["hltEle27erWPLooseGsfTrackIsoFilter", "hltEle27erWPTightGsfTrackIsoFilter", "hltEG90CaloIdLHEFilter"]
        process.TreeWriter.triggerNames = ["HLT_Ele27_eta2p1_WPLoose_Gsf_v", "HLT_Ele27_eta2p1_WPTight_Gsf_v"]
        process.TreeWriter.HT_cut = 0.
        process.TreeWriter.photon_pT_cut = 25.
        process.TreeWriter.minNumberPhotons_cut = 1
    if "Fast" in dataset: # signal scan
        process.TreeWriter.HT_cut = 0.

elif user=="jschulz" or user=="dmeuser":
    process.TreeWriter.photon_pT_cut=100
    process.TreeWriter.triggerObjectNames = ["hltEG90CaloIdLHEFilter", "hltEG165HE10Filter"]
    if "Fast" in dataset:
        process.TreeWriter.minNumberPhotons_cut=0
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
        #Lepton
        "HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v",
        "HLT_Mu17_Photon30_CaloIdL_L1ISO_v",
        "HLT_Mu38NoFilterNoVtx_Photon38_CaloIdL",
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
elif user=="swuchterl":
    process.TreeWriter.minNumberPhotons_cut=0
    process.TreeWriter.minNumberLeptons_cut=0
    process.TreeWriter.HT_cut=0
    process.TreeWriter.photon_pT_cut=5 # for leading photon
    process.TreeWriter.jet_pT_cut=30
    process.TreeWriter.triggerObjectNames = [
       #DoubleEle
        "hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter",#1
        "hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter",#2
        "hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter",#3
        "hltDiEle33CaloIdLGsfTrkIdVLMWPMS2UnseededFilter",#4
        #DoubleMu
        "hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4",#5
        "hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4",#6
        "hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2",#7
        "hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2",#8
        "hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4DzFiltered0p2",#9
        "hltDiMuonGlb27Trk8DzFiltered0p2",#10
        "hltDiMuonGlb30Trk11DzFiltered0p2",#11
        #MuEle
        "hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",#12.1
        "hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17",#12.2
        "hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",#13.1
        "hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23",#13.2
        "hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLDZFilter",#14
        "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",#15.1
        "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23",#15.2
        "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter",#16
        "hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",#17.1
        "hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8",#17.2
        "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",#18.1
        "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8",#18.2
        "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter",#19
        "hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",#20.1
        "hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered12",#20.2
        "hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter",#21
        "hltEle30CaloIdLGsfTrkIdVLDPhiUnseededFilter",#22.1
        "hltL3fL1sMu22orMu25orMu20EG15orMu5EG20L1f0L2f10QL3Filtered30Q",#22.2
        "hltEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter",#23.1
        "hltL3fL1sMu22orMu25orMu20EG15orMu5EG20L1f0L2f10QL3Filtered33Q"#23.2
    ]
    triggerNames=[
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",
        "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v",
        #DoubleMu
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        "HLT_Mu27_TkMu8_v",
        "HLT_Mu30_TkMu11_v",
        #MuEle
        "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
        "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v",
        "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v"        
    ]
    triggerNamesHT=[
        "HLT_PFHT200_v",
        "HLT_PFHT250_v",
        "HLT_PFHT300_v",
        "HLT_PFHT350_v",
        "HLT_PFHT400_v",
        "HLT_PFHT475_v",
        "HLT_PFHT600_v",
        "HLT_PFHT650_v",
        "HLT_PFHT800_v"
    ]
    triggerNamesMET=[
        "HLT_PFMET110_PFMHT110_IDTight_v",
        "HLT_PFMET120_PFMHT120_IDTight_v",
        "HLT_PFMET170_NoiseCleaned_v",
        "HLT_PFMET170_HBHECleaned_v",
        "HLT_PFMET170_JetIdCleaned_v",
        "HLT_PFMET170_NotCleaned_v",
        "HLT_PFMET300_v",
        "HLT_PFMET400_v",
        "HLT_PFMET500_v",
        "HLT_PFMET600_v"
    ]
    if (useHTTrigger):
        process.TreeWriter.triggerNames=triggerNames+triggerNamesHT+triggerNamesMET
    else:
        process.TreeWriter.triggerNames=triggerNames
    
    #process.TreeWriter.triggerNames=["HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
        #"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
        #"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
        #"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",
        #"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
        #"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        #"HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        #"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
        #"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",
        #"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v",
        #"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
        #"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
        #"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",
        #"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
        #"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
        #"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
        #"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
        #"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",
        #"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v",
        #"HLT_Mu27_TkMu8_v",
        #"HLT_Mu30_TkMu11_v",
        #"HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v",
        #"HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v",
        #"HLT_PFHT200_v",
        #"HLT_PFHT250_v",
        #"HLT_PFHT300_v",
        #"HLT_PFHT350_v",
        #"HLT_PFHT400_v",
        #"HLT_PFHT475_v",
        #"HLT_PFHT600_v",
        #"HLT_PFHT650_v",
        #"HLT_PFHT800"
    #]
else:
    print "you shall not pass!"
    print "(unkown user '%s')"%user
    exit()

for trig in process.TreeWriter.triggerPrescales:
    assert(trig in process.TreeWriter.triggerNames),"Trigger '"+trig+"' is not used, so prescale cannot be stored!"

####################
#     RUN          #
####################

############process.p = cms.Path(process.BadPFMuonFilter*process.BadChargedCandidateFilter*process.TreeWriter#)
process.p = cms.Path(process.BadPFMuonFilter + process.BadChargedCandidateFilter + seq + process.TreeWriter)
#process.p = cms.Path(process.BadPFMuonFilter + process.BadChargedCandidateFilter + seq +  process.makeGenEvt + process.TreeWriter)
############process.p = cms.Path(process.regressionApplication * process.calibratedPatElectrons * process.BadPFMuonFilter * process.BadChargedCandidateFilter  * process.TreeWriter)
