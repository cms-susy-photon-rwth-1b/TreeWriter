#!/usr/bin/env python2
from WMCore.Configuration import Configuration
import os
import subprocess

def searchUserDatasets( name ):
    return subprocess.check_output( 'das_client.py --limit=0 --query="dataset={} instance=prod/phys03"'.format(name), shell=True ).split("\n")[0:-1]

cmssw_src=os.environ['CMSSW_BASE']+'/src/'

config = Configuration()

config.section_("General")
config.General.requestName   = 'requestName'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName    = cmssw_src+'TreeWriter/TreeWriter/python/runTreeWriter.py'

config.section_("Data")
config.Data.inputDataset = '/SinglePhoton/Run2015C-PromptReco-v1/MINIAOD'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 500
config.Data.publication = False
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'V05'
config.Data.outLFNDirBase = "/store/user/kiesel/13TeV/nTuples/"

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_DE_RWTH'

datasets={}
datasets["T5Wg"] = [
    '/SMS-T5Wg_mGl-800to1000_mNLSP-0to950_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
    '/SMS-T5Wg_mGl-1150to1200_mNLSP-0to1150_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
    '/SMS-T5Wg_mGl-1350to1400_mNLSP-0to1350_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
    '/SMS-T5Wg_mGl-1550_mNLSP-0to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
]
datasets["T5gg"] = [
    '/SMS-T5gg_mGluino-1200-1350_mLSP-0to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
    '/SMS-T5gg_mGluino-1400-1550_mLSP-0to1400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
    '/SMS-T5gg_mGluino-1600-1750_mLSP-0to1600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
    '/SMS-T5gg_mGluino-1800-1850_mLSP-0to1700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
    '/SMS-T5gg_mGluino-2000_mLSP-0to1900_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
]

datasets["GJets_HT"] = [
    '/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
    '/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v4/MINIAODSIM',
    '/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
    '/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
    '/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
]
datasets["QCD_HT"] = [
    "/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM",
    # --- missing
    "/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM",
    "/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/MINIAODSIM",
    "/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v3/MINIAODSIM",
    "/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM",
]

datasets["ZJetsToNuNu_HT"] = [
    # --- missing
    "/ZJetsToNuNu_HT-600To800_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM",
    "/ZJetsToNuNu_HT-800To1200_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v3/MINIAODSIM",
    "/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM",
    "/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM",
]

datasets["WJetsToLNu_HT"] = [
    # --- missing
    "/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM",
    "/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM",
    "/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM",
    # --- missing
    "/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM",
    "/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM",
]

datasets["kiesel"] = [
    '/SinglePhoton/Run2015C_25ns-16Dec2015-v1/MINIAOD',
    '/SinglePhoton/Run2015D-16Dec2015-v1/MINIAOD',
    '/WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
    '/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM',
    '/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM',
    '/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM',
    '/JetHT/Run2015C_25ns-16Dec2015-v1/MINIAOD',
    '/SingleMuon/Run2015D-16Dec2015-v1/MINIAOD',
    '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM',
    '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
    '/SingleElectron/Run2015C_25ns-16Dec2015-v1/MINIAOD',
    '/SingleElectron/Run2015D-16Dec2015-v1/MINIAOD',
    '/SMS-T5Wg_mGl-1550_mNLSP-0to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-FastAsympt25ns_74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
] + datasets["GJets_HT"] + datasets["QCD_HT"] + datasets["ZJetsToNuNu_HT"] + datasets["WJetsToLNu_HT"]
#datasets["kiesel"] += datasets["T5Wg"]
#datasets["kiesel"] += datasets["T5gg"]
datasets["lange"] = [
    # data
    '/SinglePhoton/Run2016B-PromptReco-v2/MINIAOD',
    '/SingleElectron/Run2016B-PromptReco-v2/MINIAOD',
    '/MET/Run2016B-PromptReco-v2/MINIAOD',
    '/SingleMuon/Run2016B-PromptReco-v2/MINIAOD',
    '/MuonEG/Run2016B-PromptReco-v2/MINIAOD',
    '/JetHT/Run2016B-PromptReco-v2/MINIAOD',
    '/DoubleEG/Run2016B-PromptReco-v2/MINIAOD',
    '/DoubleMuon/Run2016B-PromptReco-v2/MINIAOD',
    # standard MC
    '/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
    '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
    # V+gamma
    '/WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
    '/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM'
    # V+gamma130
    '/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
    '/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
    '/ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCUETP8M1_13TeV-madgraph/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
    # diboson
    '/WWTo2L2Nu_13TeV-powheg/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
    '/WWToLNuQQ_13TeV-powheg/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
    '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
    '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
    # signal
    '/SMS-TChiWG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM',
] + datasets["GJets_HT"] + datasets["QCD_HT"] + datasets["ZJetsToNuNu_HT"] + datasets["WJetsToLNu_HT"]

datasets["jschulz"] = datasets["lange"]

datasets["rmeyer"] = [
    '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
]
#datasets["SMS-T5gg-private"] = searchUserDatasets( "/SMS-T5gg/kiesel-*/USER" )

# call with 'python crabConfig.py'
if __name__ == '__main__':
    import getpass
    from CRABAPI.RawCommand import crabCommand

    user=getpass.getuser()
    if user=="kiesel":
        config.Data.outputDatasetTag = 'v09'
        config.Data.outLFNDirBase = "/store/user/kiesel/13TeV/nTuples/"
    elif user=="lange":
        config.Data.outputDatasetTag = 'CMSSW_8_0_8_v01'
        config.Data.outLFNDirBase = "/store/user/jolange/run2/"
    elif user=="jschulz":
        config.Data.outputDatasetTag = 'CMSSW_8_0_8_v01'
        config.Data.outLFNDirBase = "/store/user/jschulz/run2/"
    elif user=="rmeyer":
        config.Data.outputDatasetTag = 'V01'
        config.Data.outLFNDirBase = "/store/user/rmeyer/13TeV/nTuples/"
    else:
        print "you shall not pass!"
        print "(unkown user '%s')"%user
        exit()

    for dataset in datasets[user]:

        isSim = 'SIM' in dataset
        isUser = dataset.endswith("/USER")
        isFastSim = "True" if "Fast" in dataset else "False"
        miniAODv="1" if dataset.endswith("3d7be4403ea17498be45eb057fcb0278/USER") else "2"

        config.Data.inputDBS = 'phys03' if isUser else 'global'

        if not isSim and not isUser:
            # https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2664.html
            # 2.6/fb
            config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-274443_13TeV_PromptReco_Collisions16_JSON.txt'
        else:
            try: del config.Data.lumiMask
            except: pass

        config.Data.splitting = 'FileBased' if isSim or isUser else 'LumiBased'
        config.Data.unitsPerJob = 5 if isSim or isUser else 40
        if isSim:
            config.General.requestName = dataset.split('/')[1]
        else:
            config.General.requestName = '_'.join(dataset.split('/')[1:3])
         # CRABClient.ClientExceptions.ConfigurationException: Invalid CRAB configuration: Parameter General.requestName should not be longer than 100 characters.
        config.General.requestName = config.General.requestName[:100]

        config.JobType.pyCfgParams = ["dataset="+dataset,"user="+user,"fastSim="+isFastSim,"miniAODv="+miniAODv]

        config.Data.inputDataset = dataset
        crabCommand('submit', config = config)
