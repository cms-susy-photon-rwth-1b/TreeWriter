**TreeWriter** to build a ROOT tree from MiniAOD. Photon Cut- and MVA-IDs are computed.

## Building and Running ##
Get CMSSW environment 80X

```
cmsrel CMSSW_8_0_25
cd CMSSW_8_0_25/src/
cmsenv
git cms-merge-topic shervin86:Moriond2017_JEC_energyScales
git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate
git cms-merge-topic cms-met:METRecipe_8020
git cms-merge-topic ikrav:egm_id_80X_v3_photons
git clone git@github.com:ECALELFS/ScalesSmearings.git EgammaAnalysis/ElectronTools/data/ScalesSmearings
git clone git@github.com:cms-susy-photon-rwth-1b/TreeWriter.git -b photonIdAndSmeaingTest
# get JEC and JER files
wget -qP TreeWriter/TreeWriter/data https://github.com/cms-jet/JECDatabase/raw/master/SQLiteFiles/Spring16_25nsFastSimV1_MC.db
wget -qP TreeWriter/TreeWriter/data https://github.com/cms-jet/JECDatabase/raw/master/SQLiteFiles/Summer16_23Sep2016V3_MC.db
wget -qP TreeWriter/TreeWriter/data https://github.com/cms-jet/JECDatabase/raw/master/SQLiteFiles/Summer16_23Sep2016AllV3_DATA.db
scram b -j7
cd TreeWriter
```
Create Pileup Histograms

```
make -C PUreweighting
```
Run the TreeWriter
- locally
```
voms-proxy-init -voms cms
cmsRun TreeWriter/python/runTreeWriter.py
```
- on the Grid using CRAB3
```
. /cvmfs/cms.cern.ch/crab3/crab.sh
cd crab
```
for a single dataset
```
crab submit -c crabConfig.py
```
for all datasets
```
python2 crabConfig.py
```

## Configure ##
in the python config, set
- `HT_cut`: minimum HT
- `photon_pT_cut`: minimum leading-photon pT

## Objects ##
### Photons ###
- official cut-based ID for Spring16 and Spring15 are stored [TWiki](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2)
- all photons are used. boolean flags for: loose/medium/tight

### Jets ###
- ak4PFJetsCHS
- all jets are used
- boolean flag for: loose
- boolean flags for whether a loose electron/photon is found within dR<0.4

### Muons ###
- fulfilling loose id
- tight id boolean flag
- relative isolation is stored

### Electrons ###
- fulfilling "veto" id
- boolean flags for loose/medium/tight
- recipes on [TWiki](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2)

### Generated Particles ###
- genJets collection is stored (= full slimmedGenJets)
- gen[Electrons|Photons]: status=1, pT>30
- Decay products (daughters) of W bosons
