**TreeWriter** to build a ROOT tree from MiniAOD. Photon Cut- and MVA-IDs are computed.

## Building and Running ##
Get CMSSW environment 76X

```
cmsrel CMSSW_7_6_3_patch2
cd CMSSW_7_6_3_patch2/src/
cmsenv
```
get MET Significance Recipe

```
git-cms-merge-topic -u cms-met:76X-METSignificance-patch0
```
Get and build the TreeWriter

```
git clone git@github.com:cms-susy-photon-rwth-1b/TreeWriter.git
scram b -j5
cd TreeWriter
```
Create Pilup Histograms

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

## Input files ##
located in `TreeWriter/data`:
- `Summer15_25nsV6.db` for JER (used for MET Significance), taken from [cms-met](https://github.com/cms-met/cmssw/blob/f0ac9b3e56e85d03c8dbe6e5cb101274fb356520/RecoMET/METProducers/test/Summer15_25nsV6.db) which is [this](https://github.com/cms-jet/JRDatabase/blob/aa321717d57773d074b5d328c5e71d473e7cf836/SQLiteFiles/Summer15_25nsV6_MC.db)

## Objects ##
### Photons ###
- official cut-based ID and general purpose MVA are taken from [TWiki](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2)
- all photons are used. boolean flags for: loose/medium/tight
- general purpose MVA value is stored

### Jets ###
- ak4PFJetsCHS
- all jets are used
- boolean flag for: loose
- boolean flags for whether a loose electron/photon is found within dR<0.4

### Muons ###
- fulfilling loose id
- tight id boolean flag

### Electrons ###
- fulfilling "veto" id
- boolean flags for loose/medium/tight
- recipes on [TWiki](https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2)

### Generated Particles ###
- genJets collection is stored (= full slimmedGenJets)
- gen[Electrons|Photons]: status=1, pT>30
- Decay products (daughters) of W bosons
