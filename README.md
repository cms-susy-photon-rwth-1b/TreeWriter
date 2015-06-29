**TreeWriter** to build a ROOT tree from MiniAOD. Photon Cut- and MVA-IDs are computed.

## Building and Running ##
Get CMSSW environment 74X

```
cmsrel CMSSW_7_4_5
cd CMSSW_7_4_5/src/
cmsenv
```
Get and build egamma recipes (for MVA)

```
git cms-merge-topic ikrav:egm_id_74X_v2
scram b -j 8
```
Get and build the TreeWriter

```
git clone https://github.com/cms-susy-photon-rwth-1b/TreeWriter.git
cd TreeWriter
scram b
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
crab submit -c crabConfig.py
```

## Configure ##
in the python config, set
- `HT_cut`: minimum HT
- `photon_pT_cut`: minimum leading-photon pT

## Objects ##
### Photons ###
- based on this recipe for MVA ID: [HN](https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1552.html)
- this Cut-ID is included manually: [HN](https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1541.html) (CMSSW72X)
  * for CMSSW74X: [HN](https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1573.html) [TWiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2#Recipe_for_regular_users_for_74X)
- all photons are used. boolean flags for: loose/medium/tight

### Jets ###
- ak4PFJetsCHS
- all jets are used
- boolean flag for: loose

### Muons ###
- fulfilling loose id
- tight id boolean flag

### Electrons ###
- fulfilling "veto" id
- boolean flags for loose/medium/tight
- recipes:
  * CMSSW72X [TWiki](https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_min)
  * CMSSW74X [TWiki](https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_74X)

### Generated Particles ###
- genJets collection is stored (= full slimmedGenJets)
- gen[Electrons|Photons]: status=1, pT>30

# TODO #
Replace the MVA-ID by official recipe when available.
