// -*- C++ -*-
//
// Package:    TreeWriter/TreeWriter
// Class:      TreeWriter
//

//
// Original Author:  Johannes Lange (adapted parts from Ilya Kravchenko)
//

#include "TreeWriter.hpp"

#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// compute HT using RECO objects to "reproduce" the trigger requirements
static double computeHT(const std::vector<tree::Jet>& jets)
{
   double HT=0;
   double pt=0;
   for (const tree::Jet& jet: jets){
      pt=jet.p.Pt();
      if (fabs(jet.p.Eta())<3.0 && pt>40) HT+=pt;
   }
   return HT;
}

// taken from https://github.com/manuelfs/babymaker/blob/0136340602ee28caab14e3f6b064d1db81544a0a/bmaker/plugins/bmaker_full.cc#L1268-L1295
// "recipe" https://indico.cern.ch/event/557678/contributions/2247944/attachments/1311994/1963568/16-07-19_ana_manuelf_isr.pdf
int n_isr_jets(edm::Handle<edm::View<reco::GenParticle>> const &genParticles,
                            std::vector<tree::Jet> const &jets)
{
   int nisr(0);
   bool matched;
   int momid;
   TVector3 pGen;
   for (tree::Jet const&jet: jets){
      if (jet.hasMuonMatch || jet.hasElectronMatch || jet.hasPhotonMatch) continue;
      matched=false;
      for (size_t imc(0); imc < genParticles->size(); imc++) {
         if (matched) break;
         const reco::GenParticle &mc = (*genParticles)[imc];
         if (mc.status()!=23 || abs(mc.pdgId())>5) continue;
         momid = abs(mc.mother()->pdgId());
         if(!(momid==6 || momid==23 || momid==24 || momid==25 || momid>1e6)) continue;
         //check against daughter in case of hard initial splitting
         for (size_t idau(0); idau < mc.numberOfDaughters(); idau++) {
            pGen.SetXYZ(mc.daughter(idau)->px(),mc.daughter(idau)->py(),mc.daughter(idau)->pz());
            if(jet.p.DeltaR(pGen)<0.3){
               matched = true;
               break;
            }
         }
      } // Loop over MC particles
      if(!matched) {
         nisr++;
      }
   } // Loop over jets
   return nisr;
}

template <typename T> int sign(T val) {
   return (T(0) < val) - (val < T(0));
}


//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
TreeWriter::TreeWriter(const edm::ParameterSet& iConfig)
   : dHT_cut_(iConfig.getUntrackedParameter<double>("HT_cut"))
   , dPhoton_pT_cut_(iConfig.getUntrackedParameter<double>("photon_pT_cut"))
   , dJet_pT_cut_(iConfig.getUntrackedParameter<double>("jet_pT_cut"))
   , isolatedPhotons_(iConfig.getUntrackedParameter<bool>("isolatedPhotons"))
   , minNumberPhotons_cut_(iConfig.getUntrackedParameter<unsigned>("minNumberPhotons_cut"))
   , minNumberElectrons_cut_(iConfig.getUntrackedParameter<unsigned>("minNumberElectrons_cut"))
   , newLumiBlock_(true)
   , vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
   , photonCollectionToken_  (consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photons")))
   , jetCollectionToken_     (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets")))
   , genJetCollectionToken_  (consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets")))
   , muonCollectionToken_    (consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
   , electronCollectionToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons")))
   , metCollectionToken_     (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
   , rhoToken_               (consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))
   , prunedGenToken_         (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticles")))
   , pileUpSummaryToken_     (consumes<PileupSummaryInfoCollection>(iConfig.getParameter<edm::InputTag>("pileUpSummary")))
   , LHEEventToken_          (consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProduct")))
   , METSignificance_        (consumes<double> (iConfig.getParameter<edm::InputTag>("metSig")))
   , packedCandidateToken_   (consumes<std::vector<pat::PackedCandidate>> (iConfig.getParameter<edm::InputTag>("packedCandidates")))
   // electron id
   , electronVetoIdMapToken_  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"   )))
   , electronLooseIdMapToken_ (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"  )))
   , electronMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap" )))
   , electronTightIdMapToken_ (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"  )))
   // photon id
   , photonLooseIdMapToken_  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonLooseIdMap"  )))
   , photonMediumIdMapToken_ (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonMediumIdMap" )))
   , photonTightIdMapToken_  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonTightIdMap"  )))
   , photonMvaValuesMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("photonMvaValuesMap")))
   , phoLooseIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("photonLooseIdMap" )))
   // met filters to apply
   , metFilterNames_(iConfig.getUntrackedParameter<std::vector<std::string>>("metFilterNames"))
   , phoWorstChargedIsolationToken_(consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("phoWorstChargedIsolation")))
   , pileupHistogramName_(iConfig.getUntrackedParameter<std::string>("pileupHistogramName"))
   , hardPUveto_(iConfig.getUntrackedParameter<bool>("hardPUveto"))
   , jetIdSelector(iConfig.getParameter<edm::ParameterSet>("pfJetIDSelector"))
   , triggerNames_(iConfig.getParameter<std::vector<std::string>>("triggerNames"))
   , triggerPrescales_(iConfig.getParameter<std::vector<std::string>>("triggerPrescales"))
   , storeTriggerObjects_(iConfig.getUntrackedParameter<bool>("storeTriggerObjects"))
{
   // declare consumptions that are used "byLabel" in analyze()
   mayConsume<GenLumiInfoHeader,edm::InLumi> (edm::InputTag("generator"));
   consumes<GenEventInfoProduct>(edm::InputTag("generator"));
   consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
   consumes<edm::TriggerResults>(edm::InputTag("TriggerResults",""));
   consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"));
   consumes<std::vector<pat::TriggerObjectStandAlone>>(edm::InputTag("selectedPatTrigger"));

   eventTree_ = fs_->make<TTree> ("eventTree", "event data");

   eventTree_->Branch("photons"  , &vPhotons_);
   eventTree_->Branch("jets"     , &vJets_);
   eventTree_->Branch("genJets"  , &vGenJets_);
   eventTree_->Branch("electrons", &vElectrons_);
   eventTree_->Branch("muons"    , &vMuons_);
   eventTree_->Branch("met"      , &met_);
   eventTree_->Branch("met_raw"  , &met_raw_);
   eventTree_->Branch("met_gen"  , &met_gen_);
   eventTree_->Branch("met_JESu" , &met_JESu_);
   eventTree_->Branch("met_JESd" , &met_JESd_);
   eventTree_->Branch("met_JERu" , &met_JERu_);
   eventTree_->Branch("met_JERd" , &met_JERd_);
   eventTree_->Branch("genParticles", &vGenParticles_);
   if (storeTriggerObjects_) eventTree_->Branch("triggerObjects", &vTriggerObjects_);
   eventTree_->Branch("intermediateGenParticles", &vIntermediateGenParticles_);

   eventTree_->Branch("nPV"           , &nPV_           , "nPV/I");
   eventTree_->Branch("true_nPV"      , &true_nPV_      , "true_nPV/I");
   eventTree_->Branch("nGoodVertices" , &nGoodVertices_ , "nGoodVertices/I");
   eventTree_->Branch("nTracksPV"     , &nTracksPV_     , "nTracksPV/I");
   eventTree_->Branch("rho"           , &rho_           , "rho/F");

   eventTree_->Branch("pu_weight"     , &pu_weight_     , "pu_weight/F");
   eventTree_->Branch("mc_weight"     , &mc_weight_     , "mc_weight/B");
   eventTree_->Branch("pdf_weights"   , &vPdf_weights_);

   eventTree_->Branch("genHt" , &genHt_ , "genHt/F");
   eventTree_->Branch("nISR"  , &nISR_  , "nISR/I");

   eventTree_->Branch("evtNo", &evtNo_, "evtNo/l");
   eventTree_->Branch("runNo", &runNo_, "runNo/i");
   eventTree_->Branch("lumNo", &lumNo_, "lumNo/i");
   eventTree_->Branch("modelName", &modelName_);

   // Fill trigger maps
   for( const auto& n : triggerNames_ ){
      triggerIndex_[n] = -10; // not set and not found
      triggerDecision_[n] = false;
      eventTree_->Branch( n.c_str(), &triggerDecision_[n], (n+"/O").c_str() );
   }
   // create branches for prescales
   for( std::string const& n : triggerPrescales_ ){
      std::string const name=n+"_pre";
      eventTree_->Branch( name.c_str(), &triggerPrescale_[n], (name+"/I").c_str() );
   }

   // get pileup histogram(s)
   std::string cmssw_base_src = getenv("CMSSW_BASE");
   cmssw_base_src += "/src/";
   TFile puFile(TString(cmssw_base_src+"/TreeWriter/PUreweighting/data/puWeights.root"));
   if (puFile.IsZombie() ){
      edm::LogError("File not found") << "create puWeights.root! (see README)";
      std::exit(84);
   } else {
      TH1F* hPU_ptr=(TH1F*)puFile.Get( pileupHistogramName_.c_str() );
      if (hPU_ptr){
         hPU_=*hPU_ptr;
      } else {
         edm::LogError("Pileup histogram not found") << "recreate puWeights.root! (see README)";
         std::exit(84);
      }
   }
   puFile.Close();
}

TH1F* TreeWriter::createCutFlowHist(std::string modelName)
{
   std::string const name("hCutFlow"+modelName);
   std::vector<TString> vCutBinNames{{"initial_unweighted","initial_mc_weighted","initial","METfilters","nGoodVertices", "photons","HT","final"}};
   TH1F* h = fs_->make<TH1F>(name.c_str(),name.c_str(),vCutBinNames.size(),0,vCutBinNames.size());
   for (uint i=0;i<vCutBinNames.size();i++) h->GetXaxis()->SetBinLabel(i+1,vCutBinNames.at(i));
   return h;
}


TreeWriter::~TreeWriter(){}

//
// member functions
//

// ------------ method called for each event  ------------
void
TreeWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   Bool_t  isRealData; // data or MC
   isRealData=iEvent.isRealData();

   hCutFlow_->Fill("initial_unweighted",1);

   // PileUp weights
   if (!isRealData){
      edm::Handle<PileupSummaryInfoCollection>  PupInfo;
      iEvent.getByToken(pileUpSummaryToken_, PupInfo);
      float Tnpv = -1;
      for( auto const& PVI: *PupInfo ) {
         int BX = PVI.getBunchCrossing();
         if(BX == 0) {
            Tnpv = PVI.getTrueNumInteractions();
            continue;
         }
      }
      true_nPV_=Tnpv;
      pu_weight_=hPU_.GetBinContent(hPU_.FindBin(Tnpv));
   }else{ // real data
      true_nPV_=-1;
      pu_weight_=1.;
   }

   // generator weights
   mc_weight_=1; // 1 for data
   if (!isRealData){
      edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
      iEvent.getByLabel("generator", GenEventInfoHandle);
      mc_weight_= sign(GenEventInfoHandle->weight());

      unsigned iMax=110; // these are 9 scale variations and 100 variation of the first pdf set
      if (iMax+1>GenEventInfoHandle->weights().size()) iMax=GenEventInfoHandle->weights().size()-1;
      vPdf_weights_=std::vector<float>(iMax,1.0);
      for (unsigned i=0; i<iMax; i++) {
         // 0 and 1 are the same for 80X scans
         // https://hypernews.cern.ch/HyperNews/CMS/get/susy-interpretations/242/1/1.html
         vPdf_weights_[i]=GenEventInfoHandle->weights()[i+1]/GenEventInfoHandle->weights()[1];
      }
   }

   hCutFlow_->Fill("initial_mc_weighted",mc_weight_);
   hCutFlow_->Fill("initial",mc_weight_*pu_weight_);

   edm::Handle<edm::TriggerResults> triggerBits;
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
   edm::InputTag triggerTag("TriggerResults","","HLT");
   edm::InputTag triggerPrescaleTag("patTrigger");
   iEvent.getByLabel(triggerTag, triggerBits);
   iEvent.getByLabel(triggerPrescaleTag, triggerPrescales);

   // for each lumiBlock, re-read the trigger indices (rather changes for new run)
   if( triggerIndex_.size() && newLumiBlock_ ) {
      newLumiBlock_=false;
      // set all trigger indeces to -1 as "not available"-flag
      for( auto& it : triggerIndex_ )
        it.second = -1;

      // store the indices of the trigger names that we really find
      const edm::TriggerNames &triggerNames = iEvent.triggerNames(*triggerBits);
      for( unsigned i=0; i<triggerNames.size(); i++ ) {
         for( auto& it : triggerIndex_ ) {
            if( triggerNames.triggerName(i).find( it.first ) == 0 ) {
               it.second = i;
            }
         }
      } // end trigger names
   } // found indices

   // set trigger decision
   for( auto& it : triggerIndex_ ) {
      if( it.second != -1 ) {
         triggerDecision_[it.first] = triggerBits->accept( it.second );
      }
   }
   // store prescales
   for( std::string const& n : triggerPrescales_ ){
      int const index = triggerIndex_[n];
      // if the index was not found, store '0': trigger was not run!
      triggerPrescale_[n] = index == -1 ? 0 : triggerPrescales->getPrescaleForIndex(triggerIndex_[n]);
   }

   if (storeTriggerObjects_) {
      edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
      edm::InputTag triggerObjects_("selectedPatTrigger");
      iEvent.getByLabel(triggerObjects_, triggerObjects);

      vTriggerObjects_.clear();
      tree::Particle trObj;
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
         // obj.unpackPathNames(names);
         auto ids = obj.filterIds();
         if (obj.collection() != "hltEgammaCandidates::HLT") continue;
         trObj.p.SetPtEtaPhi(obj.pt(),obj.eta(),obj.phi());
         vTriggerObjects_.push_back(trObj);
      }
   }

   // MET Filters
   edm::Handle<edm::TriggerResults> metFilterBits;
   edm::InputTag metFilterTag("TriggerResults","");
   iEvent.getByLabel(metFilterTag, metFilterBits);
   // go through the filters and check if they were passed
   const edm::TriggerNames &allFilterNames = iEvent.triggerNames(*metFilterBits);
   for (std::string const &name: metFilterNames_){
      const int index=allFilterNames.triggerIndex(name);
      if (!metFilterBits->accept(index)) return; // not passed
   }
   hCutFlow_->Fill("METfilters",mc_weight_*pu_weight_);

   // Get PV
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   if (vertices->empty()) return; // skip the event if no PV found
   nPV_ = vertices->size();

   reco::Vertex firstGoodVertex;
   nGoodVertices_=0;
   for ( const auto& vtx : *vertices ) {
      // from https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_4_14/doc/html/db/d49/GoodVertexFilter_8cc_source.html
      if( vtx.ndof() > 4
         && vtx.position().Rho()<=2.0
         && fabs(vtx.position().Z())<=24.0 )
      {
         nGoodVertices_++;
         // first one?
         if (nGoodVertices_==1) firstGoodVertex = vtx;
      }
   }
   if(!nGoodVertices_) return;
   hCutFlow_->Fill("nGoodVertices",mc_weight_*pu_weight_);

   // Get rho
   edm::Handle< double > rhoH;
   iEvent.getByToken(rhoToken_,rhoH);
   rho_ = *rhoH;

   // get gen particles before photons for the truth match
   edm::Handle<edm::View<reco::GenParticle> > prunedGenParticles;
   if (!isRealData){
      iEvent.getByToken(prunedGenToken_,prunedGenParticles);
   }

   edm::Handle<edm::ValueMap<float> > phoWorstChargedIsolationMap;
   iEvent.getByToken(phoWorstChargedIsolationToken_, phoWorstChargedIsolationMap);

   edm::Handle<edm::ValueMap<bool> > loose_id_dec;
   edm::Handle<edm::ValueMap<bool> > medium_id_dec;
   edm::Handle<edm::ValueMap<bool> > tight_id_dec;
   edm::Handle<edm::ValueMap<float>> mva_value;
   edm::Handle<edm::ValueMap<vid::CutFlowResult> > loose_id_cutflow;
   iEvent.getByToken(photonLooseIdMapToken_  ,loose_id_dec);
   iEvent.getByToken(photonMediumIdMapToken_ ,medium_id_dec);
   iEvent.getByToken(photonTightIdMapToken_  ,tight_id_dec);
   iEvent.getByToken(photonMvaValuesMapToken_,mva_value);
   iEvent.getByToken(phoLooseIdFullInfoMapToken_,loose_id_cutflow);

   // photon collection
   edm::Handle<edm::View<pat::Photon> > photonColl;
   iEvent.getByToken(photonCollectionToken_, photonColl);

   vPhotons_.clear();
   tree::Photon trPho;
   for(edm::View<pat::Photon>::const_iterator pho = photonColl->begin(); pho != photonColl->end(); pho++){
      // Kinematics
      if( pho->pt() < 15 )
         continue;

      trPho.p.SetPtEtaPhi(pho->pt(),pho->superCluster()->eta(),pho->superCluster()->phi());

      const edm::Ptr<pat::Photon> phoPtr( photonColl, pho - photonColl->begin() );

      trPho.sigmaIetaIeta=pho->full5x5_sigmaIetaIeta(); // from reco::Photon
      trPho.sigmaIphiIphi=pho->full5x5_showerShapeVariables().sigmaIphiIphi;
      trPho.hOverE=pho->hadTowOverEm() ;
      trPho.hasPixelSeed=pho->hasPixelSeed() ;
      trPho.passElectronVeto= pho->passElectronVeto() ;
      trPho.r9  = pho->r9();

      vid::CutFlowResult cutFlow = (*loose_id_cutflow)[phoPtr];
      trPho.cIso = cutFlow.getValueCutUpon(4);
      trPho.nIso = cutFlow.getValueCutUpon(5);
      trPho.pIso = cutFlow.getValueCutUpon(6);
      trPho.cIsoWorst = (*phoWorstChargedIsolationMap)[phoPtr];

      trPho.mvaValue=(*mva_value)[phoPtr];

      // MC match
      if (!isRealData){
         trPho.isTrue=matchToTruth(*pho, prunedGenParticles);
         trPho.isTrueAlternative=matchToTruthAlternative(*pho, prunedGenParticles);
      }else{
         trPho.isTrue=           UNMATCHED;
         trPho.isTrueAlternative=UNMATCHED;
      }

      // check photon working points
      trPho.isLoose = (*loose_id_dec) [phoPtr];
      trPho.isMedium= (*medium_id_dec)[phoPtr];
      trPho.isTight = (*tight_id_dec) [phoPtr];

      // write the photon to collection
      if (isolatedPhotons_ && !trPho.isLoose) continue;
      vPhotons_.push_back(trPho);
   } // photon loop

   sort(vPhotons_.begin(), vPhotons_.end(), tree::PtGreater);
   if (minNumberPhotons_cut_ && (vPhotons_.size()<minNumberPhotons_cut_|| vPhotons_.at(0).p.Pt()<dPhoton_pT_cut_)) return;
   hCutFlow_->Fill("photons",mc_weight_*pu_weight_);

   // Muons
   edm::Handle<pat::MuonCollection> muonColl;
   iEvent.getByToken(muonCollectionToken_, muonColl);

   vMuons_.clear();
   tree::Muon trMuon;
   for (const pat::Muon &mu : *muonColl) {
      if (!mu.isLooseMuon()) continue;
      trMuon.p.SetPtEtaPhi(mu.pt(),mu.eta(),mu.phi());
      trMuon.isTight=mu.isTightMuon(firstGoodVertex);
      auto const& pfIso=mu.pfIsolationR04();
      trMuon.rIso=(pfIso.sumChargedHadronPt + std::max(0., pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5*pfIso.sumPUPt))/mu.pt();
      trMuon.charge=mu.charge();
      vMuons_.push_back(trMuon);
   } // muon loop
   sort(vMuons_.begin(), vMuons_.end(), tree::PtGreater);

   // Electrons
   // Get the electron ID data from the event stream
   edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
   edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
   edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
   edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
   iEvent.getByToken(electronVetoIdMapToken_  ,veto_id_decisions);
   iEvent.getByToken(electronLooseIdMapToken_ ,loose_id_decisions);
   iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);
   iEvent.getByToken(electronTightIdMapToken_ ,tight_id_decisions);

   edm::Handle<edm::View<pat::Electron> > electronColl;
   iEvent.getByToken(electronCollectionToken_, electronColl);

   vElectrons_.clear();
   tree::Electron trEl;
   for(edm::View<pat::Electron>::const_iterator el = electronColl->begin();el != electronColl->end(); el++){
      const edm::Ptr<pat::Electron> elPtr(electronColl, el - electronColl->begin() );
      if (!(*veto_id_decisions)[elPtr]) continue; // take only 'veto' electrons
      trEl.isLoose =(*loose_id_decisions) [elPtr];
      trEl.isMedium=(*medium_id_decisions)[elPtr];
      trEl.isTight =(*tight_id_decisions) [elPtr];
      trEl.p.SetPtEtaPhi(el->pt(),el->superCluster()->eta(),el->superCluster()->phi());
      trEl.charge=el->charge();
      vElectrons_.push_back(trEl);
   }
   sort(vElectrons_.begin(), vElectrons_.end(), tree::PtGreater);
   if (vElectrons_.size()<minNumberElectrons_cut_) return;

   // Jets
   edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
   iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);
   JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
   JetCorrectionUncertainty jecUnc(JetCorPar);

   edm::Handle<pat::JetCollection> jetColl;
   iEvent.getByToken(jetCollectionToken_, jetColl);

   vJets_.clear();
   tree::Jet trJet;
   for (const pat::Jet& jet : *jetColl){
      if (jet.pt()<dJet_pT_cut_) continue;
      trJet.p.SetPtEtaPhi(jet.pt(),jet.eta(),jet.phi());
      trJet.bDiscriminator=jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      trJet.isLoose=jetIdSelector(jet);
      jecUnc.setJetEta(jet.eta());
      jecUnc.setJetPt(jet.pt());
      trJet.uncert = jecUnc.getUncertainty(true);
      trJet.chf = jet.chargedHadronEnergyFraction();
      // object matching
      trJet.hasPhotonMatch=false;
      for (tree::Photon const &ph: vPhotons_){
         if (ph.isLoose && trJet.p.DeltaR(ph.p) < 0.4){
            trJet.hasPhotonMatch=true;
            break;
         }
      }
      trJet.hasElectronMatch=false;
      for (tree::Electron const &el: vElectrons_){
         if (el.isLoose && trJet.p.DeltaR(el.p) < 0.4){
            trJet.hasElectronMatch=true;
            break;
         }
      }
      trJet.hasMuonMatch=false;
      for (tree::Muon const &mu: vMuons_){
         if (trJet.p.DeltaR(mu.p) < 0.4){
            trJet.hasMuonMatch=true;
            break;
         }
      }
      vJets_.push_back(trJet);
   } // jet loop
   sort(vJets_.begin(), vJets_.end(), tree::PtGreater);

   // number of ISR jets
   nISR_=0;
   if (!isRealData) {
      nISR_=n_isr_jets(prunedGenParticles,vJets_);
      // std::cout<<nISR_<<std::endl;
   }

   edm::Handle<reco::GenJetCollection> genJetColl;
   if (!isRealData){
     iEvent.getByToken(genJetCollectionToken_, genJetColl);
     vGenJets_.clear();
     tree::Particle trGJet;
     for (const reco::GenJet& jet: *genJetColl){
        if (jet.pt()<dJet_pT_cut_-5) continue;
        trGJet.p.SetPtEtaPhi(jet.pt(),jet.eta(),jet.phi());
        vGenJets_.push_back(trGJet);
     }
     sort(vGenJets_.begin(), vGenJets_.end(), tree::PtGreater);
   } // gen-jet loop

   if (hardPUveto_){
      for (tree::Jet const &j: vJets_){
         if (j.isLoose){
            if (j.p.Pt() > 300) return;
            break; // only check first loose jet
         }
      }
   }

   double const HT=computeHT(vJets_);
   if (HT<dHT_cut_) return;
   hCutFlow_->Fill("HT",mc_weight_*pu_weight_);

   // MET
   edm::Handle<pat::METCollection> metColl;
   iEvent.getByToken(metCollectionToken_, metColl);

   const pat::MET &met = metColl->front();
   double metPt=met.pt();
   met_.p.SetPtEtaPhi(metPt,met.eta(),met.phi());

   if( !isRealData ) {
      const reco::GenMET *genMet=met.genMET();
      met_gen_.p.SetPtEtaPhi(genMet->pt(),genMet->eta(),genMet->phi());
   }

   // jet resolution shift is set to 0 for 74X
   met_.uncertainty=0;
   // loop over all up-shifts save for last one (=NoShift)
   for (uint iShift=0; iShift<(pat::MET::METUncertaintySize-1); iShift+=2){
      // up and down shifts
      const double u=fabs(met.shiftedPt(pat::MET::METUncertainty(iShift))  -metPt);
      const double d=fabs(met.shiftedPt(pat::MET::METUncertainty(iShift+1))-metPt);
      // average
      const double a=.5*(u+d);
      // add deviations in quadrature
      met_.uncertainty+=a*a;
   }
   met_.uncertainty=TMath::Sqrt(met_.uncertainty);

   pat::MET::LorentzVector metShifted;
   metShifted=met.shiftedP4(pat::MET::NoShift, pat::MET::Raw);
   met_raw_.p.SetPtEtaPhi(metShifted.pt(),metShifted.eta(),metShifted.phi());

   metShifted=met.shiftedP4(pat::MET::JetEnUp);
   met_JESu_.p.SetPtEtaPhi(metShifted.pt(),metShifted.eta(),metShifted.phi());
   metShifted=met.shiftedP4(pat::MET::JetEnDown);
   met_JESd_.p.SetPtEtaPhi(metShifted.pt(),metShifted.eta(),metShifted.phi());

   metShifted=met.shiftedP4(pat::MET::JetResUp);
   met_JERu_.p.SetPtEtaPhi(metShifted.pt(),metShifted.eta(),metShifted.phi());
   metShifted=met.shiftedP4(pat::MET::JetResDown);
   met_JERd_.p.SetPtEtaPhi(metShifted.pt(),metShifted.eta(),metShifted.phi());

   edm::Handle<double> METSignificance;
   iEvent.getByToken(METSignificance_, METSignificance);
   met_.sig=float(*METSignificance);
   met_raw_.sig=met_.sig;
   met_JESu_.sig=met_.sig;
   met_JESd_.sig=met_.sig;
   met_JERu_.sig=met_.sig;
   met_JERd_.sig=met_.sig;

   // generated HT
   // stolen from https://github.com/Aachen-3A/PxlSkimmer/blob/master/Skimming/src/PxlSkimmer_miniAOD.cc#L590
   genHt_ = -1;
   if( !isRealData ) {

      edm::Handle<LHEEventProduct> lheInfoHandle;
      iEvent.getByToken(LHEEventToken_, lheInfoHandle);

      if (lheInfoHandle.isValid()) {
         lhef::HEPEUP lheParticleInfo = lheInfoHandle->hepeup();
         // get the five vector
         // (Px, Py, Pz, E and M in GeV)
         std::vector<lhef::HEPEUP::FiveVector> allParticles = lheParticleInfo.PUP;
         std::vector<int> statusCodes = lheParticleInfo.ISTUP;

         genHt_ = 0;
         for (unsigned int i = 0; i < statusCodes.size(); i++) {
            auto absId = abs(lheParticleInfo.IDUP[i]);
            if (statusCodes[i] == 1 && ( absId < 11 || absId > 16 ) && absId != 22  ) {
               genHt_ += sqrt(pow(allParticles[i][0], 2) + pow(allParticles[i][1], 2));
            }
         } // end paricle loop
      } else { // if no lheEventProduct is found
        genHt_ = 0;
        for (const auto& genP : *prunedGenParticles) {
          auto absId = abs(genP.pdgId());
          if (genP.status() == 23 and (absId<11 || absId > 16 ) && genP.pdgId() != 22 ) {
            genHt_ += genP.pt();
          }
        } // genParticle loop
      }
   }


   // Generated Particles
   // status flags: https://indico.cern.ch/event/459797/contribution/2/attachments/1181555/1710844/mcaod-Nov4-2015.pdf
   vGenParticles_.clear();
   vIntermediateGenParticles_.clear();
   tree::GenParticle trP;
   tree::IntermediateGenParticle trIntermP;
   if (!isRealData){
      // Get generator level info
      // Pruned particles are the one containing "important" stuff
      for (const reco::GenParticle &genP: *prunedGenParticles){
         auto absId = abs(genP.pdgId());
         if (absId==23||absId==24){ // store intermediate bosons
            int iNdaugh=genP.numberOfDaughters();
            if (iNdaugh>1){ // skip "decays" V->V
               trIntermP.pdgId = genP.pdgId();
               trIntermP.isPrompt = genP.statusFlags().isPrompt();
               trIntermP.p.SetPtEtaPhi(genP.pt(),genP.eta(),genP.phi());
               trIntermP.daughters.clear();
               for (int i=0; i<iNdaugh; i++){ // store the decay products
                  reco::Candidate const &daugh=*genP.daughter(i);
                  trP.pdgId = daugh.pdgId();
                  trP.isPrompt = false;
                  trP.p.SetPtEtaPhi(daugh.pt(),daugh.eta(),daugh.phi());
                  trIntermP.daughters.push_back(trP);
               }
               vIntermediateGenParticles_.push_back(trIntermP);
            }
         }

         // save particles
         if (genP.status() != 1) continue; // only final state particles
         if (genP.pt() < 20)     continue;
         if (absId==11 || absId==22 || absId==13       // e+-, photon, muon
             || absId==12 || absId==14 || absId==16) { // neutrino
            trP.pdgId = genP.pdgId();
            trP.isPrompt = genP.statusFlags().isPrompt();
            trP.fromHardProcess = genP.statusFlags().fromHardProcess();
            trP.p.SetPtEtaPhi(genP.pt(),genP.eta(),genP.phi());
            vGenParticles_.push_back(trP);
         }
      }
      sort(vGenParticles_.begin(), vGenParticles_.end(), tree::PtGreater);
   }

   // number of tracks
   edm::Handle<std::vector<pat::PackedCandidate>> packedCandidates;
   iEvent.getByToken(packedCandidateToken_, packedCandidates);
   nTracksPV_ = std::count_if(packedCandidates->begin(),packedCandidates->end(), []( const pat::PackedCandidate& cand ) {
      return cand.pt()>.9 && cand.charge() && cand.pvAssociationQuality() == pat::PackedCandidate::UsedInFitTight;});

   hCutFlow_->Fill("final",mc_weight_*pu_weight_);
   // store event identity
   evtNo_=iEvent.id().event();
   runNo_=iEvent.run();
   lumNo_=iEvent.luminosityBlock();
   // write the event
   eventTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
TreeWriter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TreeWriter::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  TreeWriter::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void
  TreeWriter::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------

void
TreeWriter::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&)
{
   newLumiBlock_=true;

   edm::Handle<GenLumiInfoHeader> gen_header;
   // iLumi.getByToken(genLumiHeaderToken, gen_header);
   iLumi.getByLabel("generator", gen_header);
   modelName_ = "";
   if (gen_header.isValid()) {
      modelName_ = gen_header->configDescription();
      std::cout<<modelName_<<std::endl;  // prints, e.g. T1tttt_1500_100
   }

   // create the cutflow histogram for the model if not there yet
   // (modelName="" for non-signal samples)
   if (!hCutFlowMap_.count(modelName_)) {
      hCutFlowMap_[modelName_] = createCutFlowHist(modelName_);
   }
   // point to the right cut flow histogram
   hCutFlow_ = hCutFlowMap_.at(modelName_);
}

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  TreeWriter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

int TreeWriter::matchToTruth(const pat::Photon &pho,
                                              const edm::Handle<edm::View<reco::GenParticle>>
                                              &genParticles)
{
   //
   // Explicit loop and geometric matching method
   //

   // Find the closest status 1 gen photon to the reco photon
   double dR = 999;
   const reco::Candidate *closestPhoton = 0;
   for( auto const& particle: *genParticles ){
      // Drop everything that is not photon or not status 1
      if( abs(particle.pdgId()) != 22 || particle.status() != 1 )
         continue;
      //
      double dRtmp = ROOT::Math::VectorUtil::DeltaR( pho.p4(), particle.p4() );
      if( dRtmp < dR ){
         dR = dRtmp;
         closestPhoton = &particle;
      }
   }
   // See if the closest photon (if it exists) is close enough.
   // If not, no match found.
   if( !(closestPhoton != 0 && dR < 0.1) ) {
      return UNMATCHED;
   }

   // Find ID of the parent of the found generator level photon match
   int ancestorPID = -999;
   int ancestorStatus = -999;
   findFirstNonPhotonMother(closestPhoton, ancestorPID, ancestorStatus);

   // Allowed parens: quarks pdgId 1-5, or a gluon 21
   std::vector<int> allowedParents { -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -21, 21 };
   if( !(std::find(allowedParents.begin(),
                   allowedParents.end(), ancestorPID)
         != allowedParents.end()) ){
      // So it is not from g, u, d, s, c, b. Check if it is from pi0 or not.
      if( abs(ancestorPID) == 111 )
         return MATCHED_FROM_PI0;
      else
         return MATCHED_FROM_OTHER_SOURCES;
   }
   return MATCHED_FROM_GUDSCB;

}

void TreeWriter::findFirstNonPhotonMother(const reco::Candidate *particle,
                                                           int &ancestorPID, int &ancestorStatus){

   if( particle == 0 ){
      printf("TreeWriter: ERROR! null candidate pointer, this should never happen\n");
      return;
   }

   // Is this the first non-photon parent? If yes, return, otherwise
   // go deeper into recursion
   if( abs(particle->pdgId()) == 22 ){
      findFirstNonPhotonMother(particle->mother(0), ancestorPID, ancestorStatus);
   }else{
      ancestorPID = particle->pdgId();
      ancestorStatus = particle->status();
   }

   return;
}

int TreeWriter::matchToTruthAlternative(const pat::Photon &pho,
                                                         const edm::Handle<edm::View<reco::GenParticle>>
                                                         &genParticles)
{


   //
   // Explicit loop and geometric matching method
   //

   int isMatched = UNMATCHED;

   for( auto const& particle: *genParticles ){
      int pid = particle.pdgId();
      int ancestorPID = -999;
      int ancestorStatus = -999;
      findFirstNonPhotonMother(&particle, ancestorPID, ancestorStatus);
      if( pid ==22 && TMath::Abs( ancestorPID ) <= 22 ){
         double dr = ROOT::Math::VectorUtil::DeltaR( pho.p4(), particle.p4() );
         float dpt = fabs( (pho.pt() - particle.pt() )/particle.pt());
         if (dr < 0.2 && dpt < 0.2){
            isMatched = MATCHED_FROM_GUDSCB;
            if( ancestorPID == 22 ){
               printf("Ancestor of a photon is a photon!\n");
            }
         }
      }
   }

   return isMatched;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeWriter);
