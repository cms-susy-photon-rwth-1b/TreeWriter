#ifndef PARTICLECLASSES_H
#define PARTICLECLASSES_H

#include <TLorentzVector.h>
#include <TVector3.h>


#include <math.h>
#include <regex>
#include <time.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TEfficiency.h"
#include "TRandom2.h"
#include "TTree.h"

#include "TLorentzVector.h"

#include "TreeParticles.hpp"
//#include "UserFunctions.h"
//#include "Weighter.h"
//#include "CutFlow.h"
//#include "Resolution.h"

//#include "MT2Functor.h"

#include <iostream>

//#include "config.h"

#include <chrono> //for sleeping
#include <thread> // --do--
#include <cstdlib>//for random increments 
#include <ctime>// --do--

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>


using namespace std;


struct selPhoton : public tree::Photon{
   public:
   void setAll(const tree::Photon& g){
      p=g.p;
      //sigmaIetaIeta=g.sigmaIetaIeta; // full 5x5
      //sigmaIphiIphi=g.sigmaIphiIphi;
      //hOverE=g.hOverE;
      hasPixelSeed=g.hasPixelSeed;
      passElectronVeto=g.passElectronVeto;
      //r9=g.r9;
      //sigmaPt=g.sigmaPt;
      //hasGainSwitch=g.hasGainSwitch;
//
      //cIso=g.cIso;
      //nIso=g.nIso;
      //pIso=g.pIso;
      //cIsoWorst=g.cIsoWorst;
//
      //isTrue=g.isTrue;
      //isTrueAlternative=g.isTrueAlternative;
      //pMultifit=g.pMultifit;
      //pUncorrected=g.pUncorrected;

      isLoose=g.isLoose;
      //isMedium=g.isMedium;
      //isTight=g.isTight;

      //isMediumMVA=g.isMediumMVA;
      //mvaValue=g.mvaValue;
      //mvaCategory=g.mvaCategory;
   }
   TLorentzVector vec;
   float deltaR1;
   float deltaR2;
   bool matchedToPhoton;
   bool matchedToJet;
   bool matchedToElectron;
   
   float scaleFactor;
   float scaleFactorErr;
};

struct selElectron : public tree::Electron{
   public:
   void setAll(const tree::Electron& e){
     
      p=e.p;
      
      charge=e.charge; // +/- 1
      //rIso=e.rIso;
      passImpactParameter=e.passImpactParameter;
      //d0=e.d0;
      //dZ=e.dZ;
      //SIP3D=e.SIP3D;
      miniIso=e.miniIso;
      
      //isVetoID=e.isVetoID;
      isLoose=e.isLoose;
      isMedium=e.isMedium;
      //isTight=e.isTight;
      //isMediumMVA=e.isMediumMVA;
      //isTightMVA=e.isTightMVA;
      //isTightMVASlope=e.isTightMVASlope;
      //mvaValue=e.mvaValue;
      //mvaCategory=e.mvaCategory;
      //seedCrystalE=e.seedCrystalE;
      isPassConvVeto=e.isPassConvVeto;
      //pUncorrected=e.pUncorrected;
   }
   TLorentzVector vec;
   bool matched=false;
   float deltaR1;
   float deltaR2;
   
   float scaleFactor;
   float scaleFactorErr;
   
   int chargeInt;
   
};
struct selMuon : public tree::Muon{
   public:
   void setAll(const tree::Muon& m){
     
      p=m.p;
      
      charge=m.charge; // +/- 1
      //rIso=m.rIso;
      passImpactParameter=m.passImpactParameter;
      //d0=m.d0;
      //dZ=m.dZ;
      //SIP3D=m.SIP3D;
      miniIso=m.miniIso;
      
      //isTight=m.isTight;
      isMedium=m.isMedium;
      //isMediumRun=m.isMediumRun;
      //nTrkLayers=m.nTrkLayers;
   }
   TLorentzVector vec;
   bool matched=false;
   float deltaR1;
   float deltaR2;
   
  float scaleFactor;
   float scaleFactorErr;
   
   int chargeInt;
   
};
struct selLepton : public tree::Lepton{
   public:
   void setAll(const selMuon& m){
     
      p=m.p;
      
      charge=m.charge; // +/- 1
      //rIso=m.rIso;
      passImpactParameter=m.passImpactParameter;
      //d0=m.d0;
      //dZ=m.dZ;
      //SIP3D=m.SIP3D;
      miniIso=m.miniIso;
      
      vec=m.vec;
      matched=m.matched;
      deltaR1=m.deltaR1;
      deltaR2=m.deltaR2;
     
      chargeInt=m.chargeInt;
      
   }
   void setAll(const selElectron& e){
     
      p=e.p;
      
      charge=e.charge; // +/- 1
      //rIso=e.rIso;
      passImpactParameter=e.passImpactParameter;
      //d0=e.d0;
      //dZ=e.dZ;
      //SIP3D=e.SIP3D;
      miniIso=e.miniIso;
      
      vec=e.vec;
      matched=e.matched;
      deltaR1=e.deltaR1;
      deltaR2=e.deltaR2;
     
      chargeInt=e.chargeInt;
   }
   TLorentzVector vec;
   bool matched=false;
   float deltaR1;
   float deltaR2;
   
   float scaleFactor;
   float scaleFactorErr;
   
   int chargeInt;
   
};

struct selJet : public tree::Jet{
   public:
   void setAll(const tree::Jet& j){
      p=j.p;
      isLoose = j.isLoose;
      hasPhotonMatch= j.hasPhotonMatch;
      hasElectronMatch=j.hasElectronMatch;
      hasMuonMatch=j.hasMuonMatch;
      bDiscriminator=j.bDiscriminator;
      //uncert=j.uncert;
      //chf=j.chf;
      //nhf=j.nhf;
      //cef=j.cef;
      //nef=j.nef;
      //nch=j.nch;
      //nconstituents=j.nconstituents;
      //ptRes=j.ptRes;
      //phiRes=j.phiRes;
      //sfRes=j.sfRes;
      //sfResUp=j.sfResUp;
      //sfResDn=j.sfResDn;
      //uncorJecFactor=j.uncorJecFactor; // uncorrected jet momentum over corrected jet momentum
   }
   TLorentzVector vec;
   float deltaR1;
   float deltaR2;
   bool bTag;
};


class selEvent {
    public:
    
    float totalWeight=1.;
    
    bool isDiElectron=false;
    bool isDiMuon=false;
    bool isMuonElectron=false;
    bool isElectronMuon=false;
    //TriggerDecisions(sum)
    bool trigDiEle=false;
    bool trigDiMu=false;
    bool trigMuEle=false;
    bool trigHt=false;
    //TriggerObjectMatchingResults
    bool trigDiEleMatch=false;
    bool trigDiMuMatch=false;
    bool trigMuEleMatch=false;

    //additional variables
    //leptons
    //TLorentzVector l1;
    //TLorentzVector l2;
    tree::Particle4Vector l1;
    tree::Particle4Vector l2;
    float pt1=-10000.;  
    float pt2=-10000.;  
    float phi1=-1000.;  
    float phi2=-1000;  
    float eta1=-1000.;  
    float eta2=-1000.;
    float miniIso1=1000.;  
    float miniIso2=1000.;
    //float charge;
    float chargeProduct;
    float deltaRll=0.;  
    float deltaRl1e=0.;  
    float deltaRl2e=0.;  
    float mll=-1000.0;
    
    //photon
    vector<selPhoton> selPhotons;
    vector<selJet> selJets;
    vector<selElectron> selElectrons;
    vector<selMuon> selMuons;
    //MET
    float ETmiss=-10000.;
    //TLorentzVector ETmiss_vec;
    tree::Particle4Vector ETmiss_vec;
    float MT2_val=-10000.;
    
    float calcHt=0.;
    
    bool evtHasGenPhotonVeto=false;
    
    int matchedEleSize=0;
    int matchedMuSize=0;
    int matchedLeptonSize=0;
    int selLeptonSize=0;
    int selMuonSize=0;
    int selElectronSize=0;
    int selPhotonSize=0;
    
    float lepSF_weight=1.;
    float lepSF_weightUp=1.;
    float lepSF_weightDown=1.;
   
    float photonSF_weight = 1.;
    float photonSF_weightUp=1.;
    float photonSF_weightDown=1.;
    
    float topWeight;
    float topWeightUp;
    float topWeightDown;
    
    float isrWeight;
    float isrWeightUp;
    float isrWeightDown;
    
    float ewkWeight;
    float ewkWeightUp;
    float ewkWeightDown;
    
    int nselJets;
    int nselBJets;
    
};
#endif /* PARTICLECLASSES_H */
