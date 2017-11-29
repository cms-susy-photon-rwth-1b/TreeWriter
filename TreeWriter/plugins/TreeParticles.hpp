#ifndef TREEPARTICLES_H
#define TREEPARTICLES_H

#include <TLorentzVector.h>
#include <TVector3.h>

enum PhotonMatchType {UNMATCHED = 0,
                      MATCHED_FROM_GUDSCB,
                      MATCHED_FROM_PI0,
                      MATCHED_FROM_OTHER_SOURCES};

enum PromptStatusType {
  DIRECTPROMPT, FRAGMENTPROMPT, LEPTONPROMPT, NOPROMPT
};

namespace tree
{
   struct Particle
   {
      TVector3 p;
      bool isEB() { return fabs(p.Eta())<1.4442; }
      bool isEE() {
         auto aEta = fabs(p.Eta());
         return 1.566<aEta && aEta < 2.5;
      }
   };

   struct GenParticle: public Particle
   {
      Int_t pdgId=0;
      bool isPrompt;
      bool fromHardProcess;
      UChar_t promptStatus;
   };

   struct IntermediateGenParticle: public GenParticle
   {
      std::vector<GenParticle> daughters;
   };

   struct Photon : public Particle
   {
      Float_t sigmaIetaIeta; // full 5x5
      Float_t sigmaIphiIphi;
      Float_t hOverE;
      Bool_t hasPixelSeed;
      Bool_t passElectronVeto;
      Float_t r9;
      Float_t sigmaPt;
      //Bool_t hasGainSwitch;

      //Float_t cIso;
      //Float_t nIso;
      //Float_t pIso;
      //Float_t cIsoWorst;

      //Int_t isTrue;
      //Int_t isTrueAlternative;
      //TVector3 pMultifit;
      //TVector3 pUncorrected;

      // IDs
      Bool_t  isLoose;
      Bool_t  isMedium;
      Bool_t  isTight;
      //Float_t mva;
      //int isMediumMVA;
      bool isMediumMVA;
      float mvaValue;
      //int mvaCategory;
   };

   struct Jet : public Particle
   {
      bool isLoose;
      bool hasPhotonMatch;
      bool hasElectronMatch;
      bool hasMuonMatch;
      float bDiscriminator;
      float uncert;
      float chf;
      //float nhf;
      //float cef;
      //float nef;
      //int nch;
      //int nconstituents;
      //float ptRes;
      //float phiRes;
      //float sfRes;
      //float sfResUp;
      //float sfResDn;
      //float uncorJecFactor; // uncorrected jet momentum over corrected jet momentum
   };

   struct Lepton: public Particle
   {
      Char_t charge; // +/- 1
      float rIso;
      bool passImpactParameter;
      //float d0;
      //float dZ;
      //float SIP3D;
      float miniIso;
   };

   struct Muon: public Lepton
   {
      bool isTight;
      bool isMedium;
      //bool isMediumRun;
   };

   struct Electron: public Lepton
   {
      bool isVetoID;
      bool isLoose;
      bool isMedium;
      bool isTight;
      //int isMediumMVA;
      //int isTightMVA;
      //int isTightMVASlope;
      bool isMediumMVA;
      bool isTightMVA;
      bool isTightMVASlope;
      //float mvaValue;
      bool mvaValue;
      //int mvaCategory;
      Float_t seedCrystalE;
      bool isPassConvVeto;
   };

   struct MET : public Particle
   {
      Float_t  uncertainty;
      Float_t  sig; // MET significance
   };

   inline bool PtGreater(const tree::Particle p1, const tree::Particle p2) {
      return p1.p.Pt() > p2.p.Pt();
   }

} // end namespace definition
#endif /* TREEPARTICLES_H */
