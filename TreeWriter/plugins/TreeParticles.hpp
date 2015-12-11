#ifndef TREEPARTICLES_H
#define TREEPARTICLES_H

#include <TLorentzVector.h>
#include <TVector3.h>

enum PhotonMatchType {UNMATCHED = 0,
                      MATCHED_FROM_GUDSCB,
                      MATCHED_FROM_PI0,
                      MATCHED_FROM_OTHER_SOURCES};

namespace tree
{
   struct Particle
   {
      TVector3 p;
   };

   struct GenParticle: public Particle
   {
      Int_t pdgId=0;
   };

   struct Photon : public Particle
   {
      Float_t sigmaIetaIeta; // full 5x5
      Float_t hOverE;
      Int_t hasPixelSeed;
      Int_t passElectronVeto;
      Float_t r9;

      // offline isolations
      Float_t isoChargedHadronsEA;
      Float_t isoNeutralHadronsEA;
      Float_t isoPhotonsEA;
      Float_t isoWorstChargedHadrons;

      // trigger(-like) isolations
      Float_t ecalPFClIso;
      Float_t hcalPFClIso;
      Float_t trackIso; // this is a solid cone Iso, HLT uses hollow cone

      Int_t isTrue;
      Int_t isTrueAlternative;

      // IDs
      Bool_t  isLoose;
      Bool_t  isMedium;
      Bool_t  isTight;
      Float_t mvaValue;
   };

   struct Jet : public Particle
   {
      bool isLoose;
      bool hasPhotonMatch;
      bool hasElectronMatch;
      bool hasMuonMatch;
      float bDiscriminator;
      float pileUpDiscriminator;
      float chf;
      float nhf;
      float cef;
      float nef;
      int nch;
      int nconstituents;
   };

   struct Muon: public Particle
   {
      bool isTight;
   };

   struct Electron: public Particle
   {
      bool isLoose;
      bool isMedium;
      bool isTight;
   };

   struct MET : public Particle
   {
      TVector3 p_raw;
      Float_t  uncertainty;
      Float_t  sig; // MET significance
   };

   inline bool PtGreater(const tree::Particle p1, const tree::Particle p2) {
      return p1.p.Pt() > p2.p.Pt();
   }

} // end namespace definition
#endif /* TREEPARTICLES_H */
