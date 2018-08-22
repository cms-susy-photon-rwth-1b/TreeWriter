#include "TreeParticles.hpp"
#include "particleClasses.hpp"
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class tree::Particle+;
#pragma link C++ class tree::Photon+;
#pragma link C++ class tree::Jet+;
#pragma link C++ class tree::Lepton+;
#pragma link C++ class tree::Muon+;
#pragma link C++ class tree::Electron+;
#pragma link C++ class tree::MET+;
#pragma link C++ class tree::GenParticle+;
#pragma link C++ class tree::IntermediateGenParticle+;
#pragma link C++ class std::vector<tree::Photon>+;
#pragma link C++ class std::vector<tree::Particle>+;
#pragma link C++ class std::vector<tree::Jet>+;
#pragma link C++ class std::vector<tree::Muon>+;
#pragma link C++ class std::vector<tree::Electron>+;
#pragma link C++ class std::vector<tree::GenParticle>+;
#pragma link C++ class std::vector<tree::IntermediateGenParticle>+;

#pragma link C++ class std::vector<selPhoton>+;
#pragma link C++ class std::vector<selJet>+;
#pragma link C++ class std::vector<selElectron>+;
#pragma link C++ class std::vector<selMuon>+;

#pragma link C++ class selPhoton+;
#pragma link C++ class selJet+;
#pragma link C++ class selElectron+;
#pragma link C++ class selMuon+;

#pragma link C++ class tree::Particle4Vector+;
#pragma link C++ class std::vector<tree::Particle4Vector>+;

//#pragma link C++ class TLorentzVector+;


#endif
