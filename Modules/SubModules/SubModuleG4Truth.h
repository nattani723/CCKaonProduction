#ifndef _SubModuleG4Truth_h_
#define _SubModuleG4Truth_h_

#include <string>
#include <vector>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "cetlib_except/exception.h"

#include "ubana/CCKaonProduction/Headers/ParticleTypes.h"
#include "ubana/CCKaonProduction/Headers/FV.h"
#include "ubana/CCKaonProduction/Objects/SimParticle.h"
#include "ubana/CCKaonProduction/Objects/Helpers.h"

namespace cckaon {

  // Used for comparing vertex positions
  const double _EPSILON_ = 0.0001;
  
  struct G4Truth {
    
    // Flags applying to the entire event 
    // Use for sample orthogonality
    bool EventHasKaonPScatter = false;
    bool EventHasProtonScatter = false;
    bool EventHasPionScatter = false;
    bool EventHasHyperon = false;
    bool EventHasKaon = false;
    bool EventHasKaonP = false;
    bool EventHasKaonP_NuMuP = false;
    bool EventHasKaonP_PiPPi0 = false;
    bool EventHasKaonM = false;
    bool EventHasKaon0 = false;
    
    // Flags for each MCTruth
    std::vector<bool> InActiveTPC;//add more fv cuts
    std::vector<bool> IsHyperon;
    std::vector<bool> IsKaon;
    std::vector<bool> IsKaonP;
    std::vector<bool> IsKaonP_NuMuP;
    std::vector<bool> IsKaonP_PiPPi0;
    std::vector<bool> IsKaonP_2PiPPiM;
    std::vector<bool> IsKaonP_ENuE;
    std::vector<bool> IsKaonP_2PiNPiP;
    std::vector<bool> IsKaonP_Others;
    std::vector<bool> IsKaonM;
    std::vector<bool> IsKaon0;
    std::vector<bool> IsAssociatedKaonP;

    double Weight = 1.0;
    
    std::vector<SimParticle> Lepton;
    std::vector<SimParticle> PrimaryHyperon;
    std::vector<SimParticle> PrimaryNucleon;
    std::vector<SimParticle> PrimaryPion;
    std::vector<SimParticle> PrimaryKaon;
    std::vector<SimParticle> PrimaryKaonP;
    std::vector<SimParticle> PrimaryKaonM;
    std::vector<SimParticle> PrimaryKaon0;
    std::vector<SimParticle> PrimaryNucleus;
    std::vector<SimParticle> HyperonDecay;
    std::vector<SimParticle> KaonPDecay;
    std::vector<SimParticle> KaonPDecay_NuMuP;
    std::vector<SimParticle> KaonPDecay_PiPPi0;
    std::vector<SimParticle> KaonMDecay;
    std::vector<SimParticle> Kaon0Decay;
    std::vector<SimParticle> NeutralKaonDecayK0SL;

    //TVector3 DecayVertex;

    std::vector<double> TruePrimaryVertex_X;
    std::vector<double> TruePrimaryVertex_Y;
    std::vector<double> TruePrimaryVertex_Z;
    
    std::vector<double> DecayVertex_X;
    std::vector<double> DecayVertex_Y;
    std::vector<double> DecayVertex_Z;
    
  };

  class SubModuleG4Truth {
    
  public:
    
    //SubModuleG4Truth();
    SubModuleG4Truth(art::Event const& e,std::string genlabel,std::string g4label,bool particlegunmode=false);
    SubModuleG4Truth(art::Event const& e,fhicl::ParameterSet pset,bool particlegunmode=false);

    void GetParticleLists();
    G4Truth GetG4Info();
    
    void GetPrimaryParticles();
    void GetKaonPDecay();
    void GetKaonMDecay();
    void GetKaon0Decay();
    void GetNeutralKaonDecay();
    void GetHyperonDecay();
    bool FindKaonPScatter();
    bool FindProtonScatter();
    bool FindPionScatter();
    int  GetOrigin(int trackid);
    void MCTruthMatch(SimParticle &P);
    void MCTruthMatch(SimParticle &P,int trackid);
    void SetFlags();
    
    void SetScatterThresholds(double kaonscatterthresh, double pionscatterthresh,double muonscatterthresh);
    void SetDecayThresholds(double decaymuonthresh,double decaypionthresh);

  private:
    
    art::Handle<std::vector<simb::MCTruth>> Handle_MCTruth;
    std::vector<art::Ptr<simb::MCTruth>> Vect_MCTruth;
    
    art::Handle<std::vector<simb::MCParticle>> Handle_G4;
    std::vector<art::Ptr<simb::MCParticle>> Vect_G4;
    
    std::vector<int> Primary_IDs;         // IDs of particles produced at primary vertex
    std::vector<int> Hyperon_Daughter_IDs;        // IDs of Lambda,SigmaP,SigmaM decay products

    std::vector<int> KaonP_Daughter_IDs;   // IDs of Kaon+ decay products
    std::vector<int> KaonP_Inelastic_Daughter_IDs;   // IDs of Kaon+ decay products
    std::vector<int> KaonM_Daughter_IDs;   // IDs of Kaon- decay products
    std::vector<int> Kaon0_Daughter_IDs;   // IDs of Kaon0 decay products
    std::vector<int> NeutralKaon_Daughter_IDs;
    
    std::vector<TVector3> PrimaryVertices;
    bool PosMatch(TVector3 Pos1,TVector3 Pos2);
    
    std::map<int,art::Ptr<simb::MCParticle>> partByID;
    
    std::vector<int> GetChildIDs(const art::Ptr<simb::MCParticle> &g4p,bool IsNeutron=false);
   
    double KaonScatterThresh = 0.0; 
    double MuonScatterThresh = 0.0; 
    double PionScatterThresh = 0.0; 
    double DecayPionThresh = 0.0;
    double DecayMuonThresh = 0.0;
    
    G4Truth theTruth;
    
    int NMCTruths;
    
    const bool ParticleGunMode;
    
  };
  
}

#endif
