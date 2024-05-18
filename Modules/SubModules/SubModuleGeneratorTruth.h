#ifndef _SubModuleGeneratorTruth_h_
#define _SubModuleGeneratorTruth_h_

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

#include "ubana/CCKaonProduction/Headers/FV.h"
#include "ubana/CCKaonProduction/Headers/ParticleTypes.h"
#include "ubana/CCKaonProduction/Objects/SimParticle.h"
#include "ubana/CCKaonProduction/Objects/Helpers.h"

namespace cckaon {

struct GeneratorTruth {

   double Weight = 1.0;
   int NMCTruths = 0;
   int NMCTruthsInTPC = 0;
   std::vector<std::string> Mode;
   std::vector<SimParticle> Neutrino;
   std::vector<double> TruePrimaryVertex_X;
   std::vector<double> TruePrimaryVertex_Y;
   std::vector<double> TruePrimaryVertex_Z;
   std::vector<std::string> CCNC;
   bool EventHasFinalStateProton=false;
   bool EventHasFinalStatePion=false;
   bool EventHasHyperon=false; 
   bool EventHasKaon=false; 
   bool EventHasKaonP=false; 
   bool EventHasKaonM=false; 
   bool EventHasKaon0=false; 

};

class SubModuleGeneratorTruth {

public:

   SubModuleGeneratorTruth(art::Event const& e,fhicl::ParameterSet pset,bool particlegunmode=false);
   GeneratorTruth GetGeneratorTruth();

private:

   art::Handle<std::vector<simb::MCTruth>> Handle_MCTruth;
   std::vector<art::Ptr<simb::MCTruth>> Vect_MCTruth;
   GeneratorTruth theTruth;

   std::vector<int> KaonPDGs = {321,311,310,130};
   std::vector<int> HyperonPDGs = {3122,3212,3112,3222};
   const bool ParticleGunMode;
};

}

#endif
