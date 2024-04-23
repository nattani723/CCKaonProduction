////////////////////////////////////////////////////////////////////////
// Class:       Throwaway
// Plugin Type: analyzer (art v3_03_01)
// File:        Throwaway_module.cc
//
// Generated at Mon Jan 20 06:07:14 2020 by Christopher Thorpe using cetskelgen
// from cetlib version v3_08_00.
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				

#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

//root includes
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

//local includes

//objects and helpers
#include "ubana/HyperonProduction/Objects/SimParticle.h"
#include "ubana/HyperonProduction/Objects/RecoParticle.h"
#include "ubana/HyperonProduction/Objects/Helpers.h"

//algorithms
#include "ubana/HyperonProduction/Alg/ConnectednessHelper.h"

//submodules
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleGeneratorTruth.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleG4Truth.h"
#include "ubana/HyperonProduction/Modules/SubModules/SubModuleReco.h"

namespace hyperon {
   class Throwaway;
}


class hyperon::Throwaway : public art::EDAnalyzer {
   public:
      explicit Throwaway(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      Throwaway(Throwaway const&) = delete;
      Throwaway(Throwaway&&) = delete;
      Throwaway& operator=(Throwaway const&) = delete;
      Throwaway& operator=(Throwaway&&) = delete;

      // Required functions.
      void analyze(art::Event const& e) override;

      // Selected optional functions.
      void beginJob() override;
      void endJob() override;

      void FinishEvent();

      //check if event contains a reco'd muon, proton and pion from Lambda decay
      //records their positions in track vector if they exist
      //void StoreTrackTruth();

      void beginSubRun(const art::SubRun& sr);
      void endSubRun(const art::SubRun& sr);

   private:

      // Output trees
      TTree * OutputTree;

      fhicl::ParameterSet f_G4;
      fhicl::ParameterSet f_Reco;

      int t_run,t_subrun,t_event;
      std::vector<std::vector<int>> t_TrueParticle;
      std::vector<std::vector<double>> t_TrueParticlePurity;

};

hyperon::Throwaway::Throwaway(fhicl::ParameterSet const& p)
   : EDAnalyzer{p},
   f_G4(p.get<fhicl::ParameterSet>("Geant4")),
   f_Reco(p.get<fhicl::ParameterSet>("Reco"))
{
}

void hyperon::Throwaway::analyze(art::Event const& e)
{

  t_TrueParticle.clear();
  t_TrueParticlePurity.clear();

  t_run = e.run();
  t_subrun = e.subRun();
  t_event = e.event();

  //std::cout << "Running throwaway module"  << std::endl;

  SubModuleG4Truth* G4_SM = new SubModuleG4Truth(e,f_G4,false);
  G4Truth G4T = G4_SM->GetG4Info();

  bool has_sigma = false;
  for(SimParticle hyperon : G4T.Hyperon){
    if((hyperon.PDG == 3112 || hyperon.PDG == 3222) && hyperon.ModMomentum > 0.4 && inActiveTPC(TVector3(hyperon.StartX,hyperon.StartY,hyperon.StartZ))){
      has_sigma = true;
    }
  }

  if(!has_sigma) return;

  art::Handle<std::vector<recob::PFParticle>> Handle_PFParticle;
  std::vector<art::Ptr<recob::PFParticle>> Vect_PFParticle;
  if(!e.getByLabel(f_Reco.get<std::string>("PFParticleModuleLabel"),Handle_PFParticle)) 
    throw cet::exception("LLRPIDTrainer") << "No PFParticle Data Products Found!" << std::endl;
  art::fill_ptr_vector(Vect_PFParticle,Handle_PFParticle);

  art::Handle<std::vector<recob::Track>> Handle_Tracks;
  std::vector<art::Ptr<recob::Track>> Vect_Tracks;
  if(!e.getByLabel(f_Reco.get<std::string>("TrackModuleLabel"),Handle_Tracks))  
    throw cet::exception("LLRPIDTrainer") << "No Track data product!" << std::endl;
  art::fill_ptr_vector(Vect_Tracks,Handle_Tracks);

  art::Handle<std::vector<anab::Calorimetry>> Handle_Calorimetry;
  std::vector<art::Ptr<anab::Calorimetry>> Vect_Calorimetry;
  if(!e.getByLabel(f_Reco.get<std::string>("CaloModuleLabel"),Handle_Calorimetry))  
    throw cet::exception("LLRPIDTrainer") << "No Calorimetry data product!" << std::endl;
  art::fill_ptr_vector(Vect_Calorimetry,Handle_Calorimetry);

  art::Handle<std::vector<recob::Hit>> Handle_Hit;
  std::vector<art::Ptr<recob::Hit>> Vect_Hit;
  if(!e.getByLabel(f_Reco.get<std::string>("HitModuleLabel"),Handle_Hit)) 
    throw cet::exception("SubModuleReco") << "No Hit Data Products Found!" << std::endl;
  art::fill_ptr_vector(Vect_Hit,Handle_Hit);

  art::FindManyP<recob::Track>* Assoc_PFParticleTrack = new art::FindManyP<recob::Track>(Vect_PFParticle,e,f_Reco.get<std::string>("TrackModuleLabel"));    
  art::FindManyP<recob::Hit>* Assoc_TrackHit = new art::FindManyP<recob::Hit>(Vect_Tracks,e,f_Reco.get<std::string>("TrackHitAssnLabel"));
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* ParticlesPerHit = new art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>(Handle_Hit,e,f_Reco.get<std::string>("HitTruthAssnLabel"));

  // Check if there is reco'd neutrino 
  size_t neutrinoID = 99999;
  for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){
    if(pfp->IsPrimary() && (abs(pfp->PdgCode()) == 14 || abs(pfp->PdgCode()) == 12)){
      neutrinoID = pfp->Self();
    }
  }

  if(neutrinoID == 99999) return;

  for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){
    if(pfp->Parent() != neutrinoID) continue; 

    // Grab the associated track
    std::vector<art::Ptr<recob::Track>> pfpTracks = Assoc_PFParticleTrack->at(pfp.key());
    if(pfpTracks.size() != 1) continue;
    art::Ptr<recob::Track> trk = pfpTracks.at(0);

    std::vector<art::Ptr<recob::Hit>> hits = Assoc_TrackHit->at(trk.key());

    std::unordered_map<int,double>  trkide;
    int maxhits=-1;
    int maxhits2=-1;
    int maxhits3=-1;

    simb::MCParticle const* matchedParticle = NULL;
    simb::MCParticle const* matchedParticle2 = NULL;
    simb::MCParticle const* matchedParticle3 = NULL;

    std::vector<simb::MCParticle const*> particleVec;
    std::vector<anab::BackTrackerHitMatchingData const*> matchVec;

    for(size_t i_hit=0;i_hit<hits.size();++i_hit){

      particleVec.clear();
      matchVec.clear();
      ParticlesPerHit->get(hits[i_hit].key(),particleVec,matchVec);

      for(size_t i_particle=0;i_particle<particleVec.size();++i_particle){

        trkide[particleVec[i_particle]->TrackId()]++; 

        if(trkide[particleVec[i_particle]->TrackId()] > maxhits){
          maxhits = trkide[particleVec[i_particle]->TrackId()];
          matchedParticle = particleVec[i_particle];
        }
        else if(trkide[particleVec[i_particle]->TrackId()] > maxhits2){
          maxhits2 = trkide[particleVec[i_particle]->TrackId()];
          matchedParticle2 = particleVec[i_particle];
        }
        else if(trkide[particleVec[i_particle]->TrackId()] > maxhits3){
          maxhits3 = trkide[particleVec[i_particle]->TrackId()];
          matchedParticle3 = particleVec[i_particle];
        }

      }
    }

    int pdg = 0;
    int pdg2 = 0;
    int pdg3 = 0;

    if(matchedParticle != NULL) pdg = matchedParticle->PdgCode(); 
    if(matchedParticle2 != NULL) pdg2 = matchedParticle2->PdgCode(); 
    if(matchedParticle3 != NULL) pdg3 = matchedParticle3->PdgCode(); 

    t_TrueParticle.push_back({pdg,pdg2,pdg3});
    t_TrueParticlePurity.push_back({(double)maxhits/hits.size(),(double)maxhits2/hits.size(),(double)maxhits3/hits.size()});

  }

   OutputTree->Fill();

}

void hyperon::Throwaway::beginJob(){

   art::ServiceHandle<art::TFileService> tfs;
   OutputTree=tfs->make<TTree>("OutputTree","Truth Info Tree");
   OutputTree->Branch("run",&t_run);
   OutputTree->Branch("subrun",&t_subrun);
   OutputTree->Branch("event",&t_event);
   OutputTree->Branch("TrueParticle",&t_TrueParticle);
   OutputTree->Branch("TrueParticlePurity",&t_TrueParticlePurity);

}

void hyperon::Throwaway::endJob()
{
}

void hyperon::Throwaway::beginSubRun(const art::SubRun& sr)
{
}

void hyperon::Throwaway::endSubRun(const art::SubRun& sr){}

DEFINE_ART_MODULE(hyperon::Throwaway)
