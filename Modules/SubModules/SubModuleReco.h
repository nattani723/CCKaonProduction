#ifndef _SubModuleReco_h_
#define _SubModuleReco_h_

#include <string>
#include <vector>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"	
#include "canvas/Persistency/Common/TriggerResults.h"			
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"
#include "cetlib_except/exception.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "ubana/CCKaonProduction/Headers/ParticleTypes.h"
#include "ubana/CCKaonProduction/Headers/LLR_PID.h"
#include "ubana/CCKaonProduction/Headers/LLRPID_proton_muon_lookup.h"
#include "ubana/CCKaonProduction/Headers/LLR_PID_K.h"
#include "ubana/CCKaonProduction/Headers/LLRPID_kaon_proton_lookup.h"
#include "ubana/CCKaonProduction/Objects/RecoParticle.h"
#include "ubana/CCKaonProduction/Objects/Helpers.h"
#include "ubana/CCKaonProduction/Alg/PIDManager.h"
#include "ubana/CCKaonProduction/Modules/SubModules/SubModuleG4Truth.h"
#include "ubana/CCKaonProduction/Alg/BDTHandle.h"

#include "TVector3.h"

using std::string;

namespace cckaon {

struct RecoData {

   TVector3 RecoPrimaryVertex = TVector3(-1000,-1000,-1000);

   int NPrimaryDaughters; 
   int NPrimaryTrackDaughters;
   int NPrimaryShowerDaughters;
   int NOtherTracks;
   int NOtherRebuiltTracks;
   int NOtherShowers;

   std::vector<RecoParticle> TrackPrimaryDaughters;
   std::vector<RecoParticle> ShowerPrimaryDaughters;
   std::vector<RecoParticle> TrackOthers;
   std::vector<RecoParticle> TrackRebuiltOthers;
   std::vector<RecoParticle> ShowerOthers;
  RecoParticle CCMuTrack;

   std::vector<TVector3> TrackStarts;

   size_t TrueMuonIndex = -1;
   size_t TrueKaonIndex = -1;
   size_t TrueDecayMuonIndex = -1;
   size_t TrueDecayPionIndex = -1;

   bool GoodReco = false;
   bool GoodReco_NuMuP = false;
   bool GoodReco_PiPPi0 = false;
   bool GoodPrimaryReco = false;
   bool GoodRecoAsShower = false;

};

class SubModuleReco {

   public:

      //SubModuleReco();
  SubModuleReco(art::Event const& e,bool isdata,string pfparticlelabel,string tracklabel,string trackrebuiltlabel,
		string showerlabel,string vertexlabel, string pidlabel, string calipidlabel,string calolabel,string hitlabel,
		string hittruthassnlabel,string trackhitassnlabel,string trackrebuilthitassnlabel,
		string rerunpidlabel, string reruncalipidlabel, string reruncalolabel, string showerhitassnlabel,string metadatalabel,string genlabel,
		string g4label,bool dogetpids,bool includecosmics,bool particlegunmode=false, bool withrecoalg=false);

  SubModuleReco(art::Event const& e,bool isdata,fhicl::ParameterSet pset,bool particlegunmode=false, bool withrecoalg=false);

      void PrepareInfo(); 
      TVector3 GetPrimaryVertex();
      void SetIndices(std::vector<bool> IsSignal, std::vector<bool> IsSignal_NuMuP, std::vector<bool> IsSignal_PiPPi0);

      RecoData GetInfo();
      void SetResRangeCutoff(double cutoff){ ResRangeCutoff = cutoff; }
      bool ApplyNuCCInclusiveFilter(art::Event const& e);

     

   private:

      art::Handle<std::vector<recob::PFParticle>> Handle_PFParticle;
      std::vector<art::Ptr<recob::PFParticle>> Vect_PFParticle;

      art::Handle<std::vector<recob::Track>> Handle_Track;
      std::vector<art::Ptr<recob::Track>> Vect_Track;

      art::Handle<std::vector<recob::Track>> Handle_TrackRebuilt;
      std::vector<art::Ptr<recob::Track>> Vect_TrackRebuilt;

      art::Handle<std::vector<recob::Shower>> Handle_Shower;
      std::vector<art::Ptr<recob::Shower>> Vect_Shower;

      art::Handle<std::vector<recob::Hit>> Handle_Hit;
      std::vector<art::Ptr<recob::Hit>> Vect_Hit;

      RecoParticle MakeRecoParticle(const art::Ptr<recob::PFParticle> &pfp);
      RecoParticle MakeRecoParticle(const art::Ptr<recob::Track> &trk);
      RecoParticle MakeRecoParticle(const art::Ptr<recob::Shower> &shw);

      art::FindManyP<anab::T0>* Assoc_PFPMuon;
      art::FindManyP<recob::Vertex>* Assoc_PFParticleVertex;  
      art::FindManyP<recob::Track>* Assoc_PFParticleTrack;
      //art::FindManyP<recob::Track>* Assoc_PFParticleTrackRebuilt;
      art::FindManyP<recob::Shower>* Assoc_PFParticleShower;
      art::FindManyP<larpandoraobj::PFParticleMetadata>* Assoc_PFParticleMetadata;
      art::FindManyP<recob::Hit>* Assoc_TrackHit;
      art::FindManyP<recob::Hit>* Assoc_TrackRebuiltHit;
      art::FindManyP<recob::Hit>* Assoc_ShowerHit;
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* Assoc_MCParticleBacktracker;
      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>* ParticlesPerHit;
      art::FindManyP<anab::Calorimetry>* Assoc_TrackCalo;
      art::FindManyP<anab::ParticleID>* Assoc_TrackPID;
      art::FindManyP<anab::ParticleID>* Assoc_TrackCaliPID;
      art::FindManyP<anab::Calorimetry>* Assoc_TrackRebuiltCalo;
      art::FindManyP<anab::ParticleID>* Assoc_TrackRebuiltPID;
      art::FindManyP<anab::ParticleID>* Assoc_TrackRebuiltCaliPID;

      searchingfornues::LLRPID llr_pid_calculator;
      searchingfornues::ProtonMuonLookUpParameters protonmuon_parameters;

      searchingfornuesk::LLRPIDK llr_pid_calculator_kaon;
      searchingfornuesk::KaonProtonLookUpParameters kaonproton_parameters;

      SubModuleG4Truth* G4T = nullptr;
      PIDManager PIDCalc;      

      RecoData theData;
      size_t neutrinoID = 99999;
      std::map<size_t,int> m_PFPID_TrackIndex;

      void GetPFPMetadata(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P);
      void GetTrackData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P);
      void GetTrackData(const art::Ptr<recob::Track> &trk,RecoParticle &P);
      void GetRebuiltTrackData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P);
      void GetShowerData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P);
      void GetShowerData(const art::Ptr<recob::Shower> &shw,RecoParticle &P);
      void TruthMatch(const art::Ptr<recob::Track> &trk,RecoParticle &P);
      void MergeCheck(const std::vector<art::Ptr<recob::Hit>>& hits, RecoParticle &P);
      void GetPIDs(const art::Ptr<recob::Track> &trk,RecoParticle &P);
      void GetVertexData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P);

      bool IsData;
      bool DoGetPIDs=true;
      bool WithRecoAlgorithm=false;
      double ResRangeCutoff=5; 
      const bool IncludeCosmics;
      const bool ParticleGunMode;
};

}

#endif


