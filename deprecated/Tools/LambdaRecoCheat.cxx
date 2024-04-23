#ifndef _LambdaRecoCheat_cxx_
#define _LambdaRecoCheat_cxx_

#include "ubana/HyperonProduction/Tools/LambdaRecoCheat.h" 

using namespace hyperon;

void LambdaRecoCheat::MakeHitCollections(art::Event const& e,
    std::vector<std::vector<art::Ptr<recob::Hit>>>& r_hits,
    std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>>& r_hitspacepointmap,
    std::vector<pandora::CartesianVector>& r_vertex) const {
  
  // Create map between MCParticle IDs and the particles 
  art::Handle<std::vector<simb::MCParticle>> Handle_G4;
  std::vector<art::Ptr<simb::MCParticle>> Vect_G4;

  if(!e.getByLabel(params.get<std::string>("G4ModuleLabel"),Handle_G4)) 
    throw cet::exception("LambdaRecoCheat") << "No Geant4 Data Products Found!" << std::endl;

  art::fill_ptr_vector(Vect_G4,Handle_G4);

  std::map<int,art::Ptr<simb::MCParticle>> partByID;
  for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4)
    partByID.insert(std::make_pair(g4p->TrackId(),g4p));

  // Iterate through the list of particles - find any Lambdas and their decay products
  std::vector<unsigned int> lambda_children;
  bool found_proton=false,found_pion=false;
  for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4){
    if(g4p->Mother() != 0 || g4p->PdgCode() != 3122 || g4p->EndProcess() != "Decay") continue;
    for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){
      
      if(partByID.find(g4p->Daughter(i_d)) == partByID.end()) continue;
     
      art::Ptr<simb::MCParticle> daughter = partByID.at(g4p->Daughter(i_d));
      if(daughter->PdgCode() > 10000) continue;

      lambda_children.push_back(daughter->TrackId());
      if(daughter->PdgCode() == 2212){
        found_proton = true;
      }
      if(daughter->PdgCode() == -211){
        found_pion = true;
      }
    } 
  } 

  if(lambda_children.size() != 2 || !found_proton || !found_pion) return;

  std::cout << "Attempting to get lambda daughters hits" << std::endl;

  // Find all of the hits these particles deposit energy in
  std::vector<simb::MCParticle const*> particleVec;
  std::vector<anab::BackTrackerHitMatchingData const*> matchVec;

  art::Handle<std::vector<recob::Hit>> Handle_Hit;
  std::vector<art::Ptr<recob::Hit>> Vect_Hit;
  if(!e.getByLabel(params.get<std::string>("HitModuleLabel"),Handle_Hit)) 
    throw cet::exception("LambdaRecoCheat") << "No Hit Data Products Found!" << std::endl;
  art::fill_ptr_vector(Vect_Hit,Handle_Hit);

  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> ParticlesPerHit(Handle_Hit,e,params.get<std::string>("HitTruthAssnLabel"));
  art::FindManyP<recob::SpacePoint> Assoc_HitSpacePoint(Vect_Hit,e,"pandora");

  r_hits.resize(2);
  r_hitspacepointmap.resize(2);
  r_vertex.resize(2,pandora::CartesianVector(-1000,-1000,-1000)); 
  std::vector<bool> got_start(2,false);

  for(art::Ptr<recob::Hit> hit : Vect_Hit){

    particleVec.clear();
    matchVec.clear();
    ParticlesPerHit.get(hit.key(),particleVec,matchVec);
    std::unordered_map<int,double>  trkide;
 
    // Find the MCParticle that deposited the most energy in this hit
    unsigned int id = 0;
    double maxe = 0.0;
    double e = 0.0;
    for(size_t i_particle=0;i_particle<particleVec.size();++i_particle){
      e += matchVec[i_particle]->energy;
      if(matchVec[i_particle]->energy > maxe){
        id = particleVec[i_particle]->TrackId(); 
        maxe = matchVec[i_particle]->energy;
      }
    }

    if(id == 0 || maxe < MinEnergyDep || e/hit->Integral() < MinEnergyADCRatio) continue;

    for(size_t i_d=0;i_d<lambda_children.size();i_d++){
      
      if(lambda_children.at(i_d) == id){
        std::vector<art::Ptr<recob::SpacePoint>> spacepoints = Assoc_HitSpacePoint.at(hit.key());
        if(spacepoints.size() != 1) continue;
        r_hits.at(i_d).push_back(hit);
        r_hitspacepointmap.at(i_d)[hit] = spacepoints.at(0);
        if(!got_start.at(i_d)){
          r_vertex.at(i_d) = pandora::CartesianVector(spacepoints.at(0)->XYZ()[0],spacepoints.at(0)->XYZ()[1],spacepoints.at(0)->XYZ()[2]);
          got_start.at(i_d) = true;
        }
      }
    }
  } 

  std::cout << "Hits found: " << r_hits.at(0).size() << " " << r_hits.at(1).size() << std::endl;

}


#endif
