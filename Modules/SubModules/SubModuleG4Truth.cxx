#ifndef _SubModuleG4Truth_cxx_
#define _SubModuleG4Truth_cxx_

#include "SubModuleG4Truth.h"

using namespace cckaon;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

//SubModuleG4Truth::SubModuleG4Truth(){}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleG4Truth::SubModuleG4Truth(art::Event const& e,std::string genlabel,std::string g4label,bool particlegunmode) :
ParticleGunMode(particlegunmode)
{

   if(!e.getByLabel(genlabel,Handle_MCTruth))  
      throw cet::exception("SubModuleG4Truth") << "No MC Truth data product!" << std::endl;

   if(!e.getByLabel(g4label,Handle_G4)) 
      throw cet::exception("SubModuleG4Truth") << "No Geant4 Data Products Found!" << std::endl;

   art::fill_ptr_vector(Vect_MCTruth,Handle_MCTruth);  
   art::fill_ptr_vector(Vect_G4,Handle_G4);

   // Create map between particle ID and Particles
   for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4)
      partByID.insert(std::make_pair(g4p->TrackId(),g4p));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleG4Truth::SubModuleG4Truth(art::Event const& e,fhicl::ParameterSet pset,bool particlegunmode) :
   SubModuleG4Truth(e,
         pset.get<std::string>("GeneratorModuleLabel","generator"),
         pset.get<std::string>("G4ModuleLabel","largeant"),
         particlegunmode)
{
  //SetNeutronScatterThresholds(pset.get<double>("NeutronScatterProtonThresh",0.15),pset.get<double>("NeutronScatterPionThresh",0.05));
  //SetDecayThresholds(pset.get<double>("DecayProtonThresh",0.0),pset.get<double>("DecayPionThresh",0.0));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetParticleLists(){

   Primary_IDs.clear();
   Hyperon_Daughter_IDs.clear();
   KaonP_Daughter_IDs.clear();
   KaonP_Inelastic_Daughter_IDs.clear();
   KaonM_Daughter_IDs.clear();
   Kaon0_Daughter_IDs.clear();
   NeutralKaon_Daughter_IDs.clear();

   for(const art::Ptr<simb::MCParticle> &g4p : Vect_G4){

      // Skip any particles not produced at PV
      if(g4p->Mother() != 0)  continue;

      Primary_IDs.push_back(g4p->TrackId());

      std::vector<int> IDs = GetChildIDs(g4p);

      if(isHyperon(g4p->PdgCode())){
	if(g4p->EndProcess() == "Decay")
	  Hyperon_Daughter_IDs.insert(Hyperon_Daughter_IDs.begin(),IDs.begin(),IDs.end()); 
      }		

      if(isKaon(g4p->PdgCode())){

	//if(g4p->PdgCode() == 321) std::cout << "g4p->EndProcess(): " << g4p->EndProcess() << std::endl;
        if(g4p->PdgCode() == 321 && (g4p->EndProcess() == "Decay" || g4p->EndProcess() == "FastScintillation") )
	  KaonP_Daughter_IDs.insert(KaonP_Daughter_IDs.begin(),IDs.begin(),IDs.end());

        if(g4p->PdgCode() == 321 && g4p->EndProcess() == "kaon+Inelastic")
	  KaonP_Inelastic_Daughter_IDs.insert(KaonP_Inelastic_Daughter_IDs.begin(),IDs.begin(),IDs.end());

        if(g4p->PdgCode() == -321 && g4p->EndProcess() == "Decay")
	  KaonM_Daughter_IDs.insert(KaonM_Daughter_IDs.begin(),IDs.begin(),IDs.end());

	if(g4p->PdgCode() == 311 && g4p->EndProcess() == "Decay")
	  NeutralKaon_Daughter_IDs.insert(NeutralKaon_Daughter_IDs.begin(),IDs.begin(),IDs.end());

        else if(g4p->EndProcess() == "Decay" || g4p->EndProcess() == "FastScintillation")
          Kaon0_Daughter_IDs.insert(Kaon0_Daughter_IDs.begin(),IDs.begin(),IDs.end());                     
      }
   
   }


   // Get the Neutral Kaon daughters
   // NOTE: If GENIE produces a K0 (pdg 311), it will instantly decay it into a K0S (pdg 310)/K0L (pdg 130), then propagates  
   // those. PrimaryK0SL is to capture the ID of the K0S/K0L, in order to then identify their decay products

   for(size_t i_d=0;i_d<NeutralKaon_Daughter_IDs.size();i_d++){

     if(partByID.find(NeutralKaon_Daughter_IDs[i_d]) == partByID.end()) continue;
     art::Ptr<simb::MCParticle> part = partByID[NeutralKaon_Daughter_IDs.at(i_d)];
     if(part->PdgCode() != 130 && part->PdgCode() != 310)
       throw std::invalid_argument("SubModuleG4Truth::GetParticleLists: Neutral kaon daughter with pdg code " + std::to_string(part->PdgCode()) + " expected 130 or 310");
     if(part->EndProcess() == "Decay"){
       std::vector<int> IDs = GetChildIDs(part);   
       Kaon0_Daughter_IDs.insert(Kaon0_Daughter_IDs.begin(),IDs.begin(),IDs.end());  
     }
   }


   // Set list of Primary vertices for matching to multiple MCTruths
   for(const art::Ptr<simb::MCTruth> &theMCTruth : Vect_MCTruth){         
      if(!ParticleGunMode){
         for(int k_particles=0;k_particles<theMCTruth->NParticles();k_particles++){
            simb::MCParticle Part = theMCTruth->GetParticle(k_particles);
            if(isNeutrino(Part.PdgCode()) && Part.StatusCode() == 0){ 
               PrimaryVertices.push_back(TVector3(Part.Vx(),Part.Vy(),Part.Vz()));
               theTruth.TruePrimaryVertex_X.push_back(Part.Vx());
               theTruth.TruePrimaryVertex_Y.push_back(Part.Vy());
               theTruth.TruePrimaryVertex_Z.push_back(Part.Vz());
            }
         }
      }       
      else if(theMCTruth->NParticles()){ // this is particle gun
         PrimaryVertices.push_back(TVector3(theMCTruth->GetParticle(0).Vx(),theMCTruth->GetParticle(0).Vy(),theMCTruth->GetParticle(0).Vz()));
         theTruth.TruePrimaryVertex_X.push_back(theMCTruth->GetParticle(0).Vx());
         theTruth.TruePrimaryVertex_Y.push_back(theMCTruth->GetParticle(0).Vy());
         theTruth.TruePrimaryVertex_Z.push_back(theMCTruth->GetParticle(0).Vz());
      }
      else if(theMCTruth->NParticles() == 0)
         throw cet::exception("SubModuleG4Truth") << "MCTruth made with particle gun contains no particles!" << std::endl;
   }

   NMCTruths = Vect_MCTruth.size();

   if(/*!ParticleGun &&*/ PrimaryVertices.size() != Vect_MCTruth.size())         
      throw cet::exception("SubModuleG4Truth") << "Vertex/MCTruth vector size mismatch" << std::endl;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

G4Truth SubModuleG4Truth::GetG4Info(){

   GetParticleLists();

   // Clear everything

   theTruth.InActiveTPC.resize(NMCTruths);
   theTruth.IsHyperon.resize(NMCTruths);
   theTruth.IsKaon.resize(NMCTruths);
   theTruth.IsKaonP.resize(NMCTruths);
   theTruth.IsKaonP_NuMuP.resize(NMCTruths);
   theTruth.IsKaonP_PiPPi0.resize(NMCTruths);
   theTruth.IsKaonP_2PiPPiM.resize(NMCTruths);
   theTruth.IsKaonP_ENuE.resize(NMCTruths);
   theTruth.IsKaonP_2PiNPiP.resize(NMCTruths);
   theTruth.IsKaonP_Others.resize(NMCTruths);
   theTruth.IsKaonM.resize(NMCTruths);
   theTruth.IsKaon0.resize(NMCTruths);
   theTruth.IsAssociatedKaonP.resize(NMCTruths);

   theTruth.DecayVertex_X.resize(NMCTruths);
   theTruth.DecayVertex_Y.resize(NMCTruths);
   theTruth.DecayVertex_Z.resize(NMCTruths);

   for(int i_t=0;i_t<NMCTruths;i_t++){
      theTruth.DecayVertex_X[i_t] = -1000;
      theTruth.DecayVertex_Y[i_t] = -1000;
      theTruth.DecayVertex_Z[i_t] = -1000;
   }

   theTruth.Lepton.clear();
   theTruth.PrimaryHyperon.clear();
   theTruth.PrimaryNucleon.clear();
   theTruth.PrimaryPion.clear();
   theTruth.PrimaryKaon.clear();
   theTruth.PrimaryKaonP.clear();
   theTruth.PrimaryKaonM.clear();
   theTruth.PrimaryKaon0.clear();
   theTruth.PrimaryNucleus.clear();
   theTruth.HyperonDecay.clear();
   theTruth.KaonPDecay.clear();
   theTruth.KaonPDecay_NuMuP.clear();
   theTruth.KaonPDecay_PiPPi0.clear();
   theTruth.KaonMDecay.clear();
   theTruth.Kaon0Decay.clear();
   theTruth.NeutralKaonDecayK0SL.clear();

   GetPrimaryParticles();

   if(theTruth.PrimaryHyperon.size()) GetHyperonDecay();
   if(NeutralKaon_Daughter_IDs.size()) GetNeutralKaonDecay();
   if(theTruth.PrimaryKaon.size() || theTruth.PrimaryKaon0.size()) GetKaon0Decay();
   if(theTruth.PrimaryKaonP.size()) GetKaonPDecay();
   if(theTruth.PrimaryKaonM.size()) GetKaonMDecay();

   SetFlags(); 

   return theTruth;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetPrimaryParticles(){

   for(size_t i_p=0;i_p<Primary_IDs.size();i_p++){

      // Geant does not always keep everything it simulates, make sure particle is in list of IDs (crashes otherwise!)
      if(partByID.find(Primary_IDs[i_p]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[Primary_IDs.at(i_p)];

      // Anything with very large pdg code is a nucleus, skip these
      //if(part->PdgCode() > 10000) continue;

      SimParticle P = MakeSimParticle(*part);
      P.Origin = 1;
      MCTruthMatch(P);

      //hyperon produced at primary vertex
      if(isHyperon(part->PdgCode())) theTruth.PrimaryHyperon.push_back(P);
      if(isLepton(part->PdgCode()) || isNeutrino(part->PdgCode())) theTruth.Lepton.push_back(P);
      if(isNucleon(part->PdgCode())) theTruth.PrimaryNucleon.push_back(P);
      if(isPion(part->PdgCode())) theTruth.PrimaryPion.push_back(P);
      if(isKaon(part->PdgCode())) theTruth.PrimaryKaon.push_back(P);
      if(isKaonP(part->PdgCode())) theTruth.PrimaryKaonP.push_back(P);
      if(isKaonM(part->PdgCode())) theTruth.PrimaryKaonM.push_back(P);
      if(isKaon0(part->PdgCode())) theTruth.PrimaryKaon0.push_back(P);
      if(part->PdgCode() > 10000) theTruth.PrimaryNucleus.push_back(P);
   }

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetKaonPDecay(){

   if(!KaonP_Daughter_IDs.size()) return;

   for(size_t i_d=0;i_d<KaonP_Daughter_IDs.size();i_d++){

      // Geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
      if(partByID.find(KaonP_Daughter_IDs[i_d]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[KaonP_Daughter_IDs[i_d]];

      SimParticle P = MakeSimParticle(*part);
      P.Origin = 2;

      // Check which MCTruth this decay belongs to
      for(size_t i_h=0;i_h<theTruth.PrimaryKaonP.size();i_h++){
         SimParticle H = theTruth.PrimaryKaonP.at(i_h);
         if(PosMatch(TVector3(P.StartX,P.StartY,P.StartZ),TVector3(H.EndX,H.EndY,H.EndZ))){
            P.MCTruthIndex = H.MCTruthIndex;            
            theTruth.DecayVertex_X[P.MCTruthIndex] = P.StartX;
            theTruth.DecayVertex_Y[P.MCTruthIndex] = P.StartY;
            theTruth.DecayVertex_Z[P.MCTruthIndex] = P.StartZ;
         }
      } 

      theTruth.KaonPDecay.push_back(P);
      if(part->PdgCode() == -13) {
	theTruth.KaonPDecay_NuMuP.push_back(P);
	P.Origin = 7;
      }
      else if(part->PdgCode() == 211) {
	theTruth.KaonPDecay_PiPPi0.push_back(P);
	P.Origin = 8;
      }
   }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetKaonMDecay(){

   if(!KaonM_Daughter_IDs.size()) return;

   for(size_t i_d=0;i_d<KaonM_Daughter_IDs.size();i_d++){

      // Geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
      if(partByID.find(KaonM_Daughter_IDs[i_d]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[KaonM_Daughter_IDs[i_d]];

      SimParticle P = MakeSimParticle(*part);
      P.Origin = 3;

      // Check which MCTruth this decay belongs to
      for(size_t i_h=0;i_h<theTruth.PrimaryKaonM.size();i_h++){
         SimParticle H = theTruth.PrimaryKaonM.at(i_h);
         if(PosMatch(TVector3(P.StartX,P.StartY,P.StartZ),TVector3(H.EndX,H.EndY,H.EndZ))){
            P.MCTruthIndex = H.MCTruthIndex;            
            theTruth.DecayVertex_X[P.MCTruthIndex] = P.StartX;
            theTruth.DecayVertex_Y[P.MCTruthIndex] = P.StartY;
            theTruth.DecayVertex_Z[P.MCTruthIndex] = P.StartZ;
         }
      } 

      theTruth.KaonMDecay.push_back(P);
   }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetKaon0Decay(){

   for(size_t i_d=0;i_d<Kaon0_Daughter_IDs.size();i_d++){

      // Geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
      if(partByID.find(Kaon0_Daughter_IDs[i_d]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[Kaon0_Daughter_IDs[i_d]];
       
      SimParticle P = MakeSimParticle(*part);
      P.Origin = 4;

      // Check if the parent is a K0S/K0L
      bool child_of_K0SL=false;
      if(partByID.find(part->Mother()) != partByID.end()){
        art::Ptr<simb::MCParticle> part2 = partByID[part->Mother()];
        if(part2->PdgCode() == 130 || part2->PdgCode() == 310){ 
          P.Origin = 9;
          child_of_K0SL=true;
        }
      }

      // depending on whether this particle is the child of a K+/- or K0, parents will live in different vectors
      std::vector<SimParticle>* decay_parents = child_of_K0SL ? &theTruth.NeutralKaonDecayK0SL : &theTruth.PrimaryKaon;

      for(size_t i_k=0;i_k<decay_parents->size();i_k++){
         SimParticle K = decay_parents->at(i_k);
         if(PosMatch(TVector3(P.StartX,P.StartY,P.StartZ),TVector3(K.EndX,K.EndY,K.EndZ))) P.MCTruthIndex = K.MCTruthIndex;
      } 
  
      theTruth.Kaon0Decay.push_back(P);     
   }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetNeutralKaonDecay(){

   for(size_t i_d=0;i_d<NeutralKaon_Daughter_IDs.size();i_d++){

      // Geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
      if(partByID.find(NeutralKaon_Daughter_IDs[i_d]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[NeutralKaon_Daughter_IDs[i_d]];

      SimParticle P = MakeSimParticle(*part);
      P.Origin = 5;

      // Check which MCTruth this decay belongs to
      for(size_t i_h=0;i_h<theTruth.PrimaryKaon.size();i_h++){
         SimParticle H = theTruth.PrimaryKaon.at(i_h);
         if(PosMatch(TVector3(P.StartX,P.StartY,P.StartZ),TVector3(H.EndX,H.EndY,H.EndZ))) P.MCTruthIndex = H.MCTruthIndex;
      } 

      theTruth.NeutralKaonDecayK0SL.push_back(P);

   }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::GetHyperonDecay(){

   if(!Hyperon_Daughter_IDs.size()) return;

   for(size_t i_d=0;i_d<Hyperon_Daughter_IDs.size();i_d++){

      // Geant does not always keep all particles it simulates, first check daughter is actually in list of IDs
      if(partByID.find(Hyperon_Daughter_IDs[i_d]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[Hyperon_Daughter_IDs[i_d]];

      SimParticle P = MakeSimParticle(*part);
      P.Origin = 6;

      // Check which MCTruth this decay belongs to
      for(size_t i_h=0;i_h<theTruth.PrimaryHyperon.size();i_h++){
         SimParticle H = theTruth.PrimaryHyperon.at(i_h);
         if(PosMatch(TVector3(P.StartX,P.StartY,P.StartZ),TVector3(H.EndX,H.EndY,H.EndZ))){
            P.MCTruthIndex = H.MCTruthIndex;            
            theTruth.DecayVertex_X[P.MCTruthIndex] = P.StartX;
            theTruth.DecayVertex_Y[P.MCTruthIndex] = P.StartY;
            theTruth.DecayVertex_Z[P.MCTruthIndex] = P.StartZ;
         }
      } 

      theTruth.HyperonDecay.push_back(P);     
   }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<int> SubModuleG4Truth::GetChildIDs(const art::Ptr<simb::MCParticle> &g4p,bool IsNeutron){

   std::vector<int> DecayProduct_IDs;

   //if( (g4p->EndProcess() != "Decay" || g4p->EndProcess() != "FastScintillation") && !IsNeutron && !isKaon(g4p->PdgCode())) return DecayProduct_IDs;

   for(int i_d=0;i_d<g4p->NumberDaughters();i_d++){

      if(partByID.find(g4p->Daughter(i_d)) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[g4p->Daughter(i_d)];

      if(part->PdgCode() > 10000) continue;

      double X = part->Position().X();
      double Y = part->Position().Y();
      double Z = part->Position().Z();

      double EndX = g4p->EndPosition().X();
      double EndY = g4p->EndPosition().Y();
      double EndZ = g4p->EndPosition().Z();

      if(!PosMatch(TVector3(X,Y,Z),TVector3(EndX,EndY,EndZ))) continue;

      DecayProduct_IDs.push_back(g4p->Daughter(i_d));
   }

   return DecayProduct_IDs;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

int SubModuleG4Truth::GetOrigin(int trackid){

   if(std::find(Primary_IDs.begin(),Primary_IDs.end(),trackid) != Primary_IDs.end()) return 1;
   else if(std::find(KaonP_Daughter_IDs.begin(),KaonP_Daughter_IDs.end(),trackid) != KaonP_Daughter_IDs.end()) return 2;
   else if(std::find(KaonM_Daughter_IDs.begin(),KaonM_Daughter_IDs.end(),trackid) != KaonM_Daughter_IDs.end()) return 3;
   else if(std::find(Kaon0_Daughter_IDs.begin(),Kaon0_Daughter_IDs.end(),trackid) != Kaon0_Daughter_IDs.end()) {
     if(std::find(NeutralKaon_Daughter_IDs.begin(),NeutralKaon_Daughter_IDs.end(),partByID[trackid]->Mother()) != NeutralKaon_Daughter_IDs.end()) return 9;
     return 4;
   }
   else if(std::find(Hyperon_Daughter_IDs.begin(),Hyperon_Daughter_IDs.end(),trackid) != Hyperon_Daughter_IDs.end()) return 6;
   else return 10;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

bool SubModuleG4Truth::FindKaonPScatter(){

   std::vector<int> KaonP_IDs;

   for(size_t i=0;i<Primary_IDs.size();i++){

      if(partByID.find(Primary_IDs[i]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[Primary_IDs[i]];

      if(part->PdgCode() != 321) continue;

      if(part->EndProcess() == "kaon+Inelastic") KaonP_IDs.push_back(part->TrackId()); 

   }

   for(size_t i=0;i<KaonP_IDs.size();i++){
      art::Ptr<simb::MCParticle> part = partByID[KaonP_IDs[i]];
      std::vector<int> KaonP_Inelastic_Daughter_IDs = GetChildIDs(part,true);      

      int nKaonPs = 0;

      for(size_t i_d=0;i_d<KaonP_Inelastic_Daughter_IDs.size();i_d++){

         // Look for kaons
         art::Ptr<simb::MCParticle> part2 = partByID[KaonP_Inelastic_Daughter_IDs[i_d]];
         //double P = sqrt(part2->Px()*part2->Px() + part2->Py()*part2->Py() + part2->Pz()*part2->Pz());
         if(part2->PdgCode() == 321) nKaonPs++; 
      }

      if(nKaonPs >= 1) return true;
   }

   return false;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

bool SubModuleG4Truth::FindProtonScatter(){

   std::vector<int> Proton_IDs;

   for(size_t i=0;i<Primary_IDs.size();i++){

      if(partByID.find(Primary_IDs[i]) == partByID.end()) continue;

      art::Ptr<simb::MCParticle> part = partByID[Primary_IDs[i]];

      if(part->PdgCode() != 2212) continue;

      if(part->EndProcess() == "protonInelastic") Proton_IDs.push_back(part->TrackId()); 

   }

   for(size_t i=0;i<Proton_IDs.size();i++){
      art::Ptr<simb::MCParticle> part = partByID[Proton_IDs[i]];
      std::vector<int> Proton_Inelastic_Daughter_IDs = GetChildIDs(part,true);      

      int nProtons = 0;

      for(size_t i_d=0;i_d<Proton_Inelastic_Daughter_IDs.size();i_d++){

         // Look for protons
         art::Ptr<simb::MCParticle> part2 = partByID[Proton_Inelastic_Daughter_IDs[i_d]];
         //double P = sqrt(part2->Px()*part2->Px() + part2->Py()*part2->Py() + part2->Pz()*part2->Pz());
         if(part2->PdgCode() == 2212) nProtons++; 
      }

      if(nProtons >= 1) return true;
   }

   return false;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::SetDecayThresholds(double decayprotonthresh,double decaypionthresh){

   DecayProtonThresh = decayprotonthresh;
   DecayPionThresh = decaypionthresh;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

bool SubModuleG4Truth::PosMatch(TVector3 Pos1,TVector3 Pos2){

   return (Pos1-Pos2).Mag() < _EPSILON_;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::MCTruthMatch(SimParticle &P){
 
  for(size_t i_mct=0;i_mct<PrimaryVertices.size();i_mct++)
    if(PosMatch(TVector3(P.StartX,P.StartY,P.StartZ),PrimaryVertices.at(i_mct)))
      P.MCTruthIndex = i_mct;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::MCTruthMatch(SimParticle &P,int trackid){

  // Go up the hierarchy of particles until you reach a direct child of the neutrino

  // If particle is missing from map
  if(partByID.find(trackid) == partByID.end()) {
    P.MCTruthIndex = -1;
    return;
  }
  
  art::Ptr<simb::MCParticle> mcp = partByID.at(trackid);
  while(true){
    if(partByID.find(mcp->Mother()) == partByID.end() || mcp->Mother() == 0) break;
    mcp = partByID.at(mcp->Mother());
  }

  for(size_t i_mct=0;i_mct<PrimaryVertices.size();i_mct++)
    if(PosMatch(TVector3(mcp->Position().X(),mcp->Position().Y(),mcp->Position().Z()),PrimaryVertices.at(i_mct)))
      P.MCTruthIndex = i_mct;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleG4Truth::SetFlags(){

   // Iterate over the different MCTruths in the event, set their flags
   for(int i_t=0;i_t<NMCTruths;i_t++){

      theTruth.InActiveTPC[i_t] = inActiveTPC(PrimaryVertices.at(i_t));

      theTruth.IsHyperon[i_t] = false;
      theTruth.IsKaon[i_t] = false; 
      theTruth.IsKaonP[i_t] = false;
      theTruth.IsKaonP_NuMuP[i_t] = false; 
      theTruth.IsKaonP_PiPPi0[i_t] = false;
      theTruth.IsKaonP_2PiPPiM[i_t] = false;
      theTruth.IsKaonP_ENuE[i_t] = false;
      theTruth.IsKaonP_2PiNPiP[i_t] = false;
      theTruth.IsKaonP_Others[i_t] = false;
      theTruth.IsKaonM[i_t] = false; 
      theTruth.IsKaon0[i_t] = false; 
      theTruth.IsAssociatedKaonP[i_t] = false;

      int nHyperons=0, nKaonPs=0;

      for(size_t i_k=0;i_k<theTruth.PrimaryKaonP.size();i_k++){
         if(theTruth.PrimaryKaonP.at(i_k).MCTruthIndex == i_t){
            nKaonPs++;
            theTruth.IsKaonP[i_t] = true;
         }
      }

      for(size_t i_k=0;i_k<theTruth.PrimaryKaon.size();i_k++)
         if(theTruth.PrimaryKaon.at(i_k).MCTruthIndex == i_t) theTruth.IsKaon[i_t] = true;

      for(size_t i_k=0;i_k<theTruth.PrimaryKaonM.size();i_k++)
         if(theTruth.PrimaryKaonM.at(i_k).MCTruthIndex == i_t) theTruth.IsKaonM[i_t] = true;

      for(size_t i_k=0;i_k<theTruth.PrimaryKaon.size();i_k++)
	if(theTruth.PrimaryKaon.at(i_k).MCTruthIndex == i_t) theTruth.IsKaon0[i_t] = true;


      int nProducts=0;
      bool hasMuonP=false, hasElectronP=false, hasNuMu=false, hasNuE=false, hasPionP=false, hasPionM=false, hasPion0=false;
      for(size_t i_d=0;i_d<theTruth.KaonPDecay.size();i_d++){
         if(theTruth.KaonPDecay.at(i_d).MCTruthIndex == i_t){
            nProducts++;
            if(theTruth.KaonPDecay.at(i_d).PDG ==  -13) hasMuonP     = true;  
            if(theTruth.KaonPDecay.at(i_d).PDG ==  -11) hasElectronP = true;
            if(theTruth.KaonPDecay.at(i_d).PDG ==   12) hasNuE       = true;
            if(theTruth.KaonPDecay.at(i_d).PDG ==   14) hasNuMu      = true;  
            if(theTruth.KaonPDecay.at(i_d).PDG ==  211) hasPionP     = true;
            if(theTruth.KaonPDecay.at(i_d).PDG == -211) hasPionM     = true;
            if(theTruth.KaonPDecay.at(i_d).PDG ==  111) hasPion0     = true;
         }
      }

      if(nKaonPs == 1){
	if( nProducts == 2 && hasMuonP && hasNuMu ) theTruth.IsKaonP_NuMuP[i_t] = true;
	else if( nProducts == 2 && hasPionP && hasPion0 ) theTruth.IsKaonP_PiPPi0[i_t] = true;
	else if( nProducts == 3 && hasPionP && hasPionM ) theTruth.IsKaonP_2PiPPiM[i_t] = true;
	else if( nProducts == 3 && hasElectronP && hasNuE ) theTruth.IsKaonP_ENuE[i_t] = true;
	else if( nProducts == 3 && hasPion0 && hasPionP ) theTruth.IsKaonP_2PiNPiP[i_t] = true;
	else theTruth.IsKaonP_Others[i_t] = true;

      }
     
      //if(nHyperons == 1 && theTruth.IsLambda.at(i_t) && nProducts == 2 && hasProton && hasPion) theTruth.IsLambdaCharged[i_t] = true;
      //if(nHyperons == 1 && theTruth.IsSigmaZero.at(i_t) && nProducts == 2 && hasProton && hasPion) theTruth.IsSigmaZeroCharged[i_t] = true;

      for(size_t i_h=0;i_h<theTruth.PrimaryHyperon.size();i_h++){
        if(theTruth.PrimaryHyperon.at(i_h).MCTruthIndex == i_t){
	  nHyperons++;
	  theTruth.IsHyperon[i_t] = true;
	}
      }
   
      /*
      for(SimParticle kaon : theTruth.NeutralKaonDecayK0SL)        
        if(kaon.MCTruthIndex == i_t && kaon.PDG == 310) theTruth.IsK0S[i_t] = true;

      if(theTruth.IsK0S[i_t]){
        bool hasPiP=false,hasPiM=false;
        for(SimParticle kaondecay : theTruth.KaonDecay){
          if(kaondecay.MCTruthIndex == i_t && kaondecay.Origin == 7){
            if(kaondecay.PDG == 211) hasPiP = true;  
            if(kaondecay.PDG == -211) hasPiM = true;  
          }
        }
        if(hasPiP && hasPiM) theTruth.IsK0SCharged[i_t] = true; 
      }
      */
        
      // If there are multiple hyperons/antihyperons or hyperons and kaons present together
      // flag as associated hyperon
      if(nHyperons > 1 || (nHyperons >= 1 && nKaonPs >= 1)) theTruth.IsAssociatedKaonP[i_t] = true;

   }

   theTruth.EventHasKaonPScatter = FindKaonPScatter();   
   theTruth.EventHasProtonScatter = FindProtonScatter();   
   theTruth.EventHasHyperon = std::find(theTruth.IsHyperon.begin(), theTruth.IsHyperon.end(), true) != theTruth.IsHyperon.end();
   theTruth.EventHasKaon = std::find(theTruth.IsKaon.begin(), theTruth.IsKaon.end(), true) != theTruth.IsKaon.end();
   theTruth.EventHasKaonP = std::find(theTruth.IsKaonP.begin(), theTruth.IsKaonP.end(), true) != theTruth.IsKaonP.end();
   theTruth.EventHasKaonP_NuMuP = std::find(theTruth.IsKaonP_NuMuP.begin(), theTruth.IsKaonP_NuMuP.end(), true) != theTruth.IsKaonP_NuMuP.end();
   theTruth.EventHasKaonP_PiPPi0 = std::find(theTruth.IsKaonP_PiPPi0.begin(), theTruth.IsKaonP_PiPPi0.end(), true) != theTruth.IsKaonP_PiPPi0.end();
   theTruth.EventHasKaonM = std::find(theTruth.IsKaonM.begin(), theTruth.IsKaonM.end(), true) != theTruth.IsKaonM.end();
   theTruth.EventHasKaon0 = std::find(theTruth.IsKaon0.begin(), theTruth.IsKaon0.end(), true) != theTruth.IsKaon0.end();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
