#ifndef RecoVariables_cxx_
#define _RecoVariables_cxx_

#include "SubModuleReco.h"

using namespace cckaon;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

//SubModuleReco::SubModuleReco(){}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleReco::SubModuleReco(art::Event const& e,bool isdata,fhicl::ParameterSet pset,bool particlegunmode, bool withrecoalg) :
SubModuleReco(e,isdata,
                  pset.get<std::string>("PFParticleModuleLabel"),
                  pset.get<std::string>("TrackModuleLabel"),
                  pset.get<std::string>("TrackRebuiltModuleLabel"),
                  pset.get<std::string>("ShowerModuleLabel"),
                  pset.get<std::string>("VertexModuleLabel"),
                  pset.get<std::string>("PIDModuleLabel"),
                  pset.get<std::string>("CaliPIDModuleLabel"),
                  pset.get<std::string>("CaloModuleLabel"),
                  pset.get<std::string>("HitModuleLabel"),
                  pset.get<std::string>("HitTruthAssnLabel"),
                  pset.get<std::string>("TrackHitAssnLabel"),
		  pset.get<std::string>("TrackRebuiltHitAssnLabel"),
	          pset.get<std::string>("RerunPIDModuleLabel"),
	          pset.get<std::string>("RerunCaliPIDModuleLabel"),
	          pset.get<std::string>("RerunCaloModuleLabel"),
                  pset.get<std::string>("ShowerHitAssnLabel"),
                  pset.get<std::string>("MetadataModuleLabel"),
                  pset.get<std::string>("GeneratorModuleLabel"),
                  pset.get<std::string>("G4ModuleLabel"),
                  pset.get<bool>("DoGetPIDs",true),
                  pset.get<bool>("IncludeCosmics",false),
	      particlegunmode,
	      withrecoalg)
{
  this->WithRecoAlgorithm = withrecoalg;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleReco::SubModuleReco(art::Event const& e,bool isdata,string pfparticlelabel,string tracklabel, string trackrebuiltlabel,
			     string showerlabel,string vertexlabel,string pidlabel,string calipidlabel,string calolabel,string hitlabel,
			     string hittruthassnlabel,string trackhitassnlabel,string trackrebuilthitassnlabel,
			     string rerunpidlabel, string reruncalipidlabel,string reruncalolabel, string showerhitassnlabel,string metadatalabel,string genlabel,
			     string g4label,bool dogetpids,bool includecosmics,bool particlegunmode, bool withrecoalg) :
PIDCalc(),
DoGetPIDs(dogetpids),
IncludeCosmics(includecosmics),
ParticleGunMode(particlegunmode)
{

   IsData = isdata;

   if(!e.getByLabel(pfparticlelabel,Handle_PFParticle)) 
      throw cet::exception("SubModuleReco") << "No PFParticle Data Products Found!" << std::endl;

   if(!e.getByLabel(tracklabel,Handle_Track)) 
      throw cet::exception("SubModuleReco") << "No Track Data Products Found!" << std::endl;

   if(withrecoalg){
     if(!e.getByLabel(trackrebuiltlabel,Handle_TrackRebuilt)) 
       throw cet::exception("SubModuleReco") << "No Track Data Products Found!" << std::endl;
   }

   if(!e.getByLabel(showerlabel,Handle_Shower)) 
      throw cet::exception("SubModuleReco") << "No Shower Data Products Found!" << std::endl;

   if(!e.getByLabel(hitlabel,Handle_Hit)) 
      throw cet::exception("SubModuleReco") << "No Hit Data Products Found!" << std::endl;

   art::fill_ptr_vector(Vect_PFParticle,Handle_PFParticle);
   art::fill_ptr_vector(Vect_Track,Handle_Track);
   if(withrecoalg) art::fill_ptr_vector(Vect_TrackRebuilt,Handle_TrackRebuilt);
   art::fill_ptr_vector(Vect_Shower,Handle_Shower);
   art::fill_ptr_vector(Vect_Hit,Handle_Hit);

   Assoc_PFPMuon = new art::FindManyP<anab::T0>(Vect_PFParticle,e,"NuCCproducer");
   Assoc_PFParticleVertex = new art::FindManyP<recob::Vertex>(Vect_PFParticle,e,vertexlabel);    
   Assoc_PFParticleTrack = new art::FindManyP<recob::Track>(Vect_PFParticle,e,tracklabel);
   //if(withrecoalg) Assoc_PFParticleTrackRebuilt = new art::FindManyP<recob::Track>(Vect_PFParticle,e,trackrebuiltlabel); 
   Assoc_PFParticleShower = new art::FindManyP<recob::Shower>(Vect_PFParticle,e,showerlabel);    
   Assoc_PFParticleMetadata = new art::FindManyP<larpandoraobj::PFParticleMetadata>(Vect_PFParticle,e,metadatalabel);   
   Assoc_TrackHit = new  art::FindManyP<recob::Hit>(Vect_Track,e,trackhitassnlabel);
   if(withrecoalg) Assoc_TrackRebuiltHit = new  art::FindManyP<recob::Hit>(Vect_TrackRebuilt,e,trackrebuilthitassnlabel);
   Assoc_ShowerHit = new  art::FindManyP<recob::Hit>(Vect_Shower,e,showerhitassnlabel);
   ParticlesPerHit = new art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData>(Handle_Hit,e,hittruthassnlabel);

   if(DoGetPIDs){
      Assoc_TrackCalo = new art::FindManyP<anab::Calorimetry>(Vect_Track,e,calolabel);
      if(withrecoalg) Assoc_TrackRebuiltCalo = new art::FindManyP<anab::Calorimetry>(Vect_TrackRebuilt,e,reruncalolabel);
      Assoc_TrackPID = new art::FindManyP<anab::ParticleID>(Vect_Track,e,pidlabel);
      if(withrecoalg) Assoc_TrackRebuiltPID = new art::FindManyP<anab::ParticleID>(Vect_TrackRebuilt,e,rerunpidlabel);
      Assoc_TrackCaliPID = new art::FindManyP<anab::ParticleID>(Vect_Track,e,calipidlabel);
      if(withrecoalg) Assoc_TrackRebuiltCaliPID = new art::FindManyP<anab::ParticleID>(Vect_TrackRebuilt,e,reruncalipidlabel);
   }

   llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
   llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
   llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);

   llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
   llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
   llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);

   llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
   llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
   llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);

   llr_pid_calculator_kaon.set_dedx_binning(0, kaonproton_parameters.dedx_edges_pl_0);
   llr_pid_calculator_kaon.set_par_binning(0, kaonproton_parameters.parameters_edges_pl_0);
   llr_pid_calculator_kaon.set_lookup_tables(0, kaonproton_parameters.dedx_pdf_pl_0);

   llr_pid_calculator_kaon.set_dedx_binning(1, kaonproton_parameters.dedx_edges_pl_1);
   llr_pid_calculator_kaon.set_par_binning(1, kaonproton_parameters.parameters_edges_pl_1);
   llr_pid_calculator_kaon.set_lookup_tables(1, kaonproton_parameters.dedx_pdf_pl_1);

   llr_pid_calculator_kaon.set_dedx_binning(2, kaonproton_parameters.dedx_edges_pl_2);
   llr_pid_calculator_kaon.set_par_binning(2, kaonproton_parameters.parameters_edges_pl_2);
   llr_pid_calculator_kaon.set_lookup_tables(2, kaonproton_parameters.dedx_pdf_pl_2);

   if(!IsData){
      G4T = new SubModuleG4Truth(e,genlabel,g4label,ParticleGunMode);
      G4T->GetParticleLists();
   }

   m_PFPID_TrackIndex.clear(); 
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::PrepareInfo(){

   theData.RecoPrimaryVertex = GetPrimaryVertex();

   int ipfp=0;
   int itrk=0;
   int ishw=0;
   std::vector<art::Ptr<recob::PFParticle>> pfmuons(0); 
   std::vector<int> reco_nu_daughters_id(0);
   art::Ptr<recob::PFParticle> pfmuon;
   art::Ptr<recob::Track> trkmuon;

   // Find CC muon and daughters
   for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){
     
     //if((pfp->Parent() == neutrinoID && pfp->PdgCode()==13)
     //m_PFPID_TrackIndex[pfp->Self()] = P.Index;

     // look at particles with neutrino parent and one associated track 
     if (pfp->Parent() == neutrinoID && Assoc_PFParticleTrack->at(ipfp).size()==1) {
       
       art::Ptr<recob::Track> track = Assoc_PFParticleTrack->at(ipfp).front(); 
       reco_nu_daughters_id.push_back(track->ID());
       //RecoParticle PDummy = MakeRecoParticle(pfp);
       //m_PFPID_TrackIndex[pfp->Self()] = PDummy.Index; // store index for neutrino primary tracks   
       
       // CC muon has a T0 associated  
       if (Assoc_PFPMuon->at(ipfp).size()==1){
	 RecoParticle P = MakeRecoParticle(pfp);
	 P.IsCCMu = true;
	 pfmuons.push_back(pfp);
	 theData.CCMuTrack = P;
       }
     }
     ipfp++;
   }

   if(pfmuons.size()==1){
     pfmuon = pfmuons.front(); 
     trkmuon = Assoc_PFParticleTrack->at(pfmuon.key()).front(); 
   }

   for(const art::Ptr<recob::Track> &ptrk : Vect_Track){

     const recob::Track& trk = *ptrk;
     RecoParticle P = MakeRecoParticle(ptrk);

     if(trkmuon)
       if(trk.ID() == trkmuon->ID()) continue;
     
     if( std::find(reco_nu_daughters_id.begin(), reco_nu_daughters_id.end(), trk.ID()) != reco_nu_daughters_id.end() ){
       P.Parentage = 1;
       P.InNuSlice = true;
     }else{
       P.Parentage = 2;
     }
     
     /*
     if(trk.ID() == reco_nu_daughters_id.at(itrk)){
       P.Parentage = 1;
       P.InNuSlice = true;
     }
     else{
       //if(m_PFPID_TrackIndex.find(pfp->Parent()) != m_PFPID_TrackIndex.end()){
	 if( std::find(reco_nu_daughters_id.begin(), reco_nu_daughters_id.end(), trk.ID() != reco_nu_daughters_id.end()){
	 P.Parentage = 2; 
	 //P.ParentIndex = m_PFPID_TrackIndex[pfp->Parent()]; 
       
	 //if(WithRecoAlgorithm && Vect_TrackRebuilt.size()>0) GetRebuiltTrackData(pfp,P);//override with rebuilt tracks, if any
       }
       else P.Parentage = 3;
     }
     */
     

     if(P.InNuSlice){
       theData.TrackPrimaryDaughters.push_back(P);
       //m_PFPID_TrackIndex[pfp->Self()] = P.Index; // store index for neutrino primary tracks
     }
     else theData.TrackOthers.push_back(P);
     /*
     else{
       if(P.UseRebuilt==true) theData.TrackRebuiltOthers.push_back(P);
       else theData.TrackOthers.push_back(P);
     }
     */
        
     itrk++;
   }

   if(WithRecoAlgorithm){
     for(const art::Ptr<recob::Track> &ptrk : Vect_TrackRebuilt){

       //rebuild tracks are all neutrino granddaughter
       RecoParticle P = MakeRecoParticle(ptrk);

       P.Parentage = 2;
       P.UseRebuilt = true;

       theData.TrackRebuiltOthers.push_back(P);
     }
   }

   for(art::Ptr<recob::Shower> &pshw : Vect_Shower){

     //const recob::Shower& shw = *pshw;
     RecoParticle P = MakeRecoParticle(pshw);

     if(P.InNuSlice) theData.ShowerPrimaryDaughters.push_back(P);
     else theData.ShowerOthers.push_back(P);

     ishw++;
   }

   /*
   for(art::Ptr<recob::PFParticle> &trk : Vect_Track){
     std::cout << "this is loop over Vect_Track" << std::endl;
     RecoParticle PD;
     GetPIDs(trk,PD);
   }
   */

   /*
   int ipfp=0;

   for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){

     if(!IncludeCosmics && pfp->Parent() != neutrinoID && m_PFPID_TrackIndex.find(pfp->Parent()) == m_PFPID_TrackIndex.end()) continue; 

     RecoParticle P = MakeRecoParticle(pfp);

     // look at particles with neutrino parent and one associated track
     if (pfp->Parent() == neutrinoID && Assoc_PFParticleTrack->at(ipfp).size()==1) {
         
       // CC muon has a T0 associated
       if (Assoc_PFPMuon->at(ipfp).size()==1) {
	 //PFPMuon.push_back(pfp);
	 P.IsCCMu = true;
       }
     }
     ipfp++;
  


      
     if(pfp->Parent() == neutrinoID){
       P.Parentage = 1;
       P.InNuSlice = true;         
     }
     else{ 	      
       if(m_PFPID_TrackIndex.find(pfp->Parent()) != m_PFPID_TrackIndex.end()){ // has daughter track
	 P.Parentage = 2; 
	 P.ParentIndex = m_PFPID_TrackIndex[pfp->Parent()]; 
	 if(WithRecoAlgorithm && Vect_TrackRebuilt.size()>0) GetRebuiltTrackData(pfp,P);//override with rebuilt tracks, if any
       }
       else P.Parentage = 3;
     }
     
     if(P.PDG == 13){ // This is Pandora PDG code (11 or 13)
       if(P.InNuSlice){
	 theData.TrackPrimaryDaughters.push_back(P);
	 m_PFPID_TrackIndex[pfp->Self()] = P.Index; // store index for neutrino primary tracks
       }
       else{
	 if(P.UseRebuilt==true) theData.TrackRebuiltOthers.push_back(P);
	 else theData.TrackOthers.push_back(P);
       }
     }
     else if(P.PDG == 11){
       if(P.InNuSlice) theData.ShowerPrimaryDaughters.push_back(P);
       else theData.ShowerOthers.push_back(P);
     }
   }
*/

   theData.NPrimaryDaughters = theData.TrackPrimaryDaughters.size() + theData.ShowerPrimaryDaughters.size();
   theData.NPrimaryTrackDaughters = theData.TrackPrimaryDaughters.size();
   theData.NPrimaryShowerDaughters = theData.ShowerPrimaryDaughters.size();
   theData.NOtherTracks = theData.TrackOthers.size();
   theData.NOtherRebuiltTracks = theData.TrackRebuiltOthers.size();
   theData.NOtherShowers = theData.ShowerOthers.size();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

RecoData SubModuleReco::GetInfo(){
   return theData;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 SubModuleReco::GetPrimaryVertex(){ //returns neutrino vertex

   auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

   for(const art::Ptr<recob::PFParticle> &pfp : Vect_PFParticle){

      if(pfp->IsPrimary() && isNeutrino(pfp->PdgCode())){

         neutrinoID = pfp->Self();

         std::vector<art::Ptr<recob::Vertex>> pfpVertex = Assoc_PFParticleVertex->at(pfp.key());

         for(const art::Ptr<recob::Vertex> &vtx : pfpVertex){

            geo::Point_t point = { vtx->position().X() , vtx->position().Y() , vtx->position().Z() };
            geo::Vector_t sce_corr = SCE->GetPosOffsets(point);

            return TVector3(vtx->position().X() + sce_corr.X(),vtx->position().Y()-sce_corr.Y(),vtx->position().Z()-sce_corr.Z());
         }
      }
   }

   // If there is no neutrino candidate
   return TVector3(-1000,-1000,-1000);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

RecoParticle SubModuleReco::MakeRecoParticle(const art::Ptr<recob::PFParticle> &pfp){

   RecoParticle P;

   P.PDG = pfp->PdgCode();

   std::vector<art::Ptr<recob::Track>> pfpTracks = Assoc_PFParticleTrack->at(pfp.key());
   std::vector<art::Ptr<recob::Shower>> pfpShowers = Assoc_PFParticleShower->at(pfp.key());

   if(pfp->PdgCode() == 13 && pfpTracks.size() != 1) P.PDG = 0; // how to handle scattered particles
   if(pfp->PdgCode() == 11 && pfpShowers.size() != 1) P.PDG = 0;

   GetPFPMetadata(pfp,P);

   if(pfpTracks.size() == 1){
      GetTrackData(pfp,P);
      GetVertexData(pfp,P);
   }

   if(pfpShowers.size() == 1){
     GetShowerData(pfp,P);
     GetVertexData(pfp,P);
   }

   return P;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

RecoParticle SubModuleReco::MakeRecoParticle(const art::Ptr<recob::Track> &trk){

   RecoParticle P;

   P.PDG = 13;
   GetTrackData(trk,P);

   //this is same calculation as GetVertexData(pfp,P);
   const recob::Track& track = *trk;
   TVector3 vtx(track.Vertex().X(), track.Vertex().Y(), track.Vertex().Z());
   P.SetVertex(vtx);
   P.Displacement = (vtx-theData.RecoPrimaryVertex).Mag();

   return P;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

RecoParticle SubModuleReco::MakeRecoParticle(const art::Ptr<recob::Shower> &shw){

   RecoParticle P;

   P.PDG = 11;
   GetShowerData(shw,P);

   //this is same calculation as GetVertexData(pfp,P);
   const recob::Shower& shower = *shw;
   TVector3 vtx(shower.ShowerStart().X(), shower.ShowerStart().Y(), shower.ShowerStart().Z());
   P.SetVertex(vtx);
   //P.Displacement = (vtx-theData.RecoPrimaryVertex).Mag();

   return P;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetPFPMetadata(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P){

   std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpMeta = Assoc_PFParticleMetadata->at(pfp.key());

   for(const art::Ptr<larpandoraobj::PFParticleMetadata> &meta : pfpMeta){

      const larpandoraobj::PFParticleMetadata::PropertiesMap &pfParticlePropertiesMap(meta->GetPropertiesMap());

      if (!pfParticlePropertiesMap.empty()){
         for (larpandoraobj::PFParticleMetadata::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it){
            if(it->first == "TrackScore"){
               P.TrackShowerScore = it->second;
            }
         }
      }	
   }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetTrackData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P){

  std::vector<art::Ptr<recob::Track>> pfpTracks;
  std::vector<art::Ptr<recob::Hit>> hits;
  art::Ptr<recob::Track> trk;
  
  pfpTracks = Assoc_PFParticleTrack->at(pfp.key());
  trk = pfpTracks.at(0); 
  hits = Assoc_TrackHit->at(trk.key());
  
   if(pfpTracks.size() != 1) return;

   // Sets track length/position related variables
   SetTrackVariables(P,trk);

   if(!IsData) TruthMatch(trk,P);

   if(!IsData) MergeCheck(hits,P);

   if(DoGetPIDs) GetPIDs(trk,P);
   
   theData.TrackStarts.push_back(TVector3(trk->Start().X(),trk->Start().Y(),trk->Start().Z()));
   P.Index = theData.TrackStarts.size() - 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetTrackData(const art::Ptr<recob::Track> &trk,RecoParticle &P){

  std::vector<art::Ptr<recob::Hit>> hits;
  
  hits = Assoc_TrackHit->at(trk.key());
  if(WithRecoAlgorithm) hits = Assoc_TrackRebuiltHit->at(trk.key());
  
   // Sets track length/position related variables
   SetTrackVariables(P,trk);

   if(!IsData) TruthMatch(trk,P);

   if(!IsData) MergeCheck(hits,P);

   if(DoGetPIDs) GetPIDs(trk,P);
   
   theData.TrackStarts.push_back(TVector3(trk->Start().X(),trk->Start().Y(),trk->Start().Z()));
   P.Index = theData.TrackStarts.size() - 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetRebuiltTrackData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P){

  std::vector<art::Ptr<recob::Track>> pfpTracks;
  std::vector<art::Ptr<recob::Hit>> hits;
  art::Ptr<recob::Track> trk;

  bool hasDaughterTrack=false;
  double trk_usual_length=-1;

  std::vector<art::Ptr<recob::Track>> pfpTracks_usual = Assoc_PFParticleTrack->at(pfp.key());

  if(pfpTracks_usual.size()>0){

    art::Ptr<recob::Track> ptrk_usual = pfpTracks_usual.at(0); 
    const recob::Track& trk_usual = *ptrk_usual;
    TVector3 vtx_usual(P.X, P.Y, P.Z);
    trk_usual_length=trk_usual.Length();
    hasDaughterTrack=true;

  }
    
  double maxLength=-1;
  int index=-1;
  
  if(!Vect_TrackRebuilt.size()) return;

  for (unsigned int i_rebtrk=0; i_rebtrk < Vect_TrackRebuilt.size(); i_rebtrk++){

    art::Ptr<recob::Track> ptrk_candidate(Handle_TrackRebuilt,i_rebtrk);
    const recob::Track& trk_candidate = * ptrk_candidate;

    if(trk_candidate.Length()>maxLength)
      index = i_rebtrk;
  }

  art::Ptr<recob::Track> ptrk_rebuilt(Handle_TrackRebuilt,index);
  const recob::Track& trk_rebuilt = * ptrk_rebuilt;
  TVector3 vtx_rebuilt(trk_rebuilt.Vertex().X(), trk_rebuilt.Vertex().Y(), trk_rebuilt.Vertex().Z());
 
  //double vtxDistance = (vtx_rebuilt - vtx_usual).Mag();
  //if(vtxDistance>10) return;//vtx distance
  
  //store P with rebuilt track information
  //if( trk_rebuilt.Length()>0. && ( trk_usual.Length()<40. ||  trk_usual.Length()>65.) ){

  if( (hasDaughterTrack==true && ( trk_rebuilt.Length()>0. && ( trk_usual_length<40. ||  trk_usual_length>65.))) ||
      hasDaughterTrack==false){

    trk = ptrk_rebuilt;
    hits = Assoc_TrackRebuiltHit->at(trk.key());
    P.UseRebuilt = true;

    if(!trk) return;
    P.PDG = 13;
   
    // Sets track length/position related variables
    SetTrackVariables(P,trk);

    //update vector
    P.SetVertex(vtx_rebuilt);
    P.Displacement = (vtx_rebuilt-theData.RecoPrimaryVertex).Mag(); 
    
    if(!IsData) TruthMatch(trk,P);
    
    //if(!IsData && WithRecoAlgorithm) MergeCheck(hits,P);
    
    if(DoGetPIDs) GetPIDs(trk,P);
    
    theData.TrackStarts.pop_back();
    theData.TrackStarts.push_back(TVector3(trk->Start().X(),trk->Start().Y(),trk->Start().Z()));
    P.Index = theData.TrackStarts.size() - 1;

  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetShowerData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P){

   std::vector<art::Ptr<recob::Shower>> pfpShowers = Assoc_PFParticleShower->at(pfp.key());

   if(pfpShowers.size() != 1) return;

   art::Ptr<recob::Shower> shw = pfpShowers.at(0);

   // Sets track length/position related variables
   //SetTrackVariables(P,shw);

   //if(!IsData) TruthMatch(trk,P);

   std::vector<art::Ptr<recob::Hit>> hits = Assoc_ShowerHit->at(shw.key());
   if(!IsData) MergeCheck(hits,P);

   //if(DoGetPIDs) GetPIDs(trk,P);
   
   //theData.TrackStarts.push_back(TVector3(trk->Start().X(),trk->Start().Y(),trk->Start().Z()));
   //P.Index = theData.TrackStarts.size() - 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetShowerData(const art::Ptr<recob::Shower> &shw,RecoParticle &P){

  std::vector<art::Ptr<recob::Hit>> hits;

  hits = Assoc_ShowerHit->at(shw.key());

  if(!IsData) MergeCheck(hits,P);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::TruthMatch(const art::Ptr<recob::Track> &trk,RecoParticle &P){

  std::vector<art::Ptr<recob::Hit>> hits;
  if(P.UseRebuilt==true) hits = Assoc_TrackRebuiltHit->at(trk.key());
  else hits = Assoc_TrackHit->at(trk.key());

   std::unordered_map<int,double>  trkide;
   int maxhits=-1;

   simb::MCParticle const* matchedParticle = NULL;

   std::vector<simb::MCParticle const*> particleVec;
   std::vector<anab::BackTrackerHitMatchingData const*> matchVec;

   for(size_t i_hit=0;i_hit<hits.size();++i_hit){

      particleVec.clear();
      matchVec.clear();
      ParticlesPerHit->get(hits[i_hit].key(),particleVec,matchVec);

      for(size_t i_particle=0;i_particle<particleVec.size();++i_particle){

         trkide[particleVec[i_particle]->TrackId()]++; 

         //new method - choose particle depositing energy in the most hits
         if(trkide[particleVec[i_particle]->TrackId()] > maxhits){
            maxhits = trkide[particleVec[i_particle]->TrackId()];
            matchedParticle = particleVec[i_particle];
         }
      }
   }

   if(matchedParticle != NULL){ 

      SimParticle SP = MakeSimParticle(*matchedParticle);
            
      SP.Origin = G4T->GetOrigin(matchedParticle->TrackId());
      G4T->MCTruthMatch(SP,matchedParticle->TrackId());
 
      P.HasTruth = true;
      P.MCTruthIndex = SP.MCTruthIndex;
      P.TrackTruePDG = SP.PDG;
      P.TrackTrueE = SP.E;
      P.TrackTruePx = SP.Px;
      P.TrackTruePy = SP.Py;
      P.TrackTruePz = SP.Pz;
      P.TrackTrueEndE = SP.E;
      P.TrackTrueEndPx = SP.EndPx;
      P.TrackTrueEndPy = SP.EndPy;
      P.TrackTrueEndPz = SP.EndPz;
      P.TrackTrueModMomentum = SP.ModMomentum;
      P.TrackTrueEndModMomentum = SP.EndModMomentum;
      P.TrackTrueKE = SP.KE;
      P.TrackTrueEndKE = SP.EndKE;
      P.TrackTrueLength = SP.Travel;
      P.TrackTrueOrigin = SP.Origin;
      P.TrackTruthPurity = (double)maxhits/hits.size();
   }
   else P.HasTruth = false;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::MergeCheck(const std::vector<art::Ptr<recob::Hit>>& hits, RecoParticle &P){

  std::map<int, double> trkide;
  std::map<int, double> trkidhit;
  std::map<int, int> trkidpdg;
  std::map<int, int> trkidmrgid;
  std::map<int, int> trkidmother;

  std::vector<simb::MCParticle const*> particle_vec;
  std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
  double totalenergy = 0;
  int totalhits = 0;

   for(size_t i_hit=0;i_hit<hits.size();++i_hit){

      particle_vec.clear();
      match_vec.clear();
      ParticlesPerHit->get(hits[i_hit].key(), particle_vec, match_vec);

      for(size_t i_p=0;i_p<particle_vec.size();++i_p){ 

         int trackID = particle_vec[i_p]->TrackId();
         trkide[trackID] += match_vec[i_p]->energy;
         totalenergy += match_vec[i_p]->energy;
         trkidhit[trackID]++;
	 totalhits++;
         trkidpdg[trackID] = particle_vec[i_p]->PdgCode();
         trkidmrgid[trackID] = -9;
         trkidmother[trackID] = particle_vec[i_p]->Mother();

      }
   }

   // Generate merged IDs for related tracks
  int currentMergedID = 1;
  for (auto& pdg_entry : trkidpdg) {

    if (trkidmrgid[pdg_entry.first] != -9) continue;
    trkidmrgid[pdg_entry.first] = currentMergedID;
    int currentMotherTrackID = trkidmother[pdg_entry.first];

    while (currentMotherTrackID > 0) {
      if(trkidpdg.find(currentMotherTrackID) == trkidpdg.end()) break;
      if(trkidpdg[currentMotherTrackID] != pdg_entry.second) break;
      trkidmrgid[currentMotherTrackID] = currentMergedID;
      currentMotherTrackID = trkidmother[currentMotherTrackID];
    }
    ++currentMergedID;
  }
  
  // Merge tracks with the same merged ID
  for (auto& mrgid1 : trkidmrgid) {
    for (auto& mrgid2 : trkidmrgid) {
      if (mrgid1.first == mrgid2.first || mrgid1.second != mrgid2.second) continue;
      trkide[mrgid1.first] += trkide[mrgid2.first];
      trkidhit[mrgid1.first] += trkidhit[mrgid2.first];
      trkidmrgid[mrgid2.first] = -1;  // Mark as merged
    }
  }

  // Prepare data for output
  std::vector<std::pair<int, double>> energySorted, hitsSorted;
  for (auto const& entry : trkide) {
    if (trkidmrgid[entry.first] != -1) {
      energySorted.emplace_back(trkidpdg[entry.first], entry.second);
      hitsSorted.emplace_back(trkidpdg[entry.first], trkidhit[entry.first]);
    }
  }

  std::sort(energySorted.begin(), energySorted.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
      return a.second > b.second;  // Descending order by energy
    });

  std::sort(hitsSorted.begin(), hitsSorted.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
      return a.second > b.second;  // Descending order by hits
    });

  std::vector<int> MergePDG(3, -1);
  std::vector<double> MergeEnergyPurity(3, 0.0), MergeHitPurity(3, 0.0);
 
  for (size_t i = 0; i < 3 && i < energySorted.size(); ++i) {
    MergePDG[i] = energySorted[i].first;
    MergeEnergyPurity[i] = energySorted[i].second / totalenergy;
    MergeHitPurity[i] = hitsSorted[i].second / totalhits;
  }

  P.SetMergeCheck(MergePDG, MergeEnergyPurity, MergeHitPurity);
  
 }

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetPIDs(const art::Ptr<recob::Track> &trk,RecoParticle &P){

  std::vector<art::Ptr<anab::Calorimetry>> caloFromTrack;
  std::vector<art::Ptr<anab::ParticleID>> trackPID;
  std::vector<art::Ptr<anab::ParticleID>> trackCaliPID;

   if(P.UseRebuilt==true){
	caloFromTrack = Assoc_TrackRebuiltCalo->at(trk.key());
   	trackPID = Assoc_TrackRebuiltPID->at(trk.key());
   	trackCaliPID = Assoc_TrackRebuiltCaliPID->at(trk.key());
   }
   else{
	caloFromTrack = Assoc_TrackCalo->at(trk.key());
   	trackPID = Assoc_TrackPID->at(trk.key());
   	trackCaliPID = Assoc_TrackCaliPID->at(trk.key());
   }

   std::vector<anab::sParticleIDAlgScores> AlgScoresVecPID = trackPID.at(0)->ParticleIDAlgScores();
   std::vector<anab::sParticleIDAlgScores> AlgScoresVecCaliPID = trackCaliPID.at(0)->ParticleIDAlgScores();

   PIDStore store = PIDCalc.GetPIDs(trk,caloFromTrack,AlgScoresVecPID,AlgScoresVecCaliPID);

   P.Track_LLR_PID = store.LLR;
   P.Track_LLR_PID_Kaon = store.LLR_Kaon;
   P.Track_LLR_PID_Kaon_Partial = store.LLR_Kaon_Partial;
   P.MeandEdX_Plane0 = store.MeandEdX_Plane0;
   P.MeandEdX_Plane1 = store.MeandEdX_Plane1;
   P.MeandEdX_Plane2 = store.MeandEdX_Plane2;
   P.MeandEdX_ThreePlane = store.MeandEdX_3Plane;
   P.Track_Bragg_PID_Kaon = store.Bragg_Kaon_3Plane;
   P.Track_Chi2_Kaon_Plane0 = store.Chi2_Kaon.at(0);
   P.Track_Chi2_Kaon_Plane1 = store.Chi2_Kaon.at(1);
   P.Track_Chi2_Kaon_Plane2 = store.Chi2_Kaon.at(2);
   P.Track_Chi2_Kaon_3Plane = store.Chi2_Kaon_3Plane;
   P.Track_Chi2_Proton_Plane0 = store.Chi2_Proton.at(0);
   P.Track_Chi2_Proton_Plane1 = store.Chi2_Proton.at(1);
   P.Track_Chi2_Proton_Plane2 = store.Chi2_Proton.at(2);
   P.Track_Chi2_Proton_3Plane = store.Chi2_Proton_3Plane;
   P.Track_Chi2_Pion_Plane0 = store.Chi2_Pion.at(0);
   P.Track_Chi2_Pion_Plane1 = store.Chi2_Pion.at(1);
   P.Track_Chi2_Pion_Plane2 = store.Chi2_Pion.at(2);
   P.Track_Chi2_Pion_3Plane = store.Chi2_Pion_3Plane;
   P.Track_Chi2_Muon_Plane0 = store.Chi2_Muon.at(0);
   P.Track_Chi2_Muon_Plane1 = store.Chi2_Muon.at(1);
   P.Track_Chi2_Muon_Plane2 = store.Chi2_Muon.at(2);
   P.Track_Chi2_Muon_3Plane = store.Chi2_Muon_3Plane;

/*
   // LLR PID Calculation

   double this_llr_pid=0;
   double this_llr_pid_score=0;
   double this_llr_pid_kaon=0;
   double this_llr_pid_score_kaon=0;
   
   double this_llr_pid_kaon_partial = 0;
   double this_llr_pid_score_kaon_partial = 0;

   for(auto const &calo : caloFromTrack){

      auto const &plane = calo->PlaneID().Plane;
      auto const &dedx_values = calo->dEdx();
      auto const &rr = calo->ResidualRange();
      auto const &pitch = calo->TrkPitchVec();
      std::vector<std::vector<float>> par_values;
      par_values.push_back(rr);
      par_values.push_back(pitch);

      // Get parital length PIDs
      std::vector<std::vector<float>> par_values_partial;
      std::vector<float> dedx_values_partial,rr_partial,pitch_partial;      
      if(calo->dEdx().size() != calo->ResidualRange().size() || calo->ResidualRange().size() != calo->TrkPitchVec().size())
         throw cet::exception("SubModuleReco") << "Track calo point list size mismatch" << std::endl;
      for(size_t i_p=0;i_p<calo->dEdx().size();i_p++){
         if(rr.at(i_p) > ResRangeCutoff) continue;
         dedx_values_partial.push_back(calo->dEdx().at(i_p));
         rr_partial.push_back(calo->ResidualRange().at(i_p));
         pitch_partial.push_back(calo->TrkPitchVec().at(i_p));        
      }
      par_values_partial.push_back(rr_partial);
      par_values_partial.push_back(pitch_partial);
  
      if(calo->ResidualRange().size() == 0) continue;

      float calo_energy = 0;
      for(size_t i=0;i<dedx_values.size();i++)
         calo_energy += dedx_values[i] * pitch[i];

      float llr_pid = llr_pid_calculator.LLR_many_hits_one_plane(dedx_values,par_values,plane);
      float llr_pid_kaon = llr_pid_calculator_kaon.LLR_many_hits_one_plane(dedx_values,par_values,plane);
      this_llr_pid += llr_pid;
      this_llr_pid_kaon += llr_pid_kaon;

     // Partial length calculation
     float calo_energy_partial = 0;
      for(size_t i=0;i<dedx_values_partial.size();i++)
         calo_energy_partial += dedx_values_partial[i] * pitch_partial[i];

      float llr_pid_kaon_partial = llr_pid_calculator_kaon.LLR_many_hits_one_plane(dedx_values_partial,par_values_partial,plane);
      this_llr_pid_kaon_partial += llr_pid_kaon_partial;     
   }

   this_llr_pid_score = atan(this_llr_pid/100.)*2/3.14159266;
   this_llr_pid_score_kaon = atan(this_llr_pid_kaon/100.)*2/3.14159266;
   this_llr_pid_score_kaon_partial = atan(this_llr_pid_kaon_partial/100.)*2/3.14159266;

   P.Track_LLR_PID = this_llr_pid_score;
   P.Track_LLR_PID_Kaon = this_llr_pid_score_kaon;
   P.Track_LLR_PID_Kaon_Partial = this_llr_pid_score_kaon_partial;
   */


/*
   // LLR PID Scores Calculation
   LLRPID_Result LLPIDs = LLRPIDCalc.GetScores(caloFromTrack);
   P.Track_LLR_PID = LLPIDs.Score;
   P.Track_LLR_PID_Kaon = LLPIDs.Score_Kaon;
   P.Track_LLR_PID_Kaon_Partial = LLPIDs.Score_Kaon_Partial;

   // Mean dE/dX Calculation
   dEdXStore dEdXs = dEdXCalc.ThreePlaneMeandEdX(trk,caloFromTrack);
   P.MeandEdX_Plane0 = dEdXs.Plane0;
   P.MeandEdX_Plane1 = dEdXs.Plane1;
   P.MeandEdX_Plane2 = dEdXs.Plane2;
   P.MeandEdX_ThreePlane = dEdXs.ThreePlaneAverage;
*/


/*
   // 3 Plane Proton PID (Pip Hamilton)
   std::vector<art::Ptr<anab::ParticleID>> trackPID = Assoc_TrackPID->at(trk.key());

   std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();

   for(size_t i_algscore=0;i_algscore<AlgScoresVec.size();i_algscore++){

      anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);

      if(TMath::Abs(AlgScore.fAssumedPdg) == 2212 && AlgScore.fAlgName=="ThreePlaneProtonPID" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward)
         P.TrackPID = std::log(AlgScore.fValue);

      if(TMath::Abs(AlgScore.fAssumedPdg) == 321) std::cout << AlgScore.fAlgName << std::endl;

   }
*/
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::GetVertexData(const art::Ptr<recob::PFParticle> &pfp,RecoParticle &P){ // track vertex

   std::vector<art::Ptr<recob::Vertex>> pfpVertex = Assoc_PFParticleVertex->at(pfp.key());

   auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

   for(const art::Ptr<recob::Vertex> &vtx : pfpVertex){

      geo::Point_t point = {vtx->position().X(),vtx->position().Y(),vtx->position().Z()};                
      geo::Vector_t sce_corr = SCE->GetPosOffsets(point);

      TVector3 pos(vtx->position().X()+sce_corr.X(),vtx->position().Y()-sce_corr.Y(),vtx->position().Z()-sce_corr.Z());

      P.SetVertex(pos);
      P.Displacement = (pos-theData.RecoPrimaryVertex).Mag();

   }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

bool SubModuleReco::ApplyNuCCInclusiveFilter(art::Event const& e){

  // Check if event passed the NuCC inclusive filter
  bool IsNuCCInclusive = false;
  string process(!IsData ? "OverlayFiltersPostStage2" : "DataFiltersPostStage2");
  art::InputTag trigResInputTag("TriggerResults","",process.data()); // the last is the name of process where the filters were run
  art::ValidHandle<art::TriggerResults> trigRes = e.getValidHandle<art::TriggerResults>(trigResInputTag);
	
  fhicl::ParameterSet pset;
  if (!fhicl::ParameterSetRegistry::get(trigRes->parameterSetID(), pset)) { throw cet::exception("PSet Not Found???"); }
  std::vector<std::string> trigger_path_names = pset.get<std::vector<std::string> >("trigger_paths", {});
  if (trigger_path_names.size()!=trigRes->size()) { throw cet::exception("Size mismatch???"); }
  for (size_t itp=0;itp<trigRes->size();itp++) {
    //cout << "Filter name " << trigger_path_names.at(itp) << " decision=" << trigRes->at(itp).accept() << endl;
    if (trigger_path_names.at(itp)=="NuCC") {
      IsNuCCInclusive = trigRes->at(itp).accept();
    }
  }

  return IsNuCCInclusive;
	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void SubModuleReco::SetIndices(std::vector<bool> IsSignal, std::vector<bool> IsSignal_NuMuP, std::vector<bool> IsSignal_PiPPi0){
      
   bool ContainsSignal = std::find(IsSignal.begin(),IsSignal.end(),true) == IsSignal.end() 
                      || std::find(IsSignal_NuMuP.begin(),IsSignal_NuMuP.end(),true) == IsSignal_NuMuP.end()
                      || std::find(IsSignal_PiPPi0.begin(),IsSignal_PiPPi0.end(),true) == IsSignal_PiPPi0.end();

   bool found_muon=false, found_kaon=false, found_decaymuon=false, found_decaypion=false;
   bool found_kaon_as_shower=false, found_decaymuon_as_shower=false, found_decaypion_as_shower=false;

   for(size_t i_tr=0;i_tr<theData.TrackPrimaryDaughters.size();i_tr++){

      RecoParticle P = theData.TrackPrimaryDaughters.at(i_tr);

      if(!found_muon && abs(P.TrackTruePDG) == 13 && P.TrackTrueOrigin == 1){ 
         theData.TrueMuonIndex = P.Index;
         found_muon = true; 
      }

      if(!found_kaon && abs(P.TrackTruePDG) == 321 && P.TrackTrueOrigin == 1){ 
         theData.TrueKaonIndex = P.Index;
         found_kaon = true; 
      }

      if(ContainsSignal && !found_decaymuon && P.TrackTruePDG == -13 && (P.TrackTrueOrigin == 2 || P.TrackTrueOrigin == 3)){
         theData.TrueDecayMuonIndex = P.Index;
         found_decaymuon = true;
      }

      if(ContainsSignal && !found_decaypion && P.TrackTruePDG == 211 && (P.TrackTrueOrigin == 2 || P.TrackTrueOrigin == 3)){
         theData.TrueDecayPionIndex = P.Index;
         found_decaypion = true;
      }
   }

   for(size_t i_sh=0;i_sh<theData.ShowerPrimaryDaughters.size();i_sh++){

      RecoParticle P = theData.ShowerPrimaryDaughters.at(i_sh);

      /*
      if(!found_muon && abs(P.TrackTruePDG) == 13 && P.TrackTrueOrigin == 1){ 
         theData.TrueMuonIndex = P.Index;
         found_muon_as_shower = true; 
      }
      */
      if(!found_kaon && abs(P.TrackTruePDG) == 321 && P.TrackTrueOrigin == 1){ 
         theData.TrueKaonIndex = P.Index;
         found_kaon_as_shower = true; 
      }

      if(ContainsSignal && !found_decaymuon && P.TrackTruePDG == -13 && (P.TrackTrueOrigin == 2 || P.TrackTrueOrigin == 3)){
         theData.TrueDecayMuonIndex = P.Index;
         found_decaymuon_as_shower = true;
      }

      if(ContainsSignal && !found_decaypion && P.TrackTruePDG == 211 && (P.TrackTrueOrigin == 2 || P.TrackTrueOrigin == 3)){
         theData.TrueDecayPionIndex = P.Index;
         found_decaypion_as_shower = true;
      }
   }
  

if(ContainsSignal && found_muon && found_kaon && ( found_decaymuon || found_decaypion )) theData.GoodReco = true;
if(ContainsSignal && found_muon && found_kaon && found_decaymuon) theData.GoodReco_NuMuP = true;
if(ContainsSignal && found_muon && found_kaon && found_decaypion) theData.GoodReco_PiPPi0 = true;
if(ContainsSignal && found_muon && found_kaon && ( found_decaymuon_as_shower || found_decaypion_as_shower )) theData.GoodPrimaryReco = true;
if(ContainsSignal && found_muon && found_kaon_as_shower && ( found_decaymuon_as_shower || found_decaypion_as_shower )) theData.GoodRecoAsShower = true;
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
