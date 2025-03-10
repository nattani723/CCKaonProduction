#ifndef _SubModuleGeneratorTruth_cxx_
#define _SubModuleGeneratorTruth_cxx_

#include "SubModuleGeneratorTruth.h"

using namespace cckaon;

/////////////////////////////////////////////////////////////////////////////////////////////////////////

SubModuleGeneratorTruth::SubModuleGeneratorTruth(art::Event const& e,fhicl::ParameterSet pset,bool particlegunmode) :
ParticleGunMode(particlegunmode)
{

   if(!e.getByLabel(pset.get<std::string>("GeneratorModuleLabel","generator"),Handle_MCTruth))  
      throw cet::exception("SubModuleGeneratorTruth") << "No MC Truth data product!" << std::endl;

   art::fill_ptr_vector(Vect_MCTruth,Handle_MCTruth);  

   if(!e.getByLabel(pset.get<std::string>("FluxModuleLabel","generator"),Handle_MCFlux))  
      throw cet::exception("SubModuleGeneratorTruth") << "No MC Truth data product!" << std::endl;

   art::fill_ptr_vector(Vect_MCFlux,Handle_MCFlux);  

   HyperonPDGs = pset.get<std::vector<int>>("HyperonPDGs",{3122,3212,3112,3222});
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

GeneratorTruth SubModuleGeneratorTruth::GetGeneratorTruth(){

   if(!Vect_MCTruth.size()){
      std::cout << "MCTruth vector is empty" << std::endl;
      return theTruth;
   }

   theTruth.NMCTruths = Vect_MCTruth.size();

   int i_truth=0;

   for(const art::Ptr<simb::MCFlux> &theMCFlux : Vect_MCFlux){ 
     theTruth.TrueDecayPosition_X.push_back(theMCFlux->fvx);
     theTruth.TrueDecayPosition_Y.push_back(theMCFlux->fvy);
     theTruth.TrueDecayPosition_Z.push_back(theMCFlux->fvz);
     theTruth.NuL.push_back(theMCFlux->fdk2gen + theMCFlux->fgen2vtx);
   }

   for(const art::Ptr<simb::MCTruth> &theMCTruth : Vect_MCTruth){

      simb::MCNeutrino Nu = theMCTruth->GetNeutrino();

      int mode = Nu.Mode();
      int ccnc = Nu.CCNC();

      if(ccnc == 0) theTruth.CCNC.push_back("CC");
      else theTruth.CCNC.push_back("NC");

      theTruth.NuPDG.push_back(Nu.Nu().PdgCode());
      theTruth.NuE.push_back(Nu.Nu().Trajectory().E(0));

      if(mode == 0) theTruth.Mode.push_back("QEL");
      else if(mode == 1) theTruth.Mode.push_back("RES");
      else if(mode == 2) theTruth.Mode.push_back("DIS");
      else if(mode == 3) theTruth.Mode.push_back("COH");
      else if(mode == 5) theTruth.Mode.push_back("ElectronScattering");
      else if(mode == 10) theTruth.Mode.push_back("MEC");
      else if(mode == 11) theTruth.Mode.push_back("Diffractive");
      else if(mode == 1095) theTruth.Mode.push_back("HYP");
      else theTruth.Mode.push_back("Other");	

      for(int k_particles=0;k_particles<theMCTruth->NParticles();k_particles++){

         simb::MCParticle Part = theMCTruth->GetParticle(k_particles);

         //if((isLepton(Part.PdgCode()) || isNeutrino(Part.PdgCode())) && Part.StatusCode() == 1) 
         // theTruth.TruePrimaryVertex.SetXYZ(Part.Vx(),Part.Vy(),Part.Vz());

         if(isNeutrino(Part.PdgCode()) && Part.StatusCode() == 0){
            SimParticle P = MakeSimParticle(Part);
            P.Origin = 0;
            P.MCTruthIndex = i_truth;
            theTruth.Neutrino.push_back(P);
         }

         // If there is a kaon in the final state in a RES/DIS event, change mode to KAON
         if(isKaon(Part.PdgCode()) && Part.StatusCode() == 1 && ( mode == 1 || mode == 2) ) theTruth.Mode.back() = "KAON";
         //if(isKaonP(Part.PdgCode()) && Part.StatusCode() == 1 && ( mode == 1 || mode == 2) ) theTruth.Mode.back() = "KAONP";

         if(Part.StatusCode() == 1 && Part.PdgCode() == 2212) theTruth.EventHasFinalStateProton = true;
         if(Part.StatusCode() == 1 && (Part.PdgCode() == 211 || Part.PdgCode() == -211)) theTruth.EventHasFinalStatePion = true;
         if(Part.StatusCode() == 1 && isHyperon(Part.PdgCode()) && std::find(HyperonPDGs.begin(),HyperonPDGs.end(),abs(Part.PdgCode())) != HyperonPDGs.end()){
            theTruth.EventHasHyperon = true;
         }
         if(Part.StatusCode() == 1 && isKaon(Part.PdgCode())) theTruth.EventHasKaon = true;        
         if(Part.StatusCode() == 1 && isKaonP(Part.PdgCode())) theTruth.EventHasKaonP = true;        
         if(Part.StatusCode() == 1 && isKaonM(Part.PdgCode())) theTruth.EventHasKaonM = true;        
         if(Part.StatusCode() == 1 && isKaon0(Part.PdgCode())) theTruth.EventHasKaon0 = true;        
 
         if((isLepton(Part.PdgCode()) || isNeutrino(Part.PdgCode())) && Part.StatusCode() == 1) {
            theTruth.TruePrimaryVertex_X.push_back(Part.Vx());
            theTruth.TruePrimaryVertex_Y.push_back(Part.Vy());
            theTruth.TruePrimaryVertex_Z.push_back(Part.Vz());
            if(inActiveTPC(TVector3(Part.Vx(),Part.Vy(),Part.Vz()))) theTruth.NMCTruthsInTPC++;
         }

      }

      i_truth++;
   }

   if(!ParticleGunMode && theTruth.Neutrino.size() != Vect_MCTruth.size())         
      throw cet::exception("SubModuleGeneratorTruth") << "Sim Neutrino/MCTruth vector size mismatch" << std::endl;

   return theTruth;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
