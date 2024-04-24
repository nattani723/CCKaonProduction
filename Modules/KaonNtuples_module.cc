////////////////////////////////////////////////////////////////////////
// Class:       KaonNtuples
// Plugin Type: analyzer (art v3_03_01)
// File:        KaonNtuples_module.cc
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
#include "ubana/CCKaonProduction/Objects/SimParticle.h"
#include "ubana/CCKaonProduction/Objects/RecoParticle.h"
#include "ubana/CCKaonProduction/Objects/Helpers.h"

//algorithms
#include "ubana/CCKaonProduction/Alg/ConnectednessHelper.h"

//submodules
#include "ubana/CCKaonProduction/Modules/SubModules/SubModuleGeneratorTruth.h"
#include "ubana/CCKaonProduction/Modules/SubModules/SubModuleG4Truth.h"
#include "ubana/CCKaonProduction/Modules/SubModules/SubModuleReco.h"

namespace cckaon {
   class KaonNtuples;
}


class kaon::KaonNtuples : public art::EDAnalyzer {
public:
  explicit KaonNtuples(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  
  // Plugins should not be copied or assigned.
  KaonNtuples(KaonNtuples const&) = delete;
  KaonNtuples(KaonNtuples&&) = delete;
  KaonNtuples& operator=(KaonNtuples const&) = delete;
  KaonNtuples& operator=(KaonNtuples&&) = delete;
  
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
  TTree * MetaTree;
  
  // Basic event info
  unsigned int t_EventID;
  int t_run,t_subrun,t_event;
  
  double t_Weight=1.0;
  
  // Generator/Geant4 truth info
  
  int t_NMCTruths=0;	
  int t_NMCTruthsInTPC=0;	
  
  std::vector<std::string> t_Mode;
  std::vector<std::string> t_CCNC;
  
  // Flags applying to the entire event
  bool t_EventHasKaonPScatter;
  bool t_EventHasHyperon;
  bool t_EventHasKaon;
  bool t_EventHasKaonP;
  bool t_EventHasKaonP_NuMuP;
  bool t_EventHasKaonP_PiPPi0;
  bool t_EventHasKaonM;
  bool t_EventHasKaon0;
  // bool t_EventHasK0S;
  bool t_GoodReco;
  bool t_GoodPrimaryReco;
  bool t_GoodRecoAsShower;
  
  // Flags applying to each MCTruth
  std::vector<bool> t_InActiveTPC;
  std::vector<bool> t_IsHyperon;
  std::vector<bool> t_IsKaon;
  std::vector<bool> t_IsKaonP;
  std::vector<bool> t_IsKaonP_NuMuP;
  std::vector<bool> t_IsKaonP_PiPPi0;
  std::vector<bool> t_IsKaonP_2PiPPiM;
  std::vector<bool> t_IsKaonP_ENuE;
  std::vector<bool> t_IsKaonP_2PiNPiP;
  std::vector<bool> t_IsKaonP_Others;
  std::vector<bool> t_IsKaonM;
  std::vector<bool> t_IsKaon0;
  std::vector<bool> t_IsAssociatedKaonP;
  std::vector<bool> t_IsSignal;
  std::vector<bool> t_IsSignal_NuMuP;
  std::vector<bool> t_IsSignal_PiPPi0;
  /*
      std::vector<bool> t_InActiveTPC;
      std::vector<bool> t_IsHyperon;
      std::vector<bool> t_IsLambda;
      std::vector<bool> t_IsLambdaCharged;
      std::vector<bool> t_IsSigmaZero; 		
      std::vector<bool> t_IsSigmaZeroCharged; 		
      std::vector<bool> t_IsAssociatedHyperon;
      std::vector<bool> t_IsSignal;
      std::vector<bool> t_IsSignalSigmaZero;
      std::vector<bool> t_IsKaon;
      std::vector<bool> t_IsK0S;
      std::vector<bool> t_IsK0SCharged;
  */

  bool t_EventHasFinalStateNeutron;

  std::vector<SimParticle> t_Neutrino;
  std::vector<SimParticle> t_Lepton;
  std::vector<SimParticle> t_PrimaryHyperon;
  std::vector<SimParticle> t_PrimaryNucleon;
  std::vector<SimParticle> t_PrimaryPion;
  std::vector<SimParticle> t_PrimaryKaon; 
  std::vector<SimParticle> t_PrimaryKaonP; 
  std::vector<SimParticle> t_PrimaryKaonM; 
  std::vector<SimParticle> t_PrimaryNucleus; 
  std::vector<SimParticle> t_HyperonDecay;
  std::vector<SimParticle> t_KaonPDecay;
  std::vector<SimParticle> t_KaonPDecay_NuMuP; 
  std::vector<SimParticle> t_KaonPDecay_PiPPi0;
  std::vector<SimParticle> t_KaonMDecay;
  std::vector<SimParticle> t_Kaon0Decay;
  std::vector<SimParticle> t_NeutralKaonDecayK0SL;
  
  std::vector<double> t_TruePrimaryVertex_X;
  std::vector<double> t_TruePrimaryVertex_Y;
  std::vector<double> t_TruePrimaryVertex_Z;
  //TVector3 t_TruePrimaryVertex;
  //TVector3 t_DecayVertex;
  std::vector<double> t_DecayVertex_X;
  std::vector<double> t_DecayVertex_Y;
  std::vector<double> t_DecayVertex_Z;

  int t_NPrimaryDaughters;
  int t_NPrimaryTrackDaughters;
  int t_NPrimaryShowerDaughters;
  int t_NOtherTracks;
  int t_NOtherRebuiltTracks;
  int t_NOtherShowers;

  //std::vector<RecoParticle> t_TrackPrimary;
  std::vector<RecoParticle> t_TrackPrimaryDaughters;
  std::vector<RecoParticle> t_ShowerPrimaryDaughters;
  std::vector<RecoParticle> t_TrackOthers;
  std::vector<RecoParticle> t_TrackRebuiltOthers;
  std::vector<RecoParticle> t_ShowerOthers;
  
  TVector3 t_RecoPrimaryVertex;
  
  std::vector<RecoParticle> t_RepassTrackPrimaryDaughters;
  std::vector<RecoParticle> t_RepassShowerPrimaryDaughters;   
  
  //std::vector<std::string> t_SysDials;
  //std::vector<std::vector<std::vector<double>>> t_SysWeights;
  
  std::vector<std::string> t_SysDials;
  std::vector<std::vector<double>> t_SysWeights;
  
  /////////////////////////
  // Metadata for sample //
  /////////////////////////
  
  int m_NEvents;
  int m_NHyperons;
  int m_NKaons;
  int m_NKaonPs;
  int m_NKaonPs_NuMuP;
  int m_NKaonPs_PiPPi0;
  int m_NSignal;
  int m_NSignal_NuMuP;  
  int m_NSignal_PiPPi0;  
  int m_NGoodReco;

  double m_POT = 0; //total POT of the sample

  //////////////////////////
  //   FHICL PARAMETERS   //
  //////////////////////////
  
  bool f_GetGeneratorInfo;
  bool f_GetG4Info;
  bool f_GetRecoInfo;
  bool f_GetConnInfo;
  
  fhicl::ParameterSet f_Generator;
  fhicl::ParameterSet f_G4;
  fhicl::ParameterSet f_Reco;
  std::vector<fhicl::ParameterSet> f_RecoRepass;
  std::string f_WireLabel;
  std::vector<art::InputTag> f_WeightLabels;
  std::string f_POTSummaryLabel;
  
  bool f_ParticleGun = false;
  bool f_IsData;
  bool f_Debug = false;
  
};

////////////////////////////////////////////////////
// Setup module labels/read in fhicl settings     //
////////////////////////////////////////////////////

cckaon::KaonNtuples::KaonNtuples(fhicl::ParameterSet const& p)
   : EDAnalyzer{p},
   f_GetGeneratorInfo(p.get<bool>("GetGeneratorInfo",true)),   
   f_GetG4Info(p.get<bool>("GetG4Info",true)),   
   f_GetRecoInfo(p.get<bool>("GetRecoInfo",true)),   
   //f_GetConnInfo(p.get<bool>("GetConnInfo",true)),   
   f_Generator(p.get<fhicl::ParameterSet>("Generator")),
   f_G4(p.get<fhicl::ParameterSet>("Geant4")),
   f_Reco(p.get<fhicl::ParameterSet>("Reco")),
   f_RecoRepass(p.get<std::vector<fhicl::ParameterSet>>("RecoRepass",{})),
   f_WireLabel(p.get<std::string>("WireLabel")),
   f_WeightLabels(p.get<std::vector<art::InputTag>>("WeightCalculators",{})),
   f_POTSummaryLabel(p.get<std::string>("POTSummaryLabel")),
   f_ParticleGun(p.get<bool>("ParticleGun",false)),
   f_IsData(p.get<bool>("IsData")),
   f_Debug(p.get<bool>("Debug",false)),
   Conn_Helper(p.get<bool>("DrawConnectedness",false))   // ,
{
   if(f_WeightLabels.size()){
      std::cout << "Getting weights from data products with tags:" << std::endl;
      for(size_t i=0;i<f_WeightLabels.size();i++) std::cout << f_WeightLabels.at(i) << std::endl;
   }

   /*
   if(f_Reco.get<bool>("IncludeCosmics",false) && f_GetConnInfo)
      std::cout << std::endl << "HyperonNTuples WARNING: Requesting connectedness information with cosmics included. This will take a very long time." << std::endl
                << "Set GetConnInfo fhicl parameter to false disable connectedness" << std::endl << std::endl;
   */
      
}

void cckaon::KaonNtuples::analyze(art::Event const& e)
{
   if(f_Debug) std::cout << "New Event" << std::endl;

   //begin by resetting everything

   t_Weight = 1.0;
   //t_Mode = "NONE";
   t_Mode.clear();
   t_CCNC.clear();
   t_NMCTruths = 0;
   t_NMCTruthsInTPC = 0;
   t_EventHasFinalStateNeutron = false;

   t_InActiveTPC.clear();
   t_IsHyperon.clear();
   t_IsKaon.clear();
   t_IsKaonP.clear();
   t_IsKaonP_NuMuP.clear();
   t_IsKaonP_PiPPi0.clear();
   t_IsKaonP_ENuE.clear();
   t_IsKaonP_2PiNPiP.clear();
   t_IsKaonP_2PiPPiM.clear();
   t_IsKaonP_Others.clear();
   t_IsAssociatedKaonP.clear();
   t_IsKaonM.clear();
   t_IsKaon0.clear();
   t_IsSignal.clear();
   t_IsSignal_NuMuP.clear();
   t_IsSignal_PiPPi0.clear();
   t_EventHasKaonPScatter = false;
   t_EventHasHyperon = false;
   t_EventHasKaon = false;
   t_EventHasKaonP = false;
   t_EventHasKaonP_NuMuP = false;
   t_EventHasKaonP_PiPPi0 = false;
   t_EventHasKaonM = false;
   t_EventHasKaon0 = false;
   t_GoodReco = false;

   t_Neutrino.clear();
   t_Lepton.clear();
   t_PrimaryHyperon.clear();
   t_PrimaryNucleon.clear();
   t_PrimaryPion.clear();
   t_PrimaryKaon.clear();
   t_PrimaryKaonP.clear();
   t_PrimaryKaonM.clear();
   t_PrimaryNucleus.clear();
   t_HyperonDecay.clear();
   t_KaonPDecay.clear();
   t_KaonPDecay_NuMuP.clear();
   t_KaonPDecay_PiPPi0.clear();
   t_KaonMDecay.clear();
   t_Kaon0Decay.clear();
   t_NeutralKaonDecayK0SL.clear();

   t_TruePrimaryVertex_X.clear();
   t_TruePrimaryVertex_Y.clear();
   t_TruePrimaryVertex_Z.clear();
   t_DecayVertex_X.clear();
   t_DecayVertex_Y.clear();
   t_DecayVertex_Z.clear();

   t_NPrimaryDaughters = 0; //number of primary daughters
   t_NPrimaryTrackDaughters=0; //num of track like primary daughters
   t_NPrimaryShowerDaughters=0; //num of shower like primary daughters
   t_NOtherRebuiltTracks=0;
   t_NOtherTracks=0;
   t_NOtherShowers=0;

   t_TrackPrimaryDaughters.clear();
   t_ShowerPrimaryDaughters.clear();
   t_TrackOthers.clear();
   t_TrackRebuiltOthers.clear();
   t_ShowerOthers.clear();

   t_RecoPrimaryVertex.SetXYZ(-1000,-1000,-1000); //position of reco'd primary vertex

   t_RepassTrackPrimaryDaughters.clear();
   t_RepassShowerPrimaryDaughters.clear();

   t_SysDials.clear();
   t_SysWeights.clear();

   // General Event Info

   t_EventID = e.id().event();
   t_run = e.run();
   t_subrun = e.subRun();
   t_event = e.event();

   // Event Generator Info

   if(!f_IsData && f_GetGeneratorInfo){

      if(f_Debug) std::cout << "Getting Event Generator Info" << std::endl;

      SubModuleGeneratorTruth* Generator_SM = new SubModuleGeneratorTruth(e,f_Generator,f_ParticleGun);
      GeneratorTruth GenT = Generator_SM->GetGeneratorTruth();

      t_Weight *= GenT.Weight;
      t_CCNC = GenT.CCNC;
      t_Mode = GenT.Mode;
      t_NMCTruths = GenT.NMCTruths;
      t_NMCTruthsInTPC = GenT.NMCTruthsInTPC;
      t_Neutrino = GenT.Neutrino;
      //t_TruePrimaryVertex_X = GenT.TruePrimaryVertex_X;
      //t_TruePrimaryVertex_Y = GenT.TruePrimaryVertex_Y;
      //t_TruePrimaryVertex_Z = GenT.TruePrimaryVertex_Z;

      t_EventHasFinalStateNeutron = GenT.EventHasFinalStateNeutron;

      delete Generator_SM;
   }

   // G4 Info

   if(!f_IsData && f_GetG4Info){

      if(f_Debug) std::cout << "Getting G4 Info" << std::endl;

      SubModuleG4Truth* G4_SM = new SubModuleG4Truth(e,f_G4,f_ParticleGun);
      G4Truth G4T = G4_SM->GetG4Info();

      t_InActiveTPC = G4T.InActiveTPC; // add more fv cuts
      t_IsHyperon = G4T.IsHyperon;
      t_IsKaon = G4T.IsKaon;
      t_IsKaonP = G4T.IsKaonP;
      t_IsKaonP_NuMuP = G4T.IsKaonP_NuMuP;
      t_IsKaonP_PiPPi0 = G4T.IsKaonP_PiPPi0;
      t_IsKaonP_2PiPPiM = G4T.IsKaonP_2PiPPiM;
      t_IsKaonP_ENuE = G4T.IsKaonP_ENuE;
      t_IsKaonP_2PiNPiP = G4T.IsKaonP_2PiNPiP;
      t_IsKaonP_Others = G4T.IsKaonP_Others;
      t_IsAssociatedKaonP = G4T.IsAssociatedKaonP;
      t_IsKaonM = G4T.IsKaonM;
      t_IsKaon0 = G4T.IsKaon0;

      t_EventHasKaonPScatter = G4T.EventHasKaonPScatter;
      t_EventHasHyperon = G4T.EventHasHyperon;
      t_EventHasKaon = G4T.EventHasKaon;
      t_EventHasKaonP = G4T.EventHasKaonP;
      t_EventHasKaonP_NuMuP = G4T.EventHasKaonP;
      t_EventHasKaonP_PiPPi0 = G4T.EventHasKaonP;
      t_EventHasKaonM = G4T.EventHasKaonM;
      t_EventHasKaon0 = G4T.EventHasKaon0;

      t_Weight *= G4T.Weight;

      t_Lepton = G4T.Lepton;
      t_PrimaryHyperon = G4T.Hyperon;
      t_PrimaryNucleon = G4T.PrimaryNucleon;
      t_PrimaryPion = G4T.PrimaryPion;
      t_PrimaryKaon = G4T.PrimaryKaon;
      t_PrimaryKaonP = G4T.PrimaryKaonP;
      t_PrimaryKaonP = G4T.PrimaryKaonM;
      t_PrimaryNucleus = G4T.PrimaryNucleus;
      t_HyperonDecay = G4T.Decay;
      t_KaonPDecay = G4T.KaonPDecay;
      t_KaonPDecay_NuMuP = G4T.KaonPDecay_NuMuP;
      t_KaonPDecay_PiPPi0 = G4T.KaonPDecay_PiPPi0;
      t_KaonMDecay = G4T.KaonMDecay;
      t_Kaon0Decay = G4T.Kaon0Decay;
      t_NeutralKaonDecayK0SL = G4T.NeutralKaonDecayK0SL;

      t_TruePrimaryVertex_X = G4T.TruePrimaryVertex_X;
      t_TruePrimaryVertex_Y = G4T.TruePrimaryVertex_Y;
      t_TruePrimaryVertex_Z = G4T.TruePrimaryVertex_Z;

      t_DecayVertex_X = G4T.DecayVertex_X;
      t_DecayVertex_Y = G4T.DecayVertex_Y;
      t_DecayVertex_Z = G4T.DecayVertex_Z;

      t_IsSignal.resize(t_NMCTruths);
      t_IsSignal_NuMuP.resize(t_NMCTruths);
      t_IsSignal_PiPPi0.resize(t_NMCTruths);


      for(int i_t=0;i_t<t_NMCTruths;i_t++){
         t_IsSignal[i_t] = t_Mode.at(i_t) == "KAON" && t_InActiveTPC.at(i_t) && t_Neutrino.at(i_t).PDG == -14 && t_IsKaonP.at(i_t);
         t_IsSignal_NuMuP[i_t] = t_Mode.at(i_t) == "KAON" && t_InActiveTPC.at(i_t) && t_Neutrino.at(i_t).PDG == -14 && t_IsKaonP_NuMuP.at(i_t);
         t_IsSignal_PiPPi0[i_t] = t_Mode.at(i_t) == "KAON" && t_InActiveTPC.at(i_t) && t_Neutrino.at(i_t).PDG == -14 && t_IsKaonP_PiPPi0.at(i_t);
      }

      delete G4_SM;
   }

   // Reconstructed Info

   if(f_GetRecoInfo){

      if(f_Debug) std::cout << "Getting Reconstructed Info" << std::endl;

      SubModuleReco* Reco_SM = new SubModuleReco(e,f_IsData,f_Reco,f_ParticleGun);
      Reco_SM->PrepareInfo();
      Reco_SM->SetIndices(t_IsSignal,t_IsSignalSigmaZero);
      RecoData RecoD =  Reco_SM->GetInfo();   

      t_NPrimaryDaughters = RecoD.NPrimaryDaughters;
      t_NPrimaryTrackDaughters = RecoD.NPrimaryTrackDaughters;
      t_NPrimaryShowerDaughters = RecoD.NPrimaryShowerDaughters;
      t_NOtherTracks = RecoD.NOtherTracks;
      t_NOtherRebuiltTracks = RecoD.NOtherTracks;
      t_NOtherShowers = RecoD.NOtherShowers;
  
      t_TrackPrimaryDaughters = RecoD.TrackPrimaryDaughters;
      t_ShowerPrimaryDaughters = RecoD.ShowerPrimaryDaughters;
      t_TrackOthers = RecoD.TrackOthers;
      t_TrackRebuiltOthers = RecoD.TrackRebuiltOthers;
      t_ShowerOthers = RecoD.ShowerOthers;
      t_RecoPrimaryVertex = RecoD.RecoPrimaryVertex;
      t_GoodReco = RecoD.GoodReco;

      delete Reco_SM;

      // If configured to get repass of reconstruction
      if(f_RecoRepass.size()){
        std::cout << "Getting repassed information" << std::endl;
         for(size_t i_r=0;i_r<f_RecoRepass.size();i_r++){
            SubModuleReco* Reco_SM_Repass = new SubModuleReco(e,f_IsData,f_RecoRepass.at(i_r));
            Reco_SM_Repass->PrepareInfo();
            //Reco_SM->SetIndices(t_IsSignal,t_IsSignalSigmaZero);
            RecoData RecoD_Repass =  Reco_SM_Repass->GetInfo();   
            t_RepassTrackPrimaryDaughters.insert(t_RepassTrackPrimaryDaughters.end(),RecoD_Repass.TrackPrimaryDaughters.begin(),RecoD_Repass.TrackPrimaryDaughters.end());
            t_RepassShowerPrimaryDaughters.insert(t_RepassShowerPrimaryDaughters.end(),RecoD_Repass.ShowerPrimaryDaughters.begin(),RecoD_Repass.ShowerPrimaryDaughters.end());
            delete Reco_SM_Repass;
         }
      }
   }


   // Systematics weights if requested
/*
   if(!f_IsData){

      t_SysWeights.resize(t_NMCTruths);

      for(size_t i_w=0;i_w<f_WeightLabels.size();i_w++){

         std::cout << "Getting new weight products with label " << f_WeightLabels.at(i_w) << std::endl;

         // Try to get some systematics info
         art::Handle<std::vector<evwgh::MCEventWeight>> Handle_EventWeight;
         std::vector<art::Ptr<evwgh::MCEventWeight>> Vect_EventWeight;

         if(!e.getByLabel(f_WeightLabels.at(i_w),Handle_EventWeight)) 
            throw cet::exception("KaonNtuples") << "No EventWeight Found!" << std::endl;

         art::fill_ptr_vector(Vect_EventWeight,Handle_EventWeight);

         if(!Vect_EventWeight.size())
            throw cet::exception("KaonNtuples") << "Weight vector empty!" << std::endl;

         if(Vect_EventWeight.size() != (size_t)t_NMCTruths)
            throw cet::exception("KaonNtuples") << "Weight vector size != NMCTruths" << std::endl;

         for(size_t i_tr=0;i_tr<Vect_EventWeight.size();i_tr++){       

            std::cout << "Getting weights for truth " << i_tr << std::endl;

            std::map<std::string,std::vector<double>> theWeights = Vect_EventWeight.at(i_tr)->fWeight;
            std::map<std::string,std::vector<double>>::iterator it;

            for(it = theWeights.begin(); it != theWeights.end();it++){

               if(it->first ==  "empty") continue;

               bool dial_found=false;

               // Search the dial vector for this dial         
               for(size_t i_d=0;i_d<t_SysDials.size();i_d++){
                  if(it->first == t_SysDials.at(i_d)){
                     t_SysWeights.at(i_tr).at(i_d).insert(t_SysWeights.at(i_tr).at(i_d).end(),it->second.begin(),it->second.end());
                     dial_found=true;
                  }
               }  

               if(!dial_found){
                  if(i_tr != 0)
                     throw cet::exception("KaonNtuples") << "Malformed systematics weight vectors" << std::endl;
                  t_SysDials.push_back(it->first);
                  t_SysWeights.at(0).push_back(it->second);
                  for(size_t i_tr2=0;i_tr2<t_SysWeights.size();i_tr2++)
                     t_SysWeights.at(i_tr2).resize(t_SysWeights.at(0).size());
               } // if new dial

            } // iterate over weight products       
         } // i_tr

      } // i_w
   } // if not data
*/

   if(!f_IsData){

      std::vector<std::map<std::string,std::vector<double>>> theweightmap(t_NMCTruths); 

      for(size_t i_w=0;i_w<f_WeightLabels.size();i_w++){

         std::cout << "Getting new weight products with label " << f_WeightLabels.at(i_w) << std::endl;

         art::Handle<std::vector<evwgh::MCEventWeight>> Handle_EventWeight;
         std::vector<art::Ptr<evwgh::MCEventWeight>> Vect_EventWeight;

         if(!e.getByLabel(f_WeightLabels.at(i_w),Handle_EventWeight)) 
            throw cet::exception("KaonNtuples") << "No EventWeight Found!" << std::endl;

         art::fill_ptr_vector(Vect_EventWeight,Handle_EventWeight);

         if(!Vect_EventWeight.size())
            throw cet::exception("KaonNtuples") << "Weight vector empty!" << std::endl;

         if(Vect_EventWeight.size() != (size_t)t_NMCTruths)
            throw cet::exception("KaonNtuples") << "Weight vector size != NMCTruths" << std::endl;

         for(size_t i_tr=0;i_tr<Vect_EventWeight.size();i_tr++){       

            std::cout << "Getting weights for truth " << i_tr << std::endl;

            std::map<std::string,std::vector<double>> theWeights = Vect_EventWeight.at(i_tr)->fWeight;
            std::map<std::string,std::vector<double>>::iterator it;

            for(it = theWeights.begin();it != theWeights.end();it++){

               if(it->first ==  "empty") continue;

               bool dial_found=false;

               std::map<std::string,std::vector<double>>::iterator it2;
               for(it2 = theweightmap.at(i_tr).begin();it2 != theweightmap.at(i_tr).end();it2++){
                  if(it->first == it2->first){
                     dial_found = true;
                     theweightmap.at(i_tr)[it->first].insert(theweightmap.at(i_tr)[it->first].end(),it->second.begin(),it->second.end());
                  }
               }

               if(!dial_found)
                  theweightmap.at(i_tr)[it->first] = it->second;

            }
         } // i_tr
      }

      // Organise the weights
      if(theweightmap.size()){
         std::map<std::string,std::vector<double>>::iterator it;
         for(it = theweightmap.at(0).begin();it != theweightmap.at(0).end();it++){
            std::cout << "Organising weights for dial " << it->first << std::endl;

            t_SysDials.push_back(it->first);
            t_SysWeights.push_back(it->second);

            for(size_t i_tr=1;i_tr<theweightmap.size();i_tr++){
               if(theweightmap.at(i_tr).find(it->first) == theweightmap.at(i_tr).end()) 
                  throw cet::exception("KaonNtuples") << "Dial " << it->first << " not found in weights for MC truth " << i_tr << std::endl;
               if(theweightmap.at(i_tr)[it->first].size() != t_SysWeights.back().size())
                  throw cet::exception("KaonNtuples") << "Dial " << it->first << " weight vector mismatch" << std::endl;                        
               for(size_t i_w=0;i_w<t_SysWeights.back().size();i_w++)
                  t_SysWeights.back().at(i_w) *= theweightmap.at(i_tr)[it->first].at(i_w);
            }                  
         }
      }
   }// if(!f_IsData)

   FinishEvent();
}

///////////////////////////////////////////////////////////////	
// Finished processing event - update Metadata and fill tree //
///////////////////////////////////////////////////////////////

void cckaon::KaonNtuples::FinishEvent(){

   if(f_Debug) std::cout << "Finishing Event" << std::endl;

   OutputTree->Fill();

   m_NEvents++;

   if(std::find(t_IsHyperon.begin(), t_IsHyperon.end(), true) != t_IsHyperon.end()) m_NHyperons++;
   if(std::find(t_IsKaon.begin(), t_IsKaon.end(), true) != t_IsKaon.end()) m_NKaons++;
   if(std::find(t_IsKaonP.begin(), t_IsKaonP.end(), true) != t_IsKaonP.end()) m_NKaons++;
   if(std::find(t_IsKaonP_NuMuP.begin(), t_IsKaonP_NuMuP.end(), true) != t_IsKaonP_NuMuP.end()) m_NKaons++;
   if(std::find(t_IsKaonP_PiPPi0.begin(), t_IsKaonP_PiPPi0.end(), true) != t_IsKaonP_PiPPi0.end()) m_NKaons++;
   if(std::find(t_IsSignal.begin(), t_IsSignal.end(), true) != t_IsSignal.end()) m_NSignal++;
   if(std::find(t_IsSignal_NuMuP.begin(), t_IsSignal_NuMuP.end(), true) != t_IsSignal_NuMuP.end()) m_NSignal_NuMuP++;
   if(std::find(t_IsSignal_PiPPi0.begin(), t_IsSignal_PiPPi0.end(), true) != t_IsSignal_PiPPi0.end()) m_NSignal_PiPPi0++;
   if(t_GoodReco) m_NGoodReco++;

   if(f_Debug) std::cout << "Finished event" << std::endl;

}

///////////////////////////////////////////////////////////////	

void cckaon::KaonNtuples::beginJob(){

   if(f_Debug) std::cout << "Begin job" << std::endl;

   art::ServiceHandle<art::TFileService> tfs;

   //////////////////////////////////////////
   //             Output Tree	           //
   //////////////////////////////////////////

   OutputTree=tfs->make<TTree>("OutputTree","Truth Info Tree");
   
   OutputTree->Branch("IsData",&f_IsData);
   OutputTree->Branch("EventID",&t_EventID);
   OutputTree->Branch("run",&t_run);
   OutputTree->Branch("subrun",&t_subrun);
   OutputTree->Branch("event",&t_event);

   OutputTree->Branch("Weight",&t_Weight);
   OutputTree->Branch("Mode","vector<string>",&t_Mode);
   OutputTree->Branch("CCNC","vector<string>",&t_CCNC);
   OutputTree->Branch("NMCTruths",&t_NMCTruths);
   OutputTree->Branch("NMCTruthsInTPC",&t_NMCTruthsInTPC);
   OutputTree->Branch("InActiveTPC","vector<bool>",&t_InActiveTPC);
   OutputTree->Branch("IsHyperon","vector<bool>",&t_IsHyperon);
   OutputTree->Branch("IsKaon","vector<bool>",&t_IsKaon);
   OutputTree->Branch("IsKaonP","vector<bool>",&t_IsKaonP);
   OutputTree->Branch("IsKaonP_NuMuP","vector<bool>",&t_IsKaonP_NuMuP);
   OutputTree->Branch("IsKaonP_PiPPi0","vector<bool>",&t_IsKaonP_PiPPi0);
   OutputTree->Branch("IsKaonM","vector<bool>",&t_IsKaonM);
   OutputTree->Branch("IsKaon0","vector<bool>",&t_IsKaon0);
   OutputTree->Branch("IsSignal","vector<bool>",&t_IsSignal); 
   OutputTree->Branch("IsSignal_NuMuP","vector<bool>",&t_IsSignal_NuMuP); 
   OutputTree->Branch("IsSignal_PiPPi0","vector<bool>",&t_IsSignal_PiPPi0); 
   OutputTree->Branch("GoodReco",&t_GoodReco);

   OutputTree->Branch("EventHasNeutronScatter",&t_EventHasNeutronScatter);
   OutputTree->Branch("EventHasHyperon",&t_EventHasHyperon);
   OutputTree->Branch("EventHasKaonP",&t_EventHasKaonP);
   OutputTree->Branch("EventHasKaonP_NuMuP",&t_EventHasKaonP_NuMuP);
   OutputTree->Branch("EventHasKaonP_PiPPi0",&t_EventHasKaonP_PiPPi0);
   OutputTree->Branch("EventHasKaonM",&t_EventHasKaonM);
   OutputTree->Branch("EventHasKaon0",&t_EventHasKaon0);

   OutputTree->Branch("Neutrino","vector<SimParticle>",&t_Neutrino);
   OutputTree->Branch("Lepton","vector<SimParticle>",&t_Lepton);
   OutputTree->Branch("PrimaryHyperon","vector<SimParticle>",&t_PrimaryHyperon);
   OutputTree->Branch("PrimaryNucleon","vector<SimParticle>",&t_PrimaryNucleon);
   OutputTree->Branch("PrimaryPion","vector<SimParticle>",&t_PrimaryPion);
   OutputTree->Branch("PrimaryKaon","vector<SimParticle>",&t_PrimaryKaon);
   OutputTree->Branch("PrimaryKaonP","vector<SimParticle>",&t_PrimaryKaonP);
   OutputTree->Branch("PrimaryNucleus","vector<SimParticle>",&t_PrimaryNucleus);
   OutputTree->Branch("HyperonDecay","vector<SimParticle>",&t_HyperonDecay);
   OutputTree->Branch("KaonPDecay","vector<SimParticle>",&t_KaonPDecay);
   OutputTree->Branch("KaonPDecay_NuMuP","vector<SimParticle>",&t_KaonPDecay_NuMuP);
   OutputTree->Branch("KaonPDecay_PiPPi0","vector<SimParticle>",&t_KaonPDecay_PiPPi0);
   OutputTree->Branch("NeutralKaonDecayK0SL","vector<SimParticle>",&t_NeutralKaonDecayK0SL);

   OutputTree->Branch("TruePrimaryVertex_X",&t_TruePrimaryVertex_X);
   OutputTree->Branch("TruePrimaryVertex_Y",&t_TruePrimaryVertex_Y);
   OutputTree->Branch("TruePrimaryVertex_Z",&t_TruePrimaryVertex_Z);

   OutputTree->Branch("DecayVertex_X",&t_DecayVertex_X);
   OutputTree->Branch("DecayVertex_Y",&t_DecayVertex_Y);
   OutputTree->Branch("DecayVertex_Z",&t_DecayVertex_Z);

   OutputTree->Branch("RecoPrimaryVertex","TVector3",&t_RecoPrimaryVertex);
   OutputTree->Branch("NPrimaryTracks",&t_NPrimaryTracks);
   OutputTree->Branch("NPrimaryTrackDaughters",&t_NPrimaryTrackDaughters);
   OutputTree->Branch("NPrimaryShowerDaughters",&t_NPrimaryShowerDaughters);
   OutputTree->Branch("NOtherTracks",&t_NOtherTracks);
   OutputTree->Branch("NOtherRebuiltTracks",&t_NOtherRebuiltTracks);
   OutputTree->Branch("NOtherShowers",&t_NOtherShowers);
   
   OutputTree->Branch("TracklikePrimaryDaughters","vector<RecoParticle>",&t_TrackPrimaryDaughters);
   OutputTree->Branch("ShowerlikePrimaryDaughters","vector<RecoParticle>",&t_ShowerPrimaryDaughters);
   OutputTree->Branch("TrackOthers","vector<RecoParticle>",&t_TrackOthers);
   OutputTree->Branch("TrackRebuiltOthers","vector<RecoParticle>",&t_TrackRebuiltOthers);
   OutputTree->Branch("ShowerOthers","vector<RecoParticle>",&t_ShowerOthers);
   OutputTree->Branch("RepassTracklikePrimaryDaughters","vector<RecoParticle>",&t_RepassTrackPrimaryDaughters);
   OutputTree->Branch("RepassShowerlikePrimaryDaughters","vector<RecoParticle>",&t_RepassShowerPrimaryDaughters);
   
   //OutputTree->Branch("SysDials",&t_SysDials);
   //OutputTree->Branch("SysWeights","vector<vector<vector<double>>>",&t_SysWeights);

   OutputTree->Branch("SysDials",&t_SysDials);
   OutputTree->Branch("SysWeights","vector<vector<double>>",&t_SysWeights);

   //////////////////////////////////////////
   //             Metadata Tree	           //
   //////////////////////////////////////////
   
   m_NEvents=0;
   m_NHyperons=0;
   m_NKaons=0;
   m_NKaonPs=0;
   m_NKaonPs_NuMuP=0;
   m_NKaonPs_PiPPi0=0;
   m_NSignal=0;
   m_NSignal_NuMuP=0;
   m_NSignal_PiPPi0=0;
   m_NGoodReco=0;
   m_POT=0;

   MetaTree=tfs->make<TTree>("MetaTree","Metadata Info Tree");

   MetaTree->Branch("NEvents",&m_NEvents);
   MetaTree->Branch("NHyperons",&m_NHyperons);
   MetaTree->Branch("NKaons",&m_NKaons);
   MetaTree->Branch("NKaonPs",&m_NKaonPs);
   MetaTree->Branch("NKaonPs_NuMuP",&m_NKaonPs_NuMuP);
   MetaTree->Branch("NKaonPs_PiPPi0",&m_NKaonPs_PiPPi0);
   MetaTree->Branch("NSignal",&m_NSignal);
   MetaTree->Branch("NSignal_NuMuP",&m_NSignal_NuMuP);
   MetaTree->Branch("NSignal_PiPPi0",&m_NSignal_PiPPi0);
   MetaTree->Branch("NGoodReco",&m_NGoodReco);

   MetaTree->Branch("POT",&m_POT);

   if(f_Debug) std::cout << "Finished begin job" << std::endl;

}

void cckaon::KaonNtuples::endJob()
{
   MetaTree->Fill();
}

void cckaon::KaonNtuples::beginSubRun(const art::SubRun& sr)
{
   if(f_Debug) std::cout << "Getting Subrun POT Info" << std::endl;

   art::Handle<sumdata::POTSummary> POTHandle;
   if(sr.getByLabel(f_POTSummaryLabel,POTHandle)) m_POT += POTHandle->totpot;	
}

void cckaon::KaonNtuples::endSubRun(const art::SubRun& sr){}

DEFINE_ART_MODULE(cckaon::KaonNtuples)
