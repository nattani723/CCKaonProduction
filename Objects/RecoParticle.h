#ifndef _RecoParticle_h_
#define _RecoParticle_h_

#include <iostream>

#include "TLorentzVector.h"
#include "TVector3.h"

#ifdef __MAKE_ROOT_DICT__
#include "TObject.h"
#endif

#ifdef __MAKE_ROOT_DICT__
class RecoParticle : public TObject{
#else
class RecoParticle {
#endif

public:

RecoParticle(){}
~RecoParticle(){}

int Index;
bool InNuSlice = false;
bool IsRebuilt = false;

// General reco info
int PDG; // Pandora PDG code (11 or 13)
int Parentage; // 1 - neutrino daughter, 2 - neutrino granddaughter, 3 - other
int ParentIndex=-1; // -1 - neutrino candidate or no parent
double TrackShowerScore;
double X,Y,Z;
double Displacement; // Distance from PV

// Track variables
double TrackLength=0;
double TrackDirectionX=0,TrackDirectionY=0,TrackDirectionZ=0;
double TrackStartX=0,TrackStartY=0,TrackStartZ=0;
double TrackEndX=0,TrackEndY=0,TrackEndZ=0;
double TrackPID; // 3 plane PID score
double MeandEdX_Plane0,MeandEdX_Plane1,MeandEdX_Plane2,MeandEdX_ThreePlane; // Mean dE/dX scores
double Track_LLR_PID; // LLR PID
double Track_LLR_PID_Kaon; // LLR PID with Kaon hypothesis
double Track_LLR_PID_Kaon_Partial; // LLR PID with Kaon hypothesis using last 5cm of track
double Track_Bragg_PID_Kaon;
double Track_Chi2_Kaon_Plane0, Track_Chi2_Kaon_Plane1, Track_Chi2_Kaon_Plane2, Track_Chi2_Kaon_3Plane;
double Track_Chi2_Proton_Plane0, Track_Chi2_Proton_Plane1, Track_Chi2_Proton_Plane2, Track_Chi2_Proton_3Plane;
double Track_Chi2_Pion_Plane0, Track_Chi2_Pion_Plane1, Track_Chi2_Pion_Plane2, Track_Chi2_Pion_3Plane;
double Track_Chi2_Muon_Plane0, Track_Chi2_Muon_Plane1, Track_Chi2_Muon_Plane2, Track_Chi2_Muon_3Plane;
double ProtonMomentum,MuonMomentum,KaonMomentum; // Track kinematics
double TrackWiggliness;

// Truth info
bool HasTruth; // False if reco particle has no corresponding MC particle
int MCTruthIndex=-1;
int TrackTruePDG;
double TrackTrueE,TrackTruePx,TrackTruePy,TrackTruePz;
double TrackTrueEndE,TrackTrueEndPx,TrackTrueEndPy,TrackTrueEndPz;
double TrackTrueModMomentum;
double TrackTrueEndModMomentum;
double TrackTrueKE;
double TrackTrueEndKE;
double TrackTrueLength;
int TrackTrueOrigin; // 0 = neutrino , 1 = primary , 2 = hyperon decay, 3 = other, 4 = Kaon decay, 5 = Sigma0 decay, 6 = K0S/K0L, 7 = K0S/K0L decay
double TrackTruthPurity;

// Merge Check
int MergePDG_1st, MergePDG_2nd, MergePDG_3rd;
double MergeEnergyPurity_1st, MergeEnergyPurity_2nd, MergeEnergyPurity_3rd;
double MergeHitPurity_1st, MergeHitPurity_2nd, MergeHitPurity_3rd;

inline void SetVertex(TVector3 V);
inline void SetTrackPositions(TVector3 Start,TVector3 End);
inline void Print();

#ifdef __MAKE_ROOT_DICT__
ClassDef(RecoParticle,1);
#endif

};


inline void RecoParticle::SetVertex(TVector3 V){

   X = V.X();
   Y = V.Y();
   Z = V.Z();

}

inline void RecoParticle::SetTrackPositions(TVector3 Start,TVector3 End){

   TrackStartX = Start.X();
   TrackStartY = Start.Y();
   TrackStartZ = Start.Z();

   TrackEndX = End.X();
   TrackEndY = End.Y();
   TrackEndZ = End.Z();

}

inline void RecoParticle::SetMergeCheck(std::vector<int> MergePDG, std::vector<double> MergeEnergyPurity, std::vector<double> MergeHitPurity){
   
      MergePDG_1st = MergePDG.at(0);
      MergePDG_2nd = MergePDG.at(1);
      MergePDG_3rd = MergePDG.at(2);

      MergeEnergyPurity_1st = MergeEnergyPurity.at(0);
      MergeEnergyPurity_2nd = MergeEnergyPurity.at(1);
      MergeEnergyPurity_3rd = MergeEnergyPurity.at(2);

      MergeHitPurity_1st = MergeHitPurity.at(0);
      MergeHitPurity_2nd = MergeHitPurity.at(1);
      MergeHitPurity_3rd = MergeHitPurity.at(2);
   
}

inline void RecoParticle::Print(){

   std::cout << "Reco Info:" << std::endl;
   std::cout << "PDG Code: " << PDG << "  Track/Shower score: " << TrackShowerScore << std::endl;
   std::cout << "Track length: " << TrackLength << "  PID score: " << TrackPID <<  std::endl;
   std::cout << "Truth Info:" << std::endl;
   std::cout << "PDG: " << TrackTruePDG << "  Origin: " << TrackTrueOrigin << std::endl;
   std::cout << "Length: " << TrackTrueLength << "  KE: " << TrackTrueKE << std::endl;

}

#endif
