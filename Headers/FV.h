#ifndef _FV_h_
#define _FV_h_

#include "TVector3.h"

namespace cckaon {

const std::vector<double> TPCCenter = { 126.625 , 0.97 , 518.5 }; //center of active TPC
const std::vector<double> TPCSideLengths = { 236.35 , 233.0 , 1036.8 }; //side lengths of active TPC

// Use whole TPC - decide which volumes to cut later
const double FVxmin = 0.0;
const double FVxmax = 256.35;
const double FVymin = -115.53;
const double FVymax = 117.47;
const double FVzmin = 0.1;
const double FVzmax = 1036.9;

const double FVxmin_5cm = 5.0;
const double FVxmax_5cm = 251.35;
const double FVymin_5cm = -110.53;
const double FVymax_5cm = 112.47;
const double FVzmin_5cm = 5.1;
const double FVzmax_5cm = 1031.9;

const double FVxmin_CCInclusive = 0.0;
const double FVxmax_CCInclusive = 246.35;
const double FVymin_CCInclusive = -105.53;
const double FVymax_CCInclusive = 107.47;
const double FVzmin_CCInclusive = 10.1;
const double FVzmax_CCInclusive = 986.8;


inline bool inActiveTPC(TVector3 pos){
   if(pos.X() > FVxmax || pos.X() < FVxmin) return false;
   if(pos.Y() > FVymax || pos.Y() < FVymin) return false;
   if(pos.Z() > FVzmax || pos.Z() < FVzmin) return false;
   return true;
}

inline bool inActiveTPC5cm(TVector3 pos){
   if(pos.X() > FVxmax_5cm || pos.X() < FVxmin_5cm) return false;
   if(pos.Y() > FVymax_5cm || pos.Y() < FVymin_5cm) return false;
   if(pos.Z() > FVzmax_5cm || pos.Z() < FVzmin_5cm) return false;
   return true;
}

inline bool inActiveTPCCCInclusive(TVector3 pos){
   if(pos.X() > FVxmax_CCInclusive || pos.X() < FVxmin_CCInclusive) return false;
   if(pos.Y() > FVymax_CCInclusive || pos.Y() < FVymin_CCInclusive) return false;
   if(pos.Z() > FVzmax_CCInclusive || pos.Z() < FVzmin_CCInclusive) return false;
   return true;
}


}


#endif




