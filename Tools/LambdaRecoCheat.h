#ifndef _LambdaRecoCheat_h_
#define _LambdaRecoCheat_h_

#include "HitCollectionToolBase.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

namespace hyperon {

  class LambdaRecoCheat : public HitCollectionToolBase {

  public: 

  LambdaRecoCheat(const fhicl::ParameterSet& p) : HitCollectionToolBase(p) , MinEnergyDep(p.get<double>("MinEnergyDep",0.0)), MinEnergyADCRatio(p.get<double>("MinEnergyADCRatio",0.0)) {}

  void MakeHitCollections(art::Event const& e,
                                  std::vector<std::vector<art::Ptr<recob::Hit>>>& r_hits,
                                  std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>>& r_hitspacepointmap,
                                  std::vector<pandora::CartesianVector>& r_vertex) const;

  private:

  const double MinEnergyDep;
  const double MinEnergyADCRatio;

  };

}

#endif


