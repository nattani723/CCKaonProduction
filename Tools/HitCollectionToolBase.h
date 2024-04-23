#ifndef _HitCollectionToolBase_h_
#define _HitCollectionToolBase_h_

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"				
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "ubana/KReco/TrackRebuilder/TrackRebuilder.h"



namespace hyperon {

  class HitCollectionToolBase {
 
  public: 

  HitCollectionToolBase(const fhicl::ParameterSet& p) : params(p) {}
  //virtual ~HitCollectionToolBase(); 

  // makes a vector of the collections required to generate a new track - each entry is for a separate track
  virtual void MakeHitCollections(art::Event const& e,
                                  std::vector<std::vector<art::Ptr<recob::Hit>>>& r_hits,
                                  std::vector<std::map<art::Ptr<recob::Hit>,art::Ptr<recob::SpacePoint>>>& r_hitspacepointmap,
                                  std::vector<pandora::CartesianVector>& r_vertex) const = 0;
  
  protected:
   
  const fhicl::ParameterSet params;  

  };

}

#endif
