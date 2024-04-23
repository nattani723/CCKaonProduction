////////////////////////////////////////////////////////////////////////
// Class:       DirectLambdaSignalFilter
// Plugin Type: filter (art v3_01_02)
// File:        DirectLambdaSignalFilter_module.cc
//
// Purpose: Finds events containing hyperon production (direct and associated)
//
// Generated at Thu Sep 16 09:04:20 2021 by Christopher Thorpe using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "ubana/CCKaonProduction/Modules/SubModules/SubModuleGeneratorTruth.h"

namespace hyperon {
   class DirectLambdaSignalFilter;
}


class hyperon::DirectLambdaSignalFilter : public art::EDFilter {
   public:
      explicit DirectLambdaSignalFilter(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      DirectLambdaSignalFilter(DirectLambdaSignalFilter const&) = delete;
      DirectLambdaSignalFilter(DirectLambdaSignalFilter&&) = delete;
      DirectLambdaSignalFilter& operator=(DirectLambdaSignalFilter const&) = delete;
      DirectLambdaSignalFilter& operator=(DirectLambdaSignalFilter&&) = delete;

      // Required functions.
      bool filter(art::Event& e) override;

   private:

      fhicl::ParameterSet f_G4;
      double f_DecayProtonThresh;
      double f_DecayPionThresh;
};


hyperon::DirectLambdaSignalFilter::DirectLambdaSignalFilter(fhicl::ParameterSet const& p)
   : EDFilter{p},
   f_G4(p.get<fhicl::ParameterSet>("G4"))
{
}

bool hyperon::DirectLambdaSignalFilter::filter(art::Event& e)
{
bool pass = false;
/*      
      SubModuleG4Truth* G4_SM = new SubModuleG4Truth(e,f_G4);
      G4Truth G4T = G4_SM->GetG4Info();
      //pass = G4T.EventHasNeutronScatter;       

      int signalindex = -1;
      for(size_t i_mct=0;i_mct<G4T.IsSignal.at        

      delete G4_SM;
  */
      return pass;
}

DEFINE_ART_MODULE(hyperon::DirectLambdaSignalFilter)
