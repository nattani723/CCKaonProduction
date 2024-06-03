////////////////////////////////////////////////////////////////////////
// Class:       HyperonFilter
// Plugin Type: filter (art v3_01_02)
// File:        HyperonFilter_module.cc
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

namespace cckaon {
   class KaonPFilter;
}


class cckaon::KaonPFilter : public art::EDFilter {
   public:
      explicit KaonPFilter(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      KaonPFilter(KaonPFilter const&) = delete;
      KaonPFilter(KaonPFilter&&) = delete;
      KaonPFilter& operator=(KaonPFilter const&) = delete;
      KaonPFilter& operator=(KaonPFilter&&) = delete;

      // Required functions.
      bool filter(art::Event& e) override;

   private:

      fhicl::ParameterSet f_Generator;
};


cckaon::KaonPFilter::KaonPFilter(fhicl::ParameterSet const& p)
   : EDFilter{p},
   f_Generator(p.get<fhicl::ParameterSet>("Generator"))
{
}

bool cckaon::KaonPFilter::filter(art::Event& e)
{
      
      SubModuleGeneratorTruth* Generator_SM = new SubModuleGeneratorTruth(e,f_Generator);
      GeneratorTruth GenT = Generator_SM->GetGeneratorTruth();

      bool EventHasKaonP = GenT.EventHasKaonP;       

      delete Generator_SM;
  
      return EventHasKaonP;
}

DEFINE_ART_MODULE(cckaon::KaonPFilter)
