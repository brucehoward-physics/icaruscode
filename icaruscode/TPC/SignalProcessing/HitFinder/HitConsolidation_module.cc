///////////////////////////////////////////////////////////////////////
// Class:       HitConsolidation
// Plugin Type: producer (art v3_06_03)
// File:        HitConsolidation_module.cc
//
// Bruce Howard - 2022: copied from a Gauss Hit Recovery Module and edited
//                      to just save all the hits from an input list into 1
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

#include "lardata/ArtDataHelper/HitCreator.h"

/*
#include "art/Utilities/Globals.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "art_root_io/TFileService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PCAxis.h"

#include "lardata/ArtDataHelper/HitCreator.h"

#include "cetlib_except/exception.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
*/

// C++
#include <map>
#include <vector>
#include <iostream>
#include <string>
//#include <cmath>
//#include <set>
//#include <memory>
//#include <utility>
//#include <algorithm>
//#include <tuple>

// ROOT
//#include "TVector3.h"
//#include "TMath.h"

class HitConsolidation;


class HitConsolidation : public art::EDProducer {
public:
  explicit HitConsolidation(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HitConsolidation(HitConsolidation const&) = delete;
  HitConsolidation(HitConsolidation&&) = delete;
  HitConsolidation& operator=(HitConsolidation const&) = delete;
  HitConsolidation& operator=(HitConsolidation&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  // Params
  std::vector<art::InputTag> fHitsLabelVec;
  std::string                fHitCreatorInstanceName;
};


HitConsolidation::HitConsolidation(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fHitsLabelVec           ( p.get< std::vector<art::InputTag> >("HitModuleLabels",        std::vector<art::InputTag>() = {"gausHit"}) ),
  fHitCreatorInstanceName ( p.get<std::string>("HitCreatorInstanaceName","") )
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  // TODO: is this the right thing to ask here? I think I'm assuming this isn't thread-safe but maybe we want it to be?
  // The GausHitFinder had a check related to thread-safety that I'm basing this structure on
  if (art::Globals::instance()->nthreads() > 1u) {
    throw art::Exception(art::errors::Configuration)
      << "I think it's probably safe to say this module isn't configured to run multi-thread...";
  }

  // basically from the GausHitFinder
  recob::HitCollectionCreator::declare_products(producesCollector(), fHitCreatorInstanceName, true, false);
}

void HitConsolidation::produce(art::Event& e)
{
  // As in Gauss Hit Finder code (larreco/HitFinder/GausHitFinder_module.cc)
  recob::HitCollectionCreator hitCol(e, fHitCreatorInstanceName, true, false);
  // Also as in GausHitFinder
  struct hitstruct {
    recob::Hit stHit;
    art::Ptr<recob::Wire> stWire;
  };

  for ( auto const& iLabel : fHitsLabelVec ) {
    // Get the hits
    art::Handle< std::vector<recob::Hit> > hitsHandle;
    std::vector< art::Ptr<recob::Hit> > hits;
    if ( e.getByLabel(iLabel,hitsHandle) ) {
      art::fill_ptr_vector(hits,hitsHandle);
    }
    else{
      mf::LogWarning("HitConsolidation") << "Event failed to find recob::Hit with label " << iLabel;
      return;
    }

    // And the association to the wire
    art::FindManyP<recob::Wire> fmwire(hitsHandle, e, iLabel);
    if( !fmwire.isValid() ){
      mf::LogError("HitConsolidation") << "Error in validity of fmwire. Returning.";
      return;
    }

    for ( auto const& iHitPtr : hits ) {
      recob::Hit theHit = *iHitPtr;
      // And get the associated wire
      std::vector< art::Ptr<recob::Wire> > hitWires = fmwire.at(iHitPtr.key());
      if ( hitWires.size() == 0 ) throw cet::exception("HitConsolidation") << "Hit found with no associated wires...\n";
      else if ( hitWires.size() > 1 ) mf::LogWarning("HitConsolidation") << "Hit with >1 recob::Wire associated...";

      hitCol.emplace_back(theHit,hitWires[0]);
    }
  }

  hitCol.put_into(e);
}

DEFINE_ART_MODULE(HitConsolidation)
