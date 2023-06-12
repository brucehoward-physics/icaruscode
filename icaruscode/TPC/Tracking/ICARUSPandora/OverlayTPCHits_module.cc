////////////////////////////////////////////////////////////////////////
// Class:       OverlayTPCHits
// Plugin Type: producer (Unknown Unknown)
// File:        OverlayTPCHits_module.cc
//
// Generated at Sun Jun 11 16:23:07 2023 by Bruce Howard using cetskelgen
// from  version .
//
// Useful for overlaying hits from 2 sets in the context of HARPS studies...
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

#include "cetlib_except/exception.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"

#include <memory>

class OverlayTPCHits;


class OverlayTPCHits : public art::EDProducer {
public:
  explicit OverlayTPCHits(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OverlayTPCHits(OverlayTPCHits const&) = delete;
  OverlayTPCHits(OverlayTPCHits&&) = delete;
  OverlayTPCHits& operator=(OverlayTPCHits const&) = delete;
  OverlayTPCHits& operator=(OverlayTPCHits&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  std::vector<art::InputTag> fOverlayHitModuleLabels;
  bool fTPCHitsWireAssn;
  std::string fTPCHitCreatorInstanceName;

};


OverlayTPCHits::OverlayTPCHits(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fOverlayHitModuleLabels    ( p.get< std::vector<art::InputTag> >("OverlayHitModuleLabels", {"HARPS", "gaushitTPCEE", "gaushitTPCEW"} ) ),
  fTPCHitsWireAssn           ( p.get< bool >("TPCHitsWireAssn", true) ),
  fTPCHitCreatorInstanceName ( p.get<std::string>("TPCHitCreatorInstanaceName", "") )
{
  recob::HitCollectionCreator::declare_products(producesCollector(), fTPCHitCreatorInstanceName, fTPCHitsWireAssn, false);
}

void OverlayTPCHits::produce(art::Event& e)
{
  // As in Gauss Hit Finder code (larreco/HitFinder/GausHitFinder_module.cc)
  recob::HitCollectionCreator hitCol(e, fTPCHitCreatorInstanceName, fTPCHitsWireAssn, false);

  for ( auto const& label : fOverlayHitModuleLabels ) {
    art::Handle< std::vector<recob::Hit> > hitHandle;
    std::vector< art::Ptr<recob::Hit> > hits;
    if ( e.getByLabel(label, hitHandle) ) {
      art::fill_ptr_vector(hits, hitHandle);
    }
    else {
      mf::LogError("TPCHitOverlay") << "Error pulling in Hit product... Skipping cryostat.";
    }

    // For the Wire Assns in HitCreator
    art::FindManyP<recob::Wire> fmWire(hitHandle, e, label);
    if ( !fmWire.isValid() && fTPCHitsWireAssn ){
      throw cet::exception("TPCHitOverlay") << "Module wants to use Hit-Wire Associations, but fmWire invalid." << std::endl;
    }

    for ( auto const& hit : hits ) {
      recob::Hit theHit = *hit;
      if ( fTPCHitsWireAssn ) {
        std::vector< art::Ptr<recob::Wire> > hitWires = fmWire.at(hit.key());
        if ( hitWires.size() == 0 ) throw cet::exception("TPCHitOverlay") << "Hit found with no associated wires...\n";
        else if ( hitWires.size() > 1 ) mf::LogWarning("TPCHitOverlay") << "Hit with >1 recob::Wire associated...";
        hitCol.emplace_back(theHit,hitWires[0]);
      }
      else hitCol.emplace_back(theHit);
    }
  }

  hitCol.put_into(e);
}

DEFINE_ART_MODULE(OverlayTPCHits)
