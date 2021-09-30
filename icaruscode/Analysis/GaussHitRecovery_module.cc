///////////////////////////////////////////////////////////////////////
// Class:       GaussHitRecovery
// Plugin Type: producer (art v3_06_03)
// File:        GaussHitRecovery_module.cc
//
// Generated at Wed Sep 29 14:58:50 2021 by Bruce Howard using cetskelgen
// from cetlib version v3_11_01.
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

#include "art/Utilities/Globals.h"
#include "larcore/Geometry/Geometry.h"
#include "art_root_io/TFileService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/ArtDataHelper/HitCreator.h"

// C++
#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <unordered_set>
#include <memory>
#include <utility>

class GaussHitRecovery;


class GaussHitRecovery : public art::EDProducer {
public:
  explicit GaussHitRecovery(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GaussHitRecovery(GaussHitRecovery const&) = delete;
  GaussHitRecovery(GaussHitRecovery&&) = delete;
  GaussHitRecovery& operator=(GaussHitRecovery const&) = delete;
  GaussHitRecovery& operator=(GaussHitRecovery&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  // Setup
  geo::GeometryCore const* fGeom;

  // Params
  std::vector<art::InputTag> fHitsLabelVec;
  std::vector<art::InputTag> fReducHitsLabelVec;
  std::string fHitCreatorInstanceName;
  int fNTicksSmooth;
  int fNWiresSmooth;
  int fReqNeighbors;
  float fMaxChi2NDF;
  float fLinRejTolerance;
  int fSkipRecoveryPlane;
};


GaussHitRecovery::GaussHitRecovery(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fHitsLabelVec           ( p.get< std::vector<art::InputTag> >("HitModuleLabels",        std::vector<art::InputTag>() = {"gausHit"}) ),
  fReducHitsLabelVec      ( p.get< std::vector<art::InputTag> >("ReducedHitModuleLabels", std::vector<art::InputTag>() = {"cluster3D"}) ),
  fHitCreatorInstanceName ( p.get<std::string>("HitCreatorInstanaceName","") ),
  fNTicksSmooth           ( p.get<int>("NTicksSmooth",60) ),
  fNWiresSmooth           ( p.get<int>("NWiresSmooth",12) ),
  fReqNeighbors           ( p.get<int>("ReqNeighbors",4) ),
  fMaxChi2NDF             ( p.get<float>("MaxChi2NDF",11.) ),
  fLinRejTolerance        ( p.get<float>("LinRejTolerance",1.) ),
  fSkipRecoveryPlane      ( p.get<float>("SkipRecoveryPlane",-1) )
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  // NOTE ?? - is it not okay to do what I'm doing below to get the hits to decide amongst?

  // Following exception and declaration come from the GausHitFinder, with some alterations as needed.
  if (art::Globals::instance()->nthreads() > 1u) {
    throw art::Exception(art::errors::Configuration)
      << "I think it's probably safe to say this module isn't configured to run multi-thread...";
  }

  recob::HitCollectionCreator::declare_products(producesCollector(), fHitCreatorInstanceName, false, false);
}

void GaussHitRecovery::produce(art::Event& e)
{
  // -- As in Gauss Hit Finder code (larreco/HitFinder/GausHitFinder_module.cc) but with assns turned off
  recob::HitCollectionCreator hitCol(e, fHitCreatorInstanceName, false, false);

  // -- Method explored in the HitRecoveryAnalysis module
  // Geometry
  fGeom = lar::providerFrom<geo::Geometry>();

  // Recover hits - first look at the reduced hits for linear patterns then find hits from the whole set that match this.
  // For all hits: use the hit->Channel->WireIDVector to get the possible matches
  // For reduced hits: since a choice of the hit location has already been made in the 3d match, here we just use the hit->WireID

  // TODO: right now this means recovered hits can be double-counted... Would need to add some way to force only picking one. But probably a hit being claimed in both
  //       possibilities would be pretty small?

  std::map< geo::WireID, std::vector< std::pair<float,float> > > wireToReducHitsTimeRMS;

  std::map< geo::PlaneID, std::vector< std::pair<double,double> > > linrej;

  for ( auto const& iLabel : fReducHitsLabelVec ) {
    art::Handle< std::vector<recob::Hit> > hitsHandle;
    std::vector< art::Ptr<recob::Hit> > hits;
    if ( e.getByLabel(iLabel,hitsHandle) ) {
      art::fill_ptr_vector(hits,hitsHandle);
    }
    else{
      mf::LogWarning("HitRecoverAna") << "Event failed to find recob::Hit with label " << iLabel;
      return;
    }

    for ( auto const& iHitPtr : hits ) {
      // All reduced hits should be saved in the final hit set
      recob::Hit theHit = *iHitPtr;
      hitCol.emplace_back(std::move(theHit));

      // Add the hit to this map for hopefully easier lookup later.
      wireToReducHitsTimeRMS[ iHitPtr->WireID() ].push_back( std::make_pair(iHitPtr->PeakTime(), iHitPtr->RMS()) );

      std::vector<int> wire_list;
      std::vector<float> tick_list;
      std::vector<float> rms_list;

      auto const& thisTime = iHitPtr->PeakTime();
      geo::WireID const& wID = iHitPtr->WireID();
      auto const& thisCryo = wID.Cryostat;
      auto const& thisTPC = wID.TPC;
      auto const& thisPlane = wID.Plane;
      auto const& thisWire = wID.Wire;

      if ( fSkipRecoveryPlane == (int)thisPlane ) continue;

      // Figure out how many neighbors this hit has, and if it has enough, see if they are a fairly linear grouping
      int nNeighbors = -1; // start at -1 since we will count ourselves later
      for ( auto const& jHitPtr : hits ) {
	auto const& jTime = jHitPtr->PeakTime();
	auto const& jRMS = jHitPtr->RMS();
	geo::WireID const& jWID = jHitPtr->WireID();
	auto const& jCryo = jWID.Cryostat;
	auto const& jTPC = jWID.TPC;
	auto const& jPlane = jWID.Plane;
	auto const& jWire = jWID.Wire;
	if( jCryo!=thisCryo || jTPC!=thisTPC || jPlane!=thisPlane ) continue;
	if( std::fabs(thisWire-jWire)<fNWiresSmooth && std::fabs(thisTime-jTime)<fNTicksSmooth ){
	  nNeighbors+=1;
	  wire_list.push_back(jWire);
	  tick_list.push_back(jTime);
	  rms_list.push_back(jRMS);
	}
      }

      if( nNeighbors < fReqNeighbors ) continue;

      // Do a linear regression to see if this group of hits is "line-like" (aka track-like)
      // See especially:
      //   https://en.wikipedia.org/wiki/Simple_linear_regression
      //   https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
      // Doesn't account for hit RMS's here, just uses the peak times
      // -- Get m and b
      double xave = 0.;
      double yave = 0.;
      for( unsigned int idx=0; idx<wire_list.size(); idx++ ){
	xave+=wire_list[idx];
	yave+=tick_list[idx];
      }
      xave/=double(wire_list.size());
      yave/=double(wire_list.size());
      // slope: numerator is s_{x,y}, denominator s_{x^2}
      double s_xy = 0.;
      double s_x2 = 0.;
      for( unsigned int idx=0; idx<wire_list.size(); idx++ ){
	s_xy+=((wire_list[idx]-xave)*(tick_list[idx]-yave));
	s_x2+=((wire_list[idx]-xave)*(wire_list[idx]-xave));
      }
      double slope = s_xy / s_x2;
      double intercept = yave - slope*xave;
      
      // -- Get Chi2/N.D.F.
      // assumes sigma=hit RMS
      double chi2ndf = 0.;
      for( unsigned int idx=0; idx<wire_list.size(); idx++ ){
	double expect=slope*wire_list[idx]+intercept;
	chi2ndf+=((tick_list[idx]-expect)*(tick_list[idx]-expect)/(rms_list[idx]*rms_list[idx]));
      }
      chi2ndf/=(double(wire_list.size()-2.));
      
      if( chi2ndf > fMaxChi2NDF ) continue;
      
      linrej[ wID.asPlaneID() ].push_back( std::make_pair(slope,intercept) );
    } // loop hits
  } // loop hit labels

  // All hits, save some as the recovered hits
  for ( auto const& iLabel : fHitsLabelVec ) {
    art::Handle< std::vector<recob::Hit> > hitsHandle;
    std::vector< art::Ptr<recob::Hit> > hits;
    if ( e.getByLabel(iLabel,hitsHandle) ) {
      art::fill_ptr_vector(hits,hitsHandle);
    }
    else{
      mf::LogWarning("HitRecoverAna") << "Event failed to find recob::Hit with label " << iLabel;
      return;
    }

    for ( auto const& iHitPtr : hits ) {
      auto const& thisTime = iHitPtr->PeakTime();
      auto const& thisRMS = iHitPtr->RMS();
      std::vector< geo::WireID > wires = fGeom->ChannelToWire( iHitPtr->Channel() );
      for ( auto const& iWire : wires ) {
	auto const& thisWire = iWire.Wire;
	auto const& thisPlane = iWire.Plane;
	
	if ( fSkipRecoveryPlane == (int)thisPlane ) continue;

	// if no linrej on this plane then skip
	if( linrej.find(iWire.asPlaneID()) == linrej.end() ) continue;

	// loop through the linrej for this PlaneID and see if the hit is consistent with one of them
	for ( auto const& iLine : linrej[iWire.asPlaneID()] ) {
	  double expect = iLine.first*thisWire + iLine.second;
	  if( thisTime+fLinRejTolerance*thisRMS > expect && thisTime-fLinRejTolerance*thisRMS < expect ){
	    // Check if the hit is already in the set of cluster3d hits and skip it if so
	    bool isReduc = false;
	    if ( !(wireToReducHitsTimeRMS.find(iWire) == wireToReducHitsTimeRMS.end()) ) {
	      // Now check this wire's hits for something with the same tick and RMS
	      for ( auto const& iTimeRMS : wireToReducHitsTimeRMS[iWire] ){
		if ( thisTime == iTimeRMS.first && thisRMS == iTimeRMS.second ) {
		  isReduc = true;
		  break;
		}
	      }
	    }
	    if ( isReduc ) break;
	    // Put the hit into the recovered hit vector
	    recob::Hit theHit = *iHitPtr;
	    hitCol.emplace_back(std::move(theHit));
	    break;
	  }
	}
      }
    } // loop hits
  } // loop hit labels

  hitCol.put_into(e);
}

DEFINE_ART_MODULE(GaussHitRecovery)
