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

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "art/Utilities/Globals.h"
#include "larcore/Geometry/Geometry.h"
#include "art_root_io/TFileService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PCAxis.h"

#include "lardata/ArtDataHelper/HitCreator.h"

#include "cetlib_except/exception.h"

// C++
#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <unordered_set>
#include <memory>
#include <utility>

// ROOT
#include "TVector3.h"
#include "TMath.h"

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
  std::vector<art::InputTag> fWireLabelVec;
  std::string fHitCreatorInstanceName;
  int fNTicksSmooth;         // in neighbors method, range of ticks to look around
  int fNWiresSmooth;         // in neighbors method, range of wires to look around
  int fReqNeighbors;         // in neighbors method, number of neighbors required of "pseudo-cluster"

  float fMaxChi2NDF;         // max chi2 per ndf at/above which do not save the given linear regression
  float fLinregTolerance;    // allowed tolerance (in number of RMSs) for a hit from considered line (all 3 methods) to recover
  int fSkipRecoveryPlane;    // set this value to a number (typically 0, 1, 2) to skip trying to recover hits from that plane

  bool fUseFindNeighbors;    // Default method. Look for interesting points in the reduced hit set.
  bool fUseReducHitClusters; // Instead of looking for the interesting points in the hit map, take an input set of clusters.
  bool fUseReducHitPCAxis;   // Third option - use input set of PCAxis objects.
  bool fPCAxisOnlyBetween;   // PCAxis option -- If true, only recover hits between 2 groupings on a plane, not along the extension
  float fPCAxisTolerance;    // PCAxis option -- Defines the anglular tolerance (in degrees) in PCAxis method.

  int fMinMethodsSuccess;    // how many methods does a hit need to pass to be saved
                             // This plus true/false for the three methods allows you to AND or OR the methods essentially
                             // (default "0" - a hit still has to be registed to save though)
};


GaussHitRecovery::GaussHitRecovery(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fHitsLabelVec           ( p.get< std::vector<art::InputTag> >("HitModuleLabels",        std::vector<art::InputTag>() = {"gausHit"}) ),
  fReducHitsLabelVec      ( p.get< std::vector<art::InputTag> >("ReducedHitModuleLabels", std::vector<art::InputTag>() = {"cluster3D"}) ),
  fWireLabelVec           ( p.get< std::vector<art::InputTag> >("WireLabel",              std::vector<art::InputTag>() = {"roifinder"}) ),
  fHitCreatorInstanceName ( p.get<std::string>("HitCreatorInstanaceName","") ),
  fNTicksSmooth           ( p.get<int>("NTicksSmooth",60) ),
  fNWiresSmooth           ( p.get<int>("NWiresSmooth",12) ),
  fReqNeighbors           ( p.get<int>("ReqNeighbors",4) ),
  fMaxChi2NDF             ( p.get<float>("MaxChi2NDF",11.) ),
  fLinregTolerance        ( p.get<float>("LinregTolerance",1.) ),
  fSkipRecoveryPlane      ( p.get<float>("SkipRecoveryPlane",-1) ),
  fUseFindNeighbors       ( p.get<bool>("UseFindNeighbors",true) ),
  fUseReducHitClusters    ( p.get<bool>("UseReducedHitClusters",false) ),
  fUseReducHitPCAxis      ( p.get<bool>("UseReducedHitPCAxis",false) ),
  fPCAxisOnlyBetween      ( p.get<bool>("PCAxisOnlyBetween",false) ),
  fPCAxisTolerance        ( p.get<float>("PCAxisTolerance",5.) ),
  fMinMethodsSuccess      ( p.get<int>("MinMethodsSuccess",0) )
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  const int usingNeighbors = fUseFindNeighbors    ? 1 : 0;
  const int usingClusters  = fUseReducHitClusters ? 1 : 0;
  const int usingPCAxis    = fUseReducHitPCAxis   ? 1 : 0;
  const int numMethods = usingNeighbors + usingClusters + usingPCAxis;
  if ( numMethods==0 )
    throw cet::exception("GaussHitRecovery") << "Trying to run HitRecovery with no methods activated... !\n";
  if ( numMethods < fMinMethodsSuccess )
    throw cet::exception("GaussHitRecovery") << "Too few methods run compared to ask for MinMethodsSuccess... !\n";

  // TODO: is this the right thing to ask here? I think I'm assuming this isn't thread-safe but maybe we want it to be?
  // The GausHitFinder had a check related to thread-safety that I'm basing this structure on
  if (art::Globals::instance()->nthreads() > 1u) {
    throw art::Exception(art::errors::Configuration)
      << "I think it's probably safe to say this module isn't configured to run multi-thread...";
  }

  // basically from the GausHitFinder
  recob::HitCollectionCreator::declare_products(producesCollector(), fHitCreatorInstanceName, true, false);
}

void GaussHitRecovery::produce(art::Event& e)
{
  // As in Gauss Hit Finder code (larreco/HitFinder/GausHitFinder_module.cc)
  recob::HitCollectionCreator hitCol(e, fHitCreatorInstanceName, true, false);
  // Also as in GausHitFinder
  struct hitstruct {
    recob::Hit stHit;
    art::Ptr<recob::Wire> stWire;
  };

  // Geometry
  fGeom = lar::providerFrom<geo::Geometry>();

  /////////////////////////////////////////////////////////////////////
  //
  // STEP 1: Deal with the reduced set of hits we want to build upon
  //
  /////////////////////////////////////////////////////////////////////

  std::map< geo::WireID, std::vector< std::pair<float,float> > > wireToReducHitsTimeRMS;

  for ( auto const& iLabel : fReducHitsLabelVec ) {
    // Get the various handles and fill vectors, FindMany, etc. - some checks on things. In some cases
    // these come from checks seen elsewhere (e.g. SBND analyzer module), but a note is that right now
    // the behavior is to print a warning and return in the case where a check fails. Should be considered
    // if throwing an exception would be better?
    art::Handle< std::vector<recob::Hit> > hitsHandle;
    std::vector< art::Ptr<recob::Hit> > hits;
    if ( e.getByLabel(iLabel,hitsHandle) ) {
      art::fill_ptr_vector(hits,hitsHandle);
    }
    else{
      mf::LogWarning("GaussHitRecovery") << "Event failed to find recob::Hit with label " << iLabel;
      return;
    }

    // std::vector< art::Handle< std::vector<recob::Wire> > > wiresHandleVec;
    // wiresHandleVec.resize( fWireLabelVec.size() );
    // for ( unsigned int jLabel=0; jLabel<fWireLabelVec.size(); ++jLabel ) {
    //   if ( !e.getByLabel(fWireLabelVec[jLabel],wiresHandleVec[jLabel]) ) {
    //     mf::LogWarning("GaussHitRecovery") << "Event failed to find recob::Wire with label " << fWireLabelVec[jLabel];
    //     return;
    //   }
    // }

    art::FindManyP<recob::Wire> fmwire(hitsHandle, e, iLabel);
    if( !fmwire.isValid() ){
      mf::LogError("GaussHitRecovery") << "Error in validity of fmwire. Returning.";
      return;
    }
    
    for ( auto const& iHitPtr : hits ) {
      // All reduced hits should be saved in the final hit set
      recob::Hit theHit = *iHitPtr;
      // And get the associated wire
      std::vector< art::Ptr<recob::Wire> > hitWires = fmwire.at(iHitPtr.key());
      if ( hitWires.size() == 0 ) throw cet::exception("GaussHitRecovery") << "Hit found with no associated wires...\n";
      else if ( hitWires.size() > 1 ) mf::LogWarning("GaussHitRecovery") << "Hit with >1 recob::Wire associated...";

      hitCol.emplace_back(theHit,hitWires[0]);
      // Add the hit to this map for hopefully easier lookup later.
      wireToReducHitsTimeRMS[ iHitPtr->WireID() ].push_back( std::make_pair(iHitPtr->PeakTime(), iHitPtr->RMS()) );
    }
  }

  /////////////////////////////////////////////////////////////////////
  //
  // STEP 2: Now get the properties from these reduced hits and go
  //         searching for hits to reclaim from the set of all hits
  //
  /////////////////////////////////////////////////////////////////////

  // TODO: is the RMS info needed along with the tick to guarantee we get a unique hit, i.e. not just wireid and tick?
  // Map of pair <WireID, tick> to a vector of the recovery methods (each one being a separate int put into the vector)
  std::map< std::pair< geo::WireID, float >, std::vector<int> > hitToRecoveryMethodMap;
  // Any hits recovered (will also be registered in map)
  std::vector< hitstruct >                                      hitsPossiblyRecovered;

  // TODO: I *think* the way this is being handled we shouldn't duplicate hits (e.g. if it is near the split region
  //       and has two associated wires), however I guess a "decision" on which of the locations the hit should be
  //       isn't made. Not sure if that has implications...

  // 1. Linear regression methods (a sort of "neighbors" method and input cluster set)
  // For all hits: use the hit->Channel->WireIDVector to get the possible matches
  // For reduced hits: since a choice of the hit location has already been made in the 3d match, here we just use the hit->WireID
  // TODO: is that right?

  std::map< geo::PlaneID, std::vector< std::pair<double,double> > > linregNeighbors;
  std::map< geo::PlaneID, std::vector< std::pair<double,double> > > linregClusters;

  // Get a list of linear regressions from the neighbors method if fUseFindNeighbors
  if ( fUseFindNeighbors ) {
    for ( auto const& iLabel : fReducHitsLabelVec ) {
      art::Handle< std::vector<recob::Hit> > hitsHandle;
      std::vector< art::Ptr<recob::Hit> > hits;
      if ( e.getByLabel(iLabel,hitsHandle) ) {
        art::fill_ptr_vector(hits,hitsHandle);
      }
      else{
        mf::LogWarning("GaussHitRecovery") << "Event failed to find recob::Hit with label " << iLabel;
        return;
      }

      for ( auto const& iHitPtr : hits ) {
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
        int nNeighbors = -1; // start at -1 since we will count ourself later
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
  
        linregNeighbors[ wID.asPlaneID() ].push_back( std::make_pair(slope,intercept) );
      } // loop hits
    } // loop hit labels
  } // end fUseFindNeighbors

  // if fUseReducHitClusters then use the hit assns to recob::Cluster
  if ( fUseReducHitClusters ) {
    for ( auto const& iLabel : fReducHitsLabelVec ) {
      // Hit set
      art::Handle< std::vector<recob::Hit> > hitsHandle;
      std::vector< art::Ptr<recob::Hit> > hits;
      if ( e.getByLabel(iLabel,hitsHandle) ) {
        art::fill_ptr_vector(hits,hitsHandle);
      }
      else{
        mf::LogWarning("GaussHitRecovery") << "Event failed to find recob::Hit with label " << iLabel;
        return;
      }

      // Cluster set
      art::Handle< std::vector<recob::Cluster> > clusterHandle;
      std::vector< art::Ptr<recob::Cluster> > clusters;
      if ( e.getByLabel(iLabel,clusterHandle) ) {
        art::fill_ptr_vector(clusters,clusterHandle);
      }
      else{
        mf::LogWarning("GaussHitRecovery") << "Event failed to find recob::Cluster with label " << iLabel;
        return;
      }

      // Find many for the clusters, hits
      art::FindManyP<recob::Hit> fmhit(clusterHandle, e, iLabel);
      if( !fmhit.isValid() ){
        mf::LogError("GaussHitRecovery") << "Error in validity of fmhit. Returning.";
        return;
      }

      for ( auto const& iCluster : clusters ) { // loop clusters
        std::vector< art::Ptr<recob::Hit> > clusterHits = fmhit.at(iCluster.key());

        if ( clusterHits.size()==0 ) continue;

        std::map< geo::PlaneID, std::vector<int> > wire_list;
        std::map< geo::PlaneID, std::vector<float> > tick_list;
        std::map< geo::PlaneID, std::vector<float> > rms_list;

        for( auto const& iHitPtr : clusterHits ) {
          geo::WireID const& wID = iHitPtr->WireID();
          if ( fSkipRecoveryPlane == (int)wID.Plane ) continue;

          tick_list[wID.asPlaneID()].push_back( iHitPtr->PeakTime() );
          rms_list[wID.asPlaneID()].push_back( iHitPtr->RMS() );
          wire_list[wID.asPlaneID()].push_back( wID.Wire );
        }

        for ( auto const& iPlaneID : fGeom->IteratePlaneIDs() ) {
          if ( tick_list.find(iPlaneID) == tick_list.end() ) continue;

          // Do a linear regression to see if this group of hits is "line-like" (aka track-like)
          // See especially:
          //   https://en.wikipedia.org/wiki/Simple_linear_regression
          //   https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
          // Doesn't account for hit RMS's here, just uses the peak times
          // -- Get m and b
          double xave = 0.;
          double yave = 0.;
          for( unsigned int idx=0; idx<wire_list[iPlaneID].size(); idx++ ){
            xave+=wire_list[iPlaneID][idx];
            yave+=tick_list[iPlaneID][idx];
          }
          xave/=double(wire_list[iPlaneID].size());
          yave/=double(wire_list[iPlaneID].size());
          // slope: numerator is s_{x,y}, denominator s_{x^2}
          double s_xy = 0.;
          double s_x2 = 0.;
          for( unsigned int idx=0; idx<wire_list[iPlaneID].size(); idx++ ){
            s_xy+=((wire_list[iPlaneID][idx]-xave)*(tick_list[iPlaneID][idx]-yave));
            s_x2+=((wire_list[iPlaneID][idx]-xave)*(wire_list[iPlaneID][idx]-xave));
          }
          double slope = s_xy / s_x2;
          double intercept = yave - slope*xave;
    
          // -- Get Chi2/N.D.F.
          // assumes sigma=hit RMS
          double chi2ndf = 0.;
          for( unsigned int idx=0; idx<wire_list[iPlaneID].size(); idx++ ){
            double expect=slope*wire_list[iPlaneID][idx]+intercept;
            chi2ndf+=((tick_list[iPlaneID][idx]-expect)*(tick_list[iPlaneID][idx]-expect)/(rms_list[iPlaneID][idx]*rms_list[iPlaneID][idx]));
          }
          chi2ndf/=(double(wire_list[iPlaneID].size()-2.));
    
          if( chi2ndf > fMaxChi2NDF ) continue;
    
          linregClusters[ iPlaneID ].push_back( std::make_pair(slope,intercept) );
        } // loop planeIDs
      } // loop clusters
    } // loop hit labels
  }

  // Print out how many entries we have per plane in the linreg saver:
  unsigned int totLinRegN = 0;
  unsigned int totLinRegC = 0;
  std::string printInfo = "Processing:\n ";
  for( geo::PlaneID const& thisPlaneID : fGeom->IteratePlaneIDs() ){
    printInfo += thisPlaneID.toString() + " linreg = ";
    if ( linregNeighbors.find(thisPlaneID)==linregNeighbors.end() ) {
      printInfo += "0 ";
    }
    else {
      printInfo += std::to_string(linregNeighbors[thisPlaneID].size()) + " ";
      totLinRegN += linregNeighbors[thisPlaneID].size();
    }
    if ( linregClusters.find(thisPlaneID)==linregClusters.end() ) {
      printInfo += "0";
    }
    else {
      printInfo += std::to_string(linregClusters[thisPlaneID].size());
      totLinRegC += linregClusters[thisPlaneID].size();
    }
    printInfo += "\n ";
  }
  printInfo += "Totals: " + std::to_string(totLinRegN) + " " + std::to_string(totLinRegC) + "\n";
  mf::LogInfo("GaussHitRecovery") << printInfo;

  // All hits, save some as the recovered hits. The LinReg counts should only be >0 if that method is "enabled"
  if ( totLinRegN > 0 || totLinRegC > 0 ){
    for ( auto const& iLabel : fHitsLabelVec ) {
      art::Handle< std::vector<recob::Hit> > hitsHandle;
      std::vector< art::Ptr<recob::Hit> > hits;
      if ( e.getByLabel(iLabel,hitsHandle) ) {
        art::fill_ptr_vector(hits,hitsHandle);
      }
      else{
        mf::LogWarning("GaussHitRecovery") << "Event failed to find recob::Hit with label " << iLabel;
        return;
      }

      // //art::Handle< std::vector<recob::Wire> > wiresHandle;
      // //if ( !e.getByLabel(fWireLabel,wiresHandle) ) {
      // //  mf::LogWarning("GaussHitRecovery") << "Event failed to find recob::Wire with label " << fWireLabel;
      // //	return;
      // //}
      // std::vector< art::Handle< std::vector<recob::Wire> > > wiresHandleVec;
      // wiresHandleVec.resize( fWireLabelVec.size() );
      // for ( unsigned int jLabel=0; jLabel<fWireLabelVec.size(); ++jLabel ) {
      //   if ( !e.getByLabel(fWireLabelVec[jLabel],wiresHandleVec[jLabel]) ) {
      //     mf::LogWarning("GaussHitRecovery") << "Event failed to find recob::Wire with label " << fWireLabelVec[jLabel];
      //     return;
      //   }
      // }

      art::FindManyP<recob::Wire> fmwire(hitsHandle, e, iLabel);
      if( !fmwire.isValid() ){
        mf::LogError("GaussHitRecovery") << "Error in validity of fmwire. Returning.";
        return;
      }

      // Loop through the set of all hits
      for ( auto const& iHitPtr : hits ) {
        auto const& thisTime = iHitPtr->PeakTime();
        auto const& thisRMS = iHitPtr->RMS();

        // Get the WireIDs matching to the Channel of this hit (can be >1, e.g. in the overlap region near middle)
        std::vector< geo::WireID > wires = fGeom->ChannelToWire( iHitPtr->Channel() );

        // Check if the hit is already in reduced set and/or if we're skipping this plane
        bool isReduc = false;
        bool isSkippedPlane = false;
        for ( auto const& iWire : wires ) {
          if ( fSkipRecoveryPlane == (int)iWire.Plane ) isSkippedPlane = true;
          if ( isReduc || isSkippedPlane ) break;
          if ( wireToReducHitsTimeRMS.find(iWire) == wireToReducHitsTimeRMS.end() ) continue;
          // TODO: is looping through like this and checking == on Tick and RMS the best way?
          for ( auto const& iTimeRMS : wireToReducHitsTimeRMS[iWire] ){
            if ( thisTime == iTimeRMS.first && thisRMS == iTimeRMS.second ) {
              isReduc = true;
              break;
            }
          }
        }
        if ( isReduc || isSkippedPlane ) continue; // if it's a reduced hit or on a plane we skip, don't try to save...

        for ( auto const& iWire : wires ) {
          auto const& thisWire = iWire.Wire;

          // loop through the linreg for this PlaneID and see if the hit is consistent with one of them
          // First the Clusters LinReg if any:
          if( linregClusters.find(iWire.asPlaneID()) != linregClusters.end() ) {
            for ( auto const& iLine : linregClusters[iWire.asPlaneID()] ) {
              double expect = iLine.first*thisWire + iLine.second;
              if( thisTime+fLinregTolerance*thisRMS > expect && thisTime-fLinregTolerance*thisRMS < expect ){
                // Put the hit into the recovered hit vector if it's not already registered
                if ( hitToRecoveryMethodMap.find( std::make_pair(iWire,thisTime) ) == hitToRecoveryMethodMap.end() ) {
                  std::vector< art::Ptr<recob::Wire> > hitWires = fmwire.at( iHitPtr.key() );
                  if ( hitWires.size() == 0 ) throw cet::exception("GaussHitRecovery") << "Hit found with no associated wires...\n";
                  else if ( hitWires.size() > 1 ) mf::LogWarning("GaussHitRecovery") << "Hit with >1 recob::Wire associated...";

                  hitstruct tmp{*iHitPtr,hitWires[0]};

                  hitsPossiblyRecovered.push_back( std::move(tmp) );
                  hitToRecoveryMethodMap[ std::make_pair(iWire,thisTime) ].push_back(1);
                }
                else {
                  // in case we already found this hit with this method
                  bool alreadyThisMethod = false;
                  for ( auto const& recoveryMethod : hitToRecoveryMethodMap[ std::make_pair(iWire,thisTime) ] ) {
                    if ( recoveryMethod==1 ) {
                      alreadyThisMethod = true;
                      break;
                    }
                  }
                  if ( !alreadyThisMethod ) hitToRecoveryMethodMap[ std::make_pair(iWire,thisTime) ].push_back(1);
                }
                break;
              }
            } // loop linear regressions (clusters) for this plane
          }

          // Second the Neighbors LinReg if any:
          if ( linregNeighbors.find(iWire.asPlaneID()) != linregNeighbors.end() ) {
            for ( auto const& iLine : linregNeighbors[iWire.asPlaneID()] ) {
              double expect = iLine.first*thisWire + iLine.second;
              if( thisTime+fLinregTolerance*thisRMS > expect && thisTime-fLinregTolerance*thisRMS < expect ){
                // Put the hit into the recovered hit vector if it's not already registered
                if ( hitToRecoveryMethodMap.find( std::make_pair(iWire,thisTime) ) == hitToRecoveryMethodMap.end() ) {
                  recob::Hit theHit = *iHitPtr;
                  std::vector< art::Ptr<recob::Wire> > hitWires = fmwire.at( iHitPtr.key() );
                  if ( hitWires.size() == 0 ) throw cet::exception("GaussHitRecovery") << "Hit found with no associated wires...\n";
                  else if ( hitWires.size() > 1 ) mf::LogWarning("GaussHitRecovery") << "Hit with >1 recob::Wire associated...";

                  hitstruct tmp{theHit,hitWires[0]};

                  hitsPossiblyRecovered.push_back( std::move(tmp) );
                  hitToRecoveryMethodMap[ std::make_pair(iWire,thisTime) ].push_back(0);
                }
                else {
                  // in case we already found this hit with this method
                  bool alreadyThisMethod = false;
                  for ( auto const& recoveryMethod : hitToRecoveryMethodMap[ std::make_pair(iWire,thisTime) ] ) {
                    if ( recoveryMethod==0 ) {
                      alreadyThisMethod = true;
                      break;
                    }
                  }
                  if ( !alreadyThisMethod ) hitToRecoveryMethodMap[ std::make_pair(iWire,thisTime) ].push_back(0);
                }
                break;
              }
            } // loop linear regressions (clusters) for this plane
          }
        } // loop wires matching to the channel in the hit
      } // loop hits
    } // loop hit labels
  } // only if total Linear Regressions saved is > 0


  // 2. Method using input PCAxis objects
  // Idea is to look for hits in the full set that are (between?) two clusters "pointing at" each other...
  // TODO: this method could probably be improved. Do we want to be using TVector3 or something else?

  std::map< geo::PlaneID, std::vector< std::pair<TVector3,TVector3> > >             allPcaVectors;
  std::map< geo::PlaneID, std::vector< std::pair<float,float> > >                   allPcaMagnitudes;
  // A pair of pair of floats -> wire,tick for pca1 (.first.{first,second}) & wire, tick for pca2 (.second.{first,second})
  // Then there are a vector of these per plane for the event.
  std::map< geo::PlaneID, std::vector< std::pair< std::pair<float,float>, std::pair<float,float> > > >    allPcaClusterPoint;

  // If fUseReducHitPCAxis, loop through the reduced hits and save some linear approx before checking all hits for matches
  if ( fUseReducHitPCAxis ) {
    for ( auto const& iLabel : fReducHitsLabelVec ) {
      art::Handle< std::vector<recob::Hit> > hitsHandle;
      std::vector< art::Ptr<recob::Hit> > hits;
      if ( e.getByLabel(iLabel,hitsHandle) ) {
        art::fill_ptr_vector(hits,hitsHandle);
      }
      else{
        mf::LogWarning("GaussHitRecovery") << "Event failed to find recob::Hit with label " << iLabel;
        return;
      }

      art::FindManyP<recob::SpacePoint> fmsps(hitsHandle, e, iLabel);
      if( !fmsps.isValid() ){
        mf::LogError("GaussHitRecovery") << "Error in validity of fmsps. Returning.";
        return;
      }

      art::Handle< std::vector<recob::PFParticle> > pfpsHandle;
      std::vector< art::Ptr<recob::PFParticle> > pfps;
      if ( e.getByLabel(iLabel,pfpsHandle) ) {
        art::fill_ptr_vector(pfps,pfpsHandle);
      }
      else{
        mf::LogWarning("GaussHitRecovery") << "Event failed to find recob::PFParticle with label " << iLabel;
        return;
      }

      art::FindManyP<recob::PCAxis> fmpca(pfpsHandle, e, iLabel);
      if( !fmpca.isValid() ){
        mf::LogError("GaussHitRecovery") << "Error in validity of fmpca. Returning.";
        return;
      }

      art::FindManyP<recob::Cluster> fmcl(pfpsHandle, e, iLabel);
      if( !fmcl.isValid() ){
        mf::LogError("GaussHitRecovery") << "Error in validity of fmcl. Returning.";
        return;
      }

      art::Handle< std::vector<recob::Cluster> > clusterHandle;
      std::vector< art::Ptr<recob::Cluster> > clusters;
      if ( e.getByLabel(iLabel,clusterHandle) ) {
        art::fill_ptr_vector(clusters,clusterHandle);
      }
      else{
        mf::LogWarning("GaussHitRecovery") << "Event failed to find recob::Cluster with label " << iLabel;
        return;
      }

      art::FindManyP<recob::Hit> fmhit(clusterHandle, e, iLabel);
      if( !fmhit.isValid() ){
        mf::LogError("GaussHitRecovery") << "Error in validity of fmhit. Returning.";
        return;
      }

      std::vector<TVector3> pcaVectors;
      std::vector<float> pcaMagnitudes;
      std::vector< std::map< geo::PlaneID, std::pair<float,float> > > pcaClusterPoint; // a point relating to cluster
      std::vector< std::vector< art::Ptr<recob::SpacePoint> > > pcaSpacePoint; // pick out the spacepoints relating to this point for comparison

      // Loop through PFParticles
      // Get the pairs of PCA axes pointing toward each other
      for ( auto const& iPFP : pfps ) {
        std::map< geo::PlaneID, std::pair<float,float> > thisPcaClusterPoint;

        std::vector< art::Ptr<recob::Cluster> > clst = fmcl.at(iPFP.key());
        //std::cout << "pfp matches " << clst.size() << " clusters" << std::endl;

        std::vector< art::Ptr<recob::PCAxis> > pcas = fmpca.at(iPFP.key());

        if ( pcas.size()==0 ) {
          mf::LogWarning("GaussHitRecovery") << "WARNING!! Looking at PFP without 2 PCAxis objects (it's 0!). Continuing to next PFP";
          continue;
        }
        else if ( pcas.size()!=2 ) {
          mf::LogWarning("GaussHitRecovery") << "WARNING!! Looking at PFP without 2 PCAxis objects. Curious!";
        }

        auto const& ax1 = pcas[0];

        // Principal axes:
        auto const& ax1vec = ax1->getEigenVectors();
        auto const& ax1val = ax1->getEigenValues();

        // Pick out if PCAxis suggests it is linear, at least order of magnitude larger than 2nd and 3rd components.
        if ( (ax1val[2] <= 10.*ax1val[1]) || (ax1val[2] <= 10.*ax1val[0]) ) continue;

        // Get a point in the PCAxis related cluster, map it plane by plane
        unsigned int useClstIdx=0;
        for ( auto const& iClst : clst ) {
          std::vector< art::Ptr<recob::Hit> > clstHits = fmhit.at( iClst.key() );
          if ( clstHits.size() == 0 ) {
            useClstIdx+=1;
            continue;
          }
          // consider just hit[0]
          if ( fSkipRecoveryPlane == (int)clstHits[0]->WireID().Plane ) continue;
          std::pair hitPoint = std::make_pair( clstHits[0]->WireID().Wire, clstHits[0]->PeakTime() );
          thisPcaClusterPoint[ clstHits[0]->WireID().asPlaneID() ] = hitPoint;
        }

        // if no clusters have hits, skip this pfp/pca...
        if ( useClstIdx==clst.size() ) continue;
        // get the spacepoint for the 0th hit:
        std::vector< art::Ptr<recob::SpacePoint> > useSps = fmsps.at ( (fmhit.at( clst[useClstIdx].key() ).at(0)).key() );
        // TODO: For now return all the spacepoints and see if any align when comparing potential "matches"
        //std::cout << "useSps.size() = " << useSps.size() << std::endl;
        if ( useSps.size()==0 ) continue;
        pcaSpacePoint.push_back( useSps );

        pcaClusterPoint.push_back(thisPcaClusterPoint);

        TVector3 pcaV(ax1vec[2].at(0),ax1vec[2].at(1),ax1vec[2].at(2));
        pcaVectors.push_back(pcaV);

        pcaMagnitudes.push_back( (float)ax1val[2] );
      }

      // Look at the PCAs:
      // The objects pcaClusterPoint, pcaVectors, pcaMagnitudes should all line up in terms of entries, so can be matched
      // Then, assign any matches to the overall map's vectors
      for ( unsigned int i=0; i<(pcaVectors.size()-1); ++i ) {
        auto const& pca1 = pcaVectors[i];
        // Check this against the other line-like PCAxis objects not already considered
        for ( unsigned int j=i+1; j<pcaVectors.size(); ++j ) {
          auto const& pca2 = pcaVectors[j];
          // Check if this point is within the tolerance of the PCAxis "cone"
          double vAngle = ( 180./TMath::Pi() )*pca1.Angle( pca2 );
          vAngle = std::min( vAngle, 180.-vAngle );
          if( vAngle > fPCAxisTolerance ) continue; // doesn't pass "pointing"/"angular" tolerance

          // And check if the line connecting spacepoints (any combination...) matches up with the PCAxes
          bool matchingSps = false;
          for ( auto const& spi : pcaSpacePoint[i] ) {
            for ( auto const& spj : pcaSpacePoint[j] ) {
              auto const& pointI = spi->position();
              auto const& pointJ = spj->position();
              TVector3 rPoints( pointJ.X()-pointI.X(), pointJ.Y()-pointI.Y(), pointJ.Z()-pointI.Z() );
              double rAngle = ( 180./TMath::Pi() )*pca1.Angle( rPoints );
              rAngle = std::min( rAngle, 180.-rAngle );
              if( rAngle < fPCAxisTolerance ) {
                matchingSps = true;
                break;
              }
            }
            if ( matchingSps ) break;
          }
          if ( !matchingSps ) continue;

          for ( auto const& iPlaneID : fGeom->IteratePlaneIDs() ) {
            if ( pcaClusterPoint.at(i).find(iPlaneID)==pcaClusterPoint.at(i).end() ||
                 pcaClusterPoint.at(j).find(iPlaneID)==pcaClusterPoint.at(j).end() )
              continue; // no points for cluster in this plane

            allPcaVectors[ iPlaneID ].push_back( std::make_pair(pca1,pca2) ); // store the points and such for use later
            allPcaMagnitudes[ iPlaneID ].push_back( std::make_pair(pcaMagnitudes[i],pcaMagnitudes[j]) );

            //auto wire1 = pcaClusterPoint.at(i)[iPlaneID].first;
            //auto tick1 = pcaClusterPoint.at(i)[iPlaneID].second;
            //auto wire2 = pcaClusterPoint.at(j)[iPlaneID].first;
            //auto tick2 = pcaClusterPoint.at(j)[iPlaneID].second;
            allPcaClusterPoint[ iPlaneID ].push_back( std::make_pair( pcaClusterPoint.at(i)[iPlaneID], pcaClusterPoint.at(j)[iPlaneID] ) );
          } // plane for pca1 pca2 combo
        } // end pca2 loop
      } // end pca1 loop
    } // end loop reduced hit labels

    // Print out how many entries we have per plane in the pcaxis info saver:
    std::string printInfo = "Processing:\n ";
    for( geo::PlaneID const& thisPlaneID : fGeom->IteratePlaneIDs() ){
      printInfo += thisPlaneID.toString() + " pcaxis groupings = ";
      if ( allPcaClusterPoint.find(thisPlaneID)==allPcaClusterPoint.end() )
        printInfo += "0 ";
      else
        printInfo += std::to_string(allPcaClusterPoint[ thisPlaneID ].size()) + " ";
      printInfo += "\n ";
    }
    mf::LogInfo("GaussHitRecovery") << printInfo;

    // Now loop through all hits and see if we want to recover them:
    for ( auto const& iLabel : fHitsLabelVec ) {
      art::Handle< std::vector<recob::Hit> > hitsHandle;
      std::vector< art::Ptr<recob::Hit> > hits;
      if ( e.getByLabel(iLabel,hitsHandle) ) {
        art::fill_ptr_vector(hits,hitsHandle);
      }
      else{
        mf::LogWarning("GaussHitRecovery") << "Event failed to find recob::Hit with label " << iLabel;
        return;
      }

      // std::vector< art::Handle< std::vector<recob::Wire> > > wiresHandleVec;
      // wiresHandleVec.resize( fWireLabelVec.size() );
      // for ( unsigned int jLabel=0; jLabel<fWireLabelVec.size(); ++jLabel ) {
      //   if ( !e.getByLabel(fWireLabelVec[jLabel],wiresHandleVec[jLabel]) ) {
      //     mf::LogWarning("GaussHitRecovery") << "Event failed to find recob::Wire with label " << fWireLabelVec[jLabel];
      //     return;
      //   }
      // }

      art::FindManyP<recob::Wire> fmwire(hitsHandle, e, iLabel);
      if( !fmwire.isValid() ){
        mf::LogError("GaussHitRecovery") << "Error in validity of fmwire. Returning.";
        return;
      }

      for ( auto const& iHitPtr : hits ) {
        auto const& thisTime = iHitPtr->PeakTime();
        auto const& thisRMS = iHitPtr->RMS();
        std::vector< geo::WireID > wires = fGeom->ChannelToWire( iHitPtr->Channel() );

        // Check if the hit is already in reduced set and/or if we're skipping this plane
        bool isReduc = false;
        bool isSkippedPlane = false;
        for ( auto const& iWire : wires ) {
          if ( fSkipRecoveryPlane == (int)iWire.Plane ) isSkippedPlane = true;
          if ( isReduc || isSkippedPlane ) break;
          if ( wireToReducHitsTimeRMS.find(iWire) == wireToReducHitsTimeRMS.end() ) continue;
          // TODO: is looping through like this and checking == on Tick and RMS the best way?
          for ( auto const& iTimeRMS : wireToReducHitsTimeRMS[iWire] ){
            if ( thisTime == iTimeRMS.first && thisRMS == iTimeRMS.second ) {
              isReduc = true;
              break;
            }
          }
        }
        if ( isReduc || isSkippedPlane ) continue; // if it's a reduced hit or on a plane we skip, don't try to save...

        // Consider both wire possibilities for hits in the overlap regions
        for ( auto const& iWire : wires ) {
          auto const& thisWire = iWire.Wire;
          auto thisPlaneID = iWire.asPlaneID();

          // Loop through the wire,tick point combinations and figure out 
          if ( allPcaClusterPoint.find( thisPlaneID ) == allPcaClusterPoint.end() ) continue;

          for ( unsigned int iPtPair = 0; iPtPair < allPcaClusterPoint[ thisPlaneID ].size(); ++iPtPair ) {
            bool firstIsLoWire = (allPcaClusterPoint[ thisPlaneID ].at(iPtPair).first.first <
                                  allPcaClusterPoint[ thisPlaneID ].at(iPtPair).second.first) ? true : false;
            float loWire     = firstIsLoWire ? allPcaClusterPoint[ thisPlaneID ].at(iPtPair).first.first : 
                                               allPcaClusterPoint[ thisPlaneID ].at(iPtPair).second.first;
            float loWireTick = firstIsLoWire ? allPcaClusterPoint[ thisPlaneID ].at(iPtPair).first.second : 
                                               allPcaClusterPoint[ thisPlaneID ].at(iPtPair).second.second;
            float hiWire     = firstIsLoWire ? allPcaClusterPoint[ thisPlaneID ].at(iPtPair).second.first :
                                               allPcaClusterPoint[ thisPlaneID ].at(iPtPair).first.first;
            float hiWireTick = firstIsLoWire ? allPcaClusterPoint[ thisPlaneID ].at(iPtPair).second.second :
                                               allPcaClusterPoint[ thisPlaneID ].at(iPtPair).first.second;

            // if we only want to use hits between 2 groups
            if ( fPCAxisOnlyBetween && (thisWire < loWire || thisWire > hiWire) ) continue;
            if ( (int)hiWire == (int)loWire ) continue;

            // make a simple line connecting these two:
            // m = (y2-y1) / (x2 - x1)
            float connectSlope = (hiWireTick - loWireTick) / (hiWire - loWire);
            // b = y1 - m*x1
            float connectInt = hiWireTick - (connectSlope * hiWire);
            // Is our hit within a tolerance of this line (fLinregTolerance)
            float tickExpect = (connectSlope * thisWire) + connectInt;
            if ( std::fabs(tickExpect-thisTime) > fLinregTolerance*thisRMS ) continue;

            // Put the hit into the recovered hit vector if it's not already registered
            if ( hitToRecoveryMethodMap.find( std::make_pair(iWire,thisTime) ) == hitToRecoveryMethodMap.end() ) {
              std::vector< art::Ptr<recob::Wire> > hitWires = fmwire.at( iHitPtr.key() );
              if ( hitWires.size() == 0 ) throw cet::exception("GaussHitRecovery") << "Hit found with no associated wires...\n";
              else if ( hitWires.size() > 1 ) mf::LogWarning("GaussHitRecovery") << "Hit with >1 recob::Wire associated...";

              hitstruct tmp{*iHitPtr,hitWires[0]};

              hitsPossiblyRecovered.push_back( std::move(tmp) );
              hitToRecoveryMethodMap[ std::make_pair(iWire,thisTime) ].push_back(2);
            }
            else {
              // in case we already found this hit with this method
              bool alreadyThisMethod = false;
              for ( auto const& recoveryMethod : hitToRecoveryMethodMap[ std::make_pair(iWire,thisTime) ] ) {
                if ( recoveryMethod==2 ) {
                  alreadyThisMethod = true;
                  break;
                }
              }
              if ( !alreadyThisMethod ) hitToRecoveryMethodMap[ std::make_pair(iWire,thisTime) ].push_back(2);
            }
            break;
          } // end loop of point-pairs we're testing the hit against
        } // end loop of wires for the hit
      } // end loop all hits
    } // end loop all hit labels

  } // end if fUseReducHitPCAxis


  /////////////////////////////////////////////////////////////////////
  //
  // STEP 3: Now we have all the recovery attempts, let's see what we
  //         should put back into hit set and then put that into the evt
  //
  /////////////////////////////////////////////////////////////////////

  for ( auto const& iHitStruct : hitsPossiblyRecovered ) {
    auto const& thisWireID = iHitStruct.stHit.WireID();
    float thisWireTick = iHitStruct.stHit.PeakTime();
    // Require at least fMinMethodsSuccess methods to have found the hit for it to be recovered. Allows to use
    // a combination of algorithms
    if ( fMinMethodsSuccess>=1 && 
         hitToRecoveryMethodMap.find( std::make_pair(thisWireID,thisWireTick) )==hitToRecoveryMethodMap.end() )
      continue;
    if ( fMinMethodsSuccess > (int)hitToRecoveryMethodMap[ std::make_pair(thisWireID,thisWireTick) ].size() )
      continue;
    hitCol.emplace_back(iHitStruct.stHit,iHitStruct.stWire);
  }

  hitCol.put_into(e);
}

DEFINE_ART_MODULE(GaussHitRecovery)
