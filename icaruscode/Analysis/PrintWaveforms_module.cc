////////////////////////////////////////////////////////////////////////
// Class:       PrintWaveforms
// Plugin Type: analyzer (Unknown Unknown)
// File:        PrintWaveforms_module.cc
//
// Generated at Fri Sep 23 17:04:48 2022 by Bruce Howard using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "cetlib_except/exception.h"

#include "lardataobj/RawData/RawDigit.h"

#include "TRandom3.h"

#include <iostream>
#include <memory>
#include <string>

class PrintWaveforms;

class PrintWaveforms : public art::EDAnalyzer {
public:
  explicit PrintWaveforms(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PrintWaveforms(PrintWaveforms const&) = delete;
  PrintWaveforms(PrintWaveforms&&) = delete;
  PrintWaveforms& operator=(PrintWaveforms const&) = delete;
  PrintWaveforms& operator=(PrintWaveforms&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  std::vector< art::InputTag > fTPCRawDigitVec;
  std::vector< int > fListCryostat;
  std::vector< int > fListTPC;
  std::vector< int > fListPlane;
  std::vector< int > fListWire;
  std::vector< int > fListTick;
  std::vector< int > fListLabel; // the by-hand classification 0:background, 1:signal

  std::unique_ptr< TRandom3 > rand;
  geo::GeometryCore const* fGeom;
};


PrintWaveforms::PrintWaveforms(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fTPCRawDigitVec ( p.get< std::vector< art::InputTag > >("RawDigitLabelVec",{}) ),
  fListCryostat   ( p.get< std::vector< int > >("ListCryostat",{}) ),
  fListTPC        ( p.get< std::vector< int > >("ListTPC",{}) ),
  fListPlane      ( p.get< std::vector< int > >("ListPlane",{}) ),
  fListWire       ( p.get< std::vector< int > >("ListWire",{}) ),
  fListTick       ( p.get< std::vector< int > >("ListTick",{}) ),
  fListLabel      ( p.get< std::vector< int > >("ListLabel",{}) )
  // More initializers here.
{

  // Do checks
  if ( fListCryostat.size() != fListTPC.size() ||
       fListCryostat.size() != fListPlane.size() ||
       fListCryostat.size() != fListWire.size() ||
       fListCryostat.size() != fListTick.size() || 
       fListCryostat.size() != fListLabel.size() ) {
    std::cout << "ERROR." << std::endl;
  }

  // initialize a TRandom3
  rand = std::make_unique< TRandom3 >();

}

void PrintWaveforms::analyze(art::Event const& e)
{
  // PRINT WAVEFORMS...
  fGeom = lar::providerFrom<geo::Geometry>();

  std::map < geo::WireID, std::vector< std::pair< int, std::vector<int> > > > wireidTickWvfmVecMap;
  std::map < geo::WireID, std::vector<int> >                                  wireidLabelVecMap;

  for ( unsigned int idxWvfm=0; idxWvfm<fListCryostat.size(); ++idxWvfm ) {
    geo::WireID wire( fListCryostat[idxWvfm], fListTPC[idxWvfm], fListPlane[idxWvfm], fListWire[idxWvfm] );
    if ( wireidTickWvfmVecMap.find( wire ) != wireidTickWvfmVecMap.end() ) {
      std::cout << "WARNING... Already have this WireID. Were you expecting this?" << std::endl;
      wireidTickWvfmVecMap[ wire ].push_back( std::make_pair( fListTick[idxWvfm], std::vector<int>() ) );
      wireidLabelVecMap [ wire ].push_back( fListLabel[idxWvfm] );
    }
    else {
      wireidTickWvfmVecMap[ wire ] = { std::make_pair( fListTick[idxWvfm], std::vector<int>() ) };
      wireidLabelVecMap [ wire ] = { fListLabel[idxWvfm] };
    }
  }

  for ( auto const &[wire, tickWvfmVec] : wireidTickWvfmVecMap ) {
    for ( auto const tickWvfmPair : tickWvfmVec ) {
      std::cout << wire << " " << tickWvfmPair.first << std::endl;
    }
  }

  for ( auto const& label : fTPCRawDigitVec ) {
    std::cout << std::endl;
    std::cout << "Running on " << label << std::endl;
    art::Handle< std::vector< raw::RawDigit > > digitsHandle;
    std::vector< art::Ptr< raw::RawDigit > > digits;
    if ( e.getByLabel(label, digitsHandle) ) {
      art::fill_ptr_vector( digits, digitsHandle );
    }
    else {
      mf::LogWarning( "PrintWaveforms" ) << "Event failed to find raw::RawDigit with label " << label << ".";
      return;
    }

    std::cout << "  size = " << digits.size() << std::endl;

    for ( auto const& digit : digits ) {
      //std::cout << ".";
      auto chID = digit->Channel();
      std::vector< geo::WireID > wires = fGeom->ChannelToWire( chID );

      raw::RawDigit::ADCvector_t adcVals = digit->ADCs();

      for ( auto const& wire : wires ) {
	if ( wireidTickWvfmVecMap.find( wire ) != wireidTickWvfmVecMap.end() ) {
	  for ( unsigned int idxROI = 0; idxROI < wireidTickWvfmVecMap[ wire ].size(); ++idxROI ) {
	    int possibleTickStart = int( wireidTickWvfmVecMap[ wire ].at(idxROI).first - 100 ) + int( (rand->Rndm()*150.)-75. );
	    if ( possibleTickStart + 200 > 4095 ) possibleTickStart = 4095-200;
	    unsigned int tickStart = possibleTickStart < 0 ? 0 : (unsigned int)possibleTickStart;
	    std::cout << "    tick start = " << tickStart << std::endl;

	    for ( unsigned int idxTick=0; idxTick<200; ++idxTick ) {
	      wireidTickWvfmVecMap[ wire ].at(idxROI).second.push_back( adcVals[idxTick + tickStart] );
	    } // ADCs
	  } // ROIs on wire of interest
	} // is wire of interest
      } // WireIDs 
    } // RawDigits
  } // RawDigit Labels

  // Print waveforms in useful format
  std::string strInfo = "np.array( [";
  std::string strSgBk = "np.array( [";
  std::string strWvfm = "np.array( [";

  for ( auto const &[wire, tickWvfmVec] : wireidTickWvfmVecMap ) {
    for ( unsigned int idxROI = 0; idxROI < tickWvfmVec.size(); ++idxROI ) {
      strInfo+="["+std::to_string(wire.Cryostat)+", "+std::to_string(wire.TPC)+", "+std::to_string(wire.Plane)+", "+std::to_string(wire.Wire)+"], ";
      strSgBk+=std::to_string(wireidLabelVecMap[ wire ].at(idxROI))+", ";
      std::vector<int> adcVals = tickWvfmVec.at(idxROI).second;
      strWvfm+="[";
      for ( unsigned int idxTick = 0; idxTick < adcVals.size(); ++idxTick ) {
	strWvfm+=std::to_string(adcVals[idxTick]);
	if ( idxTick < adcVals.size()-1 ) strWvfm+=", ";
      }
      strWvfm+="], ";
    }
  }

  strInfo+="] )";
  strSgBk+="] )";
  strWvfm+="] )";

  std::cout << "RESULTS" << std::endl;
  std::cout << std::endl;
  std::cout << strInfo << std::endl;
  std::cout << std::endl;
  std::cout << strSgBk << std::endl;
  std::cout << std::endl;
  std::cout << strWvfm << std::endl;

}

DEFINE_ART_MODULE(PrintWaveforms)
