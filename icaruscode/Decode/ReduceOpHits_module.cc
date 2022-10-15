////////////////////////////////////////////////////////////////////////
// Class:       ReduceOpHits
// Plugin Type: producer (Unknown Unknown)
// File:        ReduceOpHits_module.cc
//
// Generated at Thu Oct 13 23:24:19 2022 by Bruce Howard using cetskelgen
// from  version .
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
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"

#include <memory>
#include <string>
#include <map>

class ReduceOpHits;


class ReduceOpHits : public art::EDProducer {
public:
  explicit ReduceOpHits(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ReduceOpHits(ReduceOpHits const&) = delete;
  ReduceOpHits(ReduceOpHits&&) = delete;
  ReduceOpHits& operator=(ReduceOpHits const&) = delete;
  ReduceOpHits& operator=(ReduceOpHits&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.

  bool fCheatMode; // In a file with waveforms fully saved (not just during the beam window) consider only waveforms with >10.5k samples to try to use just the beam waveforms (todo: better way?)
  art::InputTag fPMTWaveInputLabel;
  art::InputTag fPMTHitInputLabel;
};


ReduceOpHits::ReduceOpHits(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fCheatMode         ( p.get< bool >("CheatMode", false) ),
  fPMTWaveInputLabel ( p.get< art::InputTag >("PMTWaveInputLabel", "") ),
  fPMTHitInputLabel  ( p.get< art::InputTag >("PMTHitInputLabel", "") )
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  produces< std::vector<recob::OpHit> >();
}

void ReduceOpHits::produce(art::Event& e)
{
  // Implementation of required member function here.

  double deltaT = 0.002; // us -- todo: better way to automate PMT deltaT (is it in DetProp?)

  // Make a map of time ranges covered by a given channel:
  std::map < int, std::vector< std::pair<double,double> > > channelWaveformTimesMap;

  art::Handle< std::vector<raw::OpDetWaveform> > waveformsHandle;
  std::vector< art::Ptr<raw::OpDetWaveform> > waveforms;
  if ( e.getByLabel(fPMTWaveInputLabel,waveformsHandle) ) {
    art::fill_ptr_vector(waveforms,waveformsHandle);
  }
  else{
    mf::LogWarning("ReduceOpHits") << "Event failed to find raw::OpDetWaveform with label " << fPMTWaveInputLabel << ".";
    return;
  }

  for ( auto const& pdWaveform : waveforms ) {
    const raw::Channel_t chID = pdWaveform->ChannelNumber();
    const raw::TimeStamp_t time = pdWaveform->TimeStamp();
    const unsigned int wvfmLen = pdWaveform->Waveform().size();

    if ( fCheatMode && wvfmLen < 10500 ) continue;

    if ( channelWaveformTimesMap.find(chID) == channelWaveformTimesMap.end() ) channelWaveformTimesMap[chID] = {};
    channelWaveformTimesMap[ chID ].push_back( std::make_pair( time, time+(double(wvfmLen)*deltaT) ) );
  } // loop waveforms

  // Now loop through hits and remove the ones within the time ranges in the map:
  std::unique_ptr< std::vector< recob::OpHit > > opHitVec = std::make_unique< std::vector< recob::OpHit > >();

  art::Handle< std::vector<recob::OpHit> > hitsHandle;
  std::vector< art::Ptr<recob::OpHit> > hits;
  if ( e.getByLabel(fPMTHitInputLabel,hitsHandle) ) {
    art::fill_ptr_vector(hits,hitsHandle);
  }
  else{
    mf::LogWarning("ReduceOpHits") << "Event failed to find recob::OpHit with label " << fPMTHitInputLabel << ".";
    return;
  }

  for ( auto const& pdHit : hits ) {
    const int chID = pdHit->OpChannel();
    //const double pdTime = pdHit->PeakTime();
    const double pdTimeAbs = pdHit->PeakTimeAbs();

    bool shouldSaveOpHit = true;
    if ( channelWaveformTimesMap.find( chID ) != channelWaveformTimesMap.end() ) {
      for ( auto const& timePair : channelWaveformTimesMap[chID] ) {
	//std::cout << "Hit Times " << pdTime << ", abs = " << pdTimeAbs << ", Waveform Ranges " << timePair.first << " to " << timePair.second << std::endl;
	if ( pdTimeAbs >= timePair.first && pdTimeAbs <= timePair.second ) {
	  shouldSaveOpHit = false;
	  break;
	}
      }
    }

    if ( shouldSaveOpHit ) {
      recob::OpHit saveOpHit( pdHit->OpChannel(), pdHit->PeakTime(), pdHit->PeakTimeAbs(), pdHit->Frame(),
			      pdHit->Width(), pdHit->Area(), pdHit->Amplitude(), pdHit->PE(), pdHit->FastToTotal() );
      opHitVec->push_back( saveOpHit );
    }
  }

  e.put( std::move(opHitVec) );
}

DEFINE_ART_MODULE(ReduceOpHits)
