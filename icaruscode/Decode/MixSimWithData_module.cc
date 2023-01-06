// Bruce Howard
// Based very heavily on MixSimEvents_module.cc from NOvA,
//   see e.g. /cvmfs/nova.opensciencegrid.org/novasoft/slf7/novasoft/releases/S22-08-10/NovaSimMixer/MixSimEvents_module.cc

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardataobj/RawData/OpDetWaveform.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/PtrRemapper.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"
#include "art/Framework/Modules/MixFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/RootIOPolicy.h"
#include "art/Persistency/Common/CollectionUtilities.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Utilities/ScheduleID.h"

#include "fhiclcpp/ParameterSet.h"

#include <memory>
#include <string>
#include <vector>

namespace mix {
  class MixSimWithData {
  public:
    explicit MixSimWithData ( const fhicl::ParameterSet& p,
			      art::MixHelper& helper );

    // N secondaries per primary
    size_t nSecondaries();

    // Mixing functions:
    typedef std::vector<raw::RawDigit> RawDigitVector;
    bool mixRawDigit ( const std::vector<const RawDigitVector*>& in,
		       RawDigitVector& out,
		       const art::PtrRemapper& prm );

    typedef std::vector<recob::Hit> HitVector;
    bool mixHit ( const std::vector<const HitVector*>& in,
		  HitVector& out,
		  const art::PtrRemapper& prm );

    typedef std::vector<raw::OpDetWaveform> OpDetWaveformVector;
    bool mixOpDetWaveform ( const std::vector<const OpDetWaveformVector*>& in,
			    OpDetWaveformVector& out,
			    const art::PtrRemapper& prm );

  private:
    art::InputTag fRawDigitProducerLabel;
    art::InputTag fHitProducerLabel;
    art::InputTag fOpDetWaveformProducerLabel;

    bool          fMixRawDigits;
    bool          fMixHits;
    bool          fMixOpDetWaveforms;

    // Per the NOvA module, these are "returned from flattening the GenParticle collections."
    // (Also "These are needed in Ptr remapping operations")
    std::vector<size_t> fGenOffsets;
  };

  MixSimWithData::MixSimWithData ( fhicl::ParameterSet const &p,
				   art::MixHelper &helper ):
    fRawDigitProducerLabel      ( p.get< art::InputTag >("RawDigitProducerLabel","") ),
    fHitProducerLabel           ( p.get< art::InputTag >("HitProducerLabel","") ),
    fOpDetWaveformProducerLabel ( p.get< art::InputTag >("OpDetWaveformProducerLabel","") ),
    fMixRawDigits               ( p.get< bool >("MixRawDigits",false) ),
    fMixHits                    ( p.get< bool >("MixHits",false) ),
    fMixOpDetWaveforms          ( p.get< bool >("MixOpDetWaveforms",false) )
  {
    if ( !fMixRawDigits && !fMixHits && !fMixOpDetWaveforms )
      throw cet::exception("MixSimWithData") << "Error... need to be doing SOME mixing." << std::endl;

    // Registering mix-ops (same order as declared -- following the NOvA style)
    if ( fMixRawDigits ) {
      helper.declareMixOp ( fRawDigitProducerLabel, &MixSimWithData::mixRawDigit, *this );
    }
    if ( fMixHits ) {
      helper.declareMixOp ( fHitProducerLabel, &MixSimWithData::mixHit, *this );
    }
    if ( fMixOpDetWaveforms ) {
      helper.declareMixOp ( fOpDetWaveformProducerLabel, &MixSimWithData::mixOpDetWaveform, *this );
    }
  }

  size_t MixSimWithData::nSecondaries()
  {
    return 1; // need to understand if this needs to change
  }

  bool MixSimWithData::mixRawDigit ( const std::vector<const MixSimWithData::RawDigitVector*>& in,
				     MixSimWithData::RawDigitVector& out,
				     const art::PtrRemapper& prm )
  {
    // Just flattening to one product
    art::flattenCollections(in, out, fGenOffsets);
    return true;
  }

  bool MixSimWithData::mixHit ( const std::vector<const MixSimWithData::HitVector*>& in,
				MixSimWithData::HitVector& out,
				const art::PtrRemapper& prm )
  {
    // Just flattening to one product
    art::flattenCollections(in, out, fGenOffsets);
    return true;
  }

  bool MixSimWithData::mixOpDetWaveform ( const std::vector<const MixSimWithData::OpDetWaveformVector*>& in,
					  MixSimWithData::OpDetWaveformVector& out,
					  const art::PtrRemapper& prm )
  {
    // Just flattening to one product
    art::flattenCollections(in, out, fGenOffsets);
    return true;
  }

  using ModuleType = art::MixFilter<MixSimWithData, art::RootIOPolicy>;
  DEFINE_ART_MODULE(ModuleType)
} // namespace mix
