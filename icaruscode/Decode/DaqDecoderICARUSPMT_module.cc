/**
 * @file    DaqDecoderICARUSPMT_module.cc
 * @brief   Produces `raw::OpDetWaveform` from V1730 artDAQ data fragments.
 * @authors Andrea Scarpelli, Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date    May 21, 2021
 * 
 */


// ICARUS/SBN libraries
#include "icaruscode/Decode/DecoderTools/details/PMTDecoderUtils.h"
#include "icaruscode/Decode/ChannelMapping/IICARUSChannelMap.h"
#include "icarusalg/Utilities/FHiCLutils.h" // util::fhicl::getOptionalValue()
#include "icarusalg/Utilities/BinaryDumpUtils.h" // icarus::ns::util::bin()

#include "sbnobj/Common/PMT/Data/PMTconfiguration.h" // sbn::PMTconfiguration
#include "sbndaq-artdaq-core/Overlays/Common/CAENV1730Fragment.hh"

// LArSoft libraries
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataalg/DetectorInfo/DetectorTimings.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/Utilities/quantities/spacetime.h" // nanoseconds
#include "lardataalg/Utilities/intervals_fhicl.h" // for nanoseconds in FHiCL
#include "larcorealg/CoreUtils/enumerate.h"
#include "larcorealg/CoreUtils/counter.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/ExternalTrigger.h"

// artDAQ
#include "artdaq-core/Data/Fragment.hh"

// framework libraries
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/types/TableAs.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h"

// ROOT libraries
#include "TTree.h"

// C/C++ standard libraries
#include <memory>
#include <ostream>
#include <vector>
#include <string>
#include <optional>
#include <cassert>


//------------------------------------------------------------------------------
using namespace util::quantities::time_literals;

//------------------------------------------------------------------------------
namespace icarus { class DaqDecoderICARUSPMT; }
/**
 * @brief Produces `raw::OpDetWaveform` from V1730 artDAQ data fragments.
 * 
 * The module can read fragments from CAEN V1730 readout boards delivered by
 * artDAQ.
 * 
 * This decoder must support both a off-line mode (for storage and downstream
 * processing) and a on-line mode (for monitoring).
 * In particular, the on-line workflow is such that it may not be possible to
 * access the FHiCL configuration of the job and therefore the PMT configuration
 * data (see `icarus::PMTconfigurationExtraction` module).
 * 
 * 
 * Configuration
 * --------------
 * 
 * The set of supported parameters can be seen on command line by running
 * `lar --print-description DaqDecoderICARUSPMT`.
 * 
 * Description of the configuration parameters:
 * * `DiagnosticOutput` (flag, default: `false`): enables additional console
 *     output, including dumping of the fragments (that is huge output).
 * * `PMTconfigTag` (data product tag, optional): if specified, the pre-trigger
 *     buffer duration is read from there; although optional, it is strongly
 *     recommended that this information be provided, since it is essential for
 *     the correct timing of the PMT waveforms (see
 *     @ref icarus_PMTDecoder_timestamps "the discussion on time stamps below").
 * * `BoardSetup` (list of board setup information): each entry specifies some
 *     information about a specific readout board; the boards are identified by
 *     their name; if a board is found in input that has no setup information,
 *     some time corrections are not applied (see
 *     @ref icarus_PMTDecoder_timestamps "the discussion on time stamps below").
 *     Each entry is in the form of a table:
 *     * `Name` (string, mandatory): the name of the board
 *        (e.g. `"icaruspmtwwtop01"`); this is use to match the setup
 *        information to a fragment ID in the PMT configuration.
 *     * `FragmentID` (integral, optional): if specified, allows the corrections
 *       using setup information to be applied even when no PMT configuration is
 *       provided (if neither PMT configuration nor setup information including
 *       `FragmentID` is available, no time correction is applied).
 *     * `TriggerDelay` (nanoseconds, default: 0 ns): measured delay from the
 *       primitive trigger time to the execution of the PMT trigger; specify
 *       the unit! (e.g. `"43 ns"`).
 * * `RequireKnownBoards` (flag, default: `true`): if set, the readout boards
 *     in input must each have a setup configuration (`BoardSetup`) *and* must
 *     be present in the PMT DAQ configuration, or an exception is thrown;
 *     if not set, readout boards in input will be processed even when their
 *     setup or DAQ configuration is not known, at the cost of fewer corrections
 *     on the timestamps (which should then be considered unreliable).
 * * `RequireBoardConfig` (flag, default: `true`): if set, the readout boards
 *     which have a setup (`BoardSetup`) are required to be included in the DAQ
 *     configuration of the input file, or an exception is thrown; if not set,
 *     missing readout boards are unnoticed.
 * * `TriggerTag` (data product tag, mandatory): tag for the information
 *     (currently required to be a collection of `raw::ExternalTrigger`,
 *     in the future it should become `raw::Trigger`);
 * * `TTTresetEverySecond` (optional): if set, the decoder will take advantage
 *     of the assumption that the Trigger Time Tag of all PMT readout boards is
 *     synchronised with the global trigger time and reset at every change of
 *     second of the timescale of the latter; this is currently the only
 *     implementation supporting multiple PMT readout windows on the same board;
 *     if this option is set to `false`, all PMT readout boards are assumed to
 *     have been triggered at the time of the global trigger. By default, this
 *     option is set to `true` unless `TriggerTag` is specified empty.
 * * `DataTrees` (list of strings, default: none): list of data trees to be
 *     produced; if none (default), then `TFileService` is not required.
 * * `LogCategory` (string, default: `DaqDecoderICARUSPMT`): name of the message
 *     facility category where the output is sent.
 * 
 * 
 * Requirements
 * -------------
 * 
 * Services required include:
 * 
 * * `IICARUSChannelMap` for the association of fragments to LArSoft channel ID;
 * * `DetectorClocksService` for the correct decoding of the time stamps
 *   (always required, even when dumbed-down timestamp decoding is requested);
 * * `TFileService` only if the production of trees or plots is requested.
 * 
 * 
 * Waveform time stamp
 * --------------------
 * 
 * @anchor icarus_PMTDecoder_timestamps
 * 
 * All waveforms on the same readout board share the same timestamp.
 * 
 * The time stamp of the waveform is defined as the time when the first sample
 * of the waveform started (that is, if the sample represent the value of the
 * signal in an interval of 2 ns, the time stamp is pointing at the beginning
 * of those 2 ns). Whether we can honour that definition, though, is a different
 * matter.
 * The representation of the time stamp is in the
 * @ref DetectorClocksElectronicsTime "electronics time scale".
 * 
 * There are two "types" of waveforms: the ones acquired at global trigger time,
 * and the ones acquired because of a "local trigger" which was not promoted to
 * global (likely because not in coincidence with the beam gate).
 * In both cases, it is the same type of signal: a trigger primitive from
 * the NI7820 FPGA, which initializes the acquisition of the waveform.
 * Every delay between when that signal is emitted and when the PMT trigger is
 * executed shifts the time stamp of the waveform backward.
 * 
 * We assign the the time stamp of the waveforms matching the global trigger
 * as follow:
 * * the base time is the global trigger time; this effectively defines the
 *   electronics time scale, so its representation is a fixed number that is
 *   configured in LArSoft and can be accessed with
 *   `DetectorClocksData::TriggerTime()`;
 * * the delay of the propagation from the trigger board to the readout board
 *   is subtracted to the timestamp; this value must be independently measured
 *   and provided to this decoder via tool configuration as setup information
 *   (`TriggerDelay`); if not present in the setup, this delay is not added;
 * * upon receiving the trigger, the readout board will keep some of the samples
 *   already digitized, in what we call pre-trigger buffer; the size of this
 *   buffer is a fixed number of samples which is specified in DAQ as a fraction
 *   of the complete buffer that is _post-trigger_; this amount, converted in
 *   time, is subtracted to the trigger time to point back to the beginning of
 *   the waveform instead that to the trigger primitive time. The necessary
 *   information is read from the PMT configuration (`PMTconfigTag`); if no
 *   configuration is available, this offset is not subtracted; note that this
 *   is a major shift (typically, a few microseconds) that should be always
 *   included.
 * 
 * Each V1730 event record includes a trigger time tag (TTT), which is the value
 * of an internal counter of the board at the time the board received a trigger.
 * This can be used to relate the various waveforms (and the various fragments)
 * in the _art_ event.
 * 
 * 
 * Data trees
 * -----------
 * 
 * The module supports the following ROOT trees production on demand:
 * 
 * * `PMTfragments`: data pertaining a single fragment; each entry is about a
 *   single fragment, and it includes full event ID, event time stamp (from
 *   _art_, i.e. the one assigned to the event by artDAQ), the ID of the
 *   fragment the entry is describing, and then the rest of the data, including:
 *     * `TTT` (64-bit integer): (Extended) Trigger Time Tag value, in readout
 *       board ticks (each worth 8 ns) from the last board reset;
 *       currently the value includes the 31-bit counter and in addition the
 *       overflow bit as MSB; the overflow bit is set by the readout board
 *       the first time the counter passes its limit (2^31) and wraps, and never
 *       cleared until the next board reset.
 *       While the tree has room for the 48-bit version (ETTT), the rest of the
 *       decoder does not yet support it.
 *     * `trigger` (64-bit signed integer): global trigger time (from the
 *       trigger decoder), in nanoseconds from The Epoch.
 *     * `triggerSec`, `triggerNS` (32-bit integer each): same time as `trigger`
 *       branch, split into second and nanosecond components.
 *     * `fragTime` (64-bit signed integer), `fragTimeSec` (32-bit signed
 *       integer): the timestamp of the PMT fragment, assigned by the board
 *       reader constructing the fragment.
 * 
 * 
 * 
 * Technical notes
 * ----------------
 * 
 * In order to correctly reconstruct the time stamp, this module needs several
 * pieces of information.
 * These include the size of the pre-trigger buffer, which is set by the readout
 * board configuration, and the delay between the global trigger and the time
 * that trigger is received and acted upon in the readout board, which needs to
 * be measured.
 * The first category of information, from readout board configuration, are read
 * from the input file (`sbn::PMTconfiguration`), while the second category 
 * needs to be specified in the tool FHiCL configuration.
 * 
 * PMT configuration is optional, in the sense that it can be omitted; in that
 * case, some standard values will be used for it.
 * For a board to be served, an entry of that board must be present in the
 * tool configuration (`BoardSetup`). It is an error for a fragment in input not
 * to have an entry for the corresponding board setup.
 * 
 * The module code extract the needed information and matches it into a
 * sort-of-database keyed by fragment ID, so that it can be quickly applied
 * when decoding a fragment. The matching is performed by board name.
 * 
 * 
 * Glossary
 * ---------
 * 
 * * **setup**, **[PMT] configuration**: this is jargon specific to this tool.
 *     Information about a readout board can come from two sources: the "setup"
 *     is information included in the `BoardSetup` configuration list of this
 *     tool; the "PMT configuration" is information included in the DAQ
 *     configuration that is delivered via `PMTconfigTag`.
 * * **TTT**: trigger time tag, from the V1730 event record (31 bits); may be:
 * * **ETTT**: extended trigger time tag, from the V1730 event record (48 bits).
 * * **trigger delay**: time point when a V1730 board processes a (PMT) trigger
 *     signal (and increments the TTT register) with respect to the time of the
 *     time stamp of the (SPEXi) global trigger that acquired the event.
 * 
 * 
 * @todo Merge contiguous waveforms on the same channel
 * @todo Add interface for fragment containers
 * 
 * 
 */
class icarus::DaqDecoderICARUSPMT: public art::EDProducer {
  
  // --- BEGIN -- some debugging tree declarations -----------------------------
  
  /// Enumerate the supported data trees.
  enum class DataTrees: std::size_t {
    Fragments, ///< Information about fragments
    N          ///< Counter.
  };
  using TreeNameList_t
    = std::array<std::string, static_cast<std::size_t>(DataTrees::N)>;
  static TreeNameList_t const TreeNames;
  
  /// Returns a string with all supported tree names.
  static std::string listTreeNames(std::string const& sep = "\n");
  
  // --- END ---- some debugging tree declarations -----------------------------
  
    public:
  
  // --- BEGIN -- public data types --------------------------------------------
  
  using nanoseconds = util::quantities::intervals::nanoseconds; ///< Alias.
  
  /// Data structure for trigger time.
  struct SplitTimestamp_t {
    static_assert(sizeof(int) >= 4U);
    struct Split_t {
        int seconds = std::numeric_limits<int>::min(); ///< The ongoing second.
        /// Nanoseconds from the start of current second.
        unsigned int nanoseconds = std::numeric_limits<unsigned int>::max();
    }; // Split_t
    
    /// Trigger time in nanoseconds from The Epoch.
    long long int time = std::numeric_limits<long long int>::min();
    /// Trigger time in nanoseconds from The Epoch (in components).
    Split_t split;
    
    constexpr SplitTimestamp_t() = default;
    constexpr SplitTimestamp_t(int sec, unsigned int ns);
    constexpr SplitTimestamp_t(long long int triggerTime);
  }; // SplitTimestamp_t
  
  // --- END ---- public data types --------------------------------------------
  
  
  // --- BEGIN -- FHiCL configuration ------------------------------------------
  
  /// Configuration of the V1730 readout board setup.
  struct BoardSetupConfig {
    
    fhicl::Atom<std::string> Name {
      fhicl::Name("Name"),
      fhicl::Comment("board name, as specified in the DAQ configuration")
      };
    
    fhicl::OptionalAtom<unsigned int> FragmentID {
      fhicl::Name("FragmentID"),
      fhicl::Comment("ID of the fragments associated with the board")
      };
    
    fhicl::Atom<nanoseconds> TriggerDelay {
      fhicl::Name("TriggerDelay"),
      fhicl::Comment
        ("from delay from the trigger timestamp to the PMT trigger [ns]"),
      0_ns
      };
    
    fhicl::Atom<nanoseconds> TTTresetDelay {
      fhicl::Name("TTTresetDelay"),
      fhicl::Comment
        ("assume that V1730 counter (Trigger Time Tag) is reset every second"),
      0_ns
      };
    
  }; // struct BoardSetupConfig
  
  
  /// Main tool configuration.
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<art::InputTag> FragmentsLabel {
      Name("FragmentsLabel"),
      Comment("data product with the PMT fragments from DAQ"),
      "daq:CAEN1730" // default
      };
    
    fhicl::Atom<bool> SurviveExceptions {
      Name("SurviveExceptions"),
      Comment
        ("when the decoding tool throws an exception, print a message and move on"),
      true // default
      };
    
    fhicl::Atom<bool> DiagnosticOutput {
      Name("DiagnosticOutput"),
      Comment("enable additional console output"),
      false // default
      };
    
    fhicl::Atom<bool> RequireKnownBoards {
      Name("RequireKnownBoards"),
      Comment
        ("all readout boards in input must be known (setup+PMT configuration)"),
      true
      };
    
    fhicl::Atom<bool> RequireBoardConfig {
      Name("RequireBoardConfig"),
      Comment
        ("all readout boards in setup must have a matching PMT configuration"),
      true
      };
    
    fhicl::OptionalAtom<art::InputTag> PMTconfigTag {
      Name("PMTconfigTag"),
      Comment("input tag for the PMT readout board configuration information")
      };
    
    fhicl::Sequence
      <fhicl::TableAs<daq::details::BoardSetup_t, BoardSetupConfig>>
    BoardSetup {
      Name("BoardSetup"),
      Comment("list of the setup settings for all relevant V1730 boards")
      };
    
    fhicl::Atom<art::InputTag> TriggerTag {
      Name("TriggerTag"),
      Comment("input tag for the global trigger object (raw::ExternalTrigger)")
      };
    
    fhicl::OptionalAtom<bool> TTTresetEverySecond {
      Name("TTTresetEverySecond"),
      Comment
        ("assume that V1730 counter (Trigger Time Tag) is reset every second")
      };
    
    fhicl::Sequence<std::string> DataTrees {
      fhicl::Name("DataTrees"),
      fhicl::Comment
        ("produces the specified ROOT trees (" + listTreeNames(",") + ")"),
      std::vector<std::string>{} // default
      };
    
    fhicl::Atom<std::string> LogCategory {
      Name("LogCategory"),
      Comment("name of the category for message stream"),
      "PMTDecoder" // default
      };
    
  }; // Config
  
  using Parameters = art::EDProducer::Table<Config>;
  
  
  static constexpr detinfo::timescales::electronics_time NoTimestamp
    = std::numeric_limits<detinfo::timescales::electronics_time>::min();
  
  
  // --- END ---- FHiCL configuration ------------------------------------------
  
  
  /// Constructor.
  explicit DaqDecoderICARUSPMT(Parameters const& params);
  
  /// On a new run: cache PMT configuration information.
  void beginRun(art::Run& run) override;
  
  /// Processes the event.
  void produce(art::Event& event) override;
  
  
    private:
  
  /// Information used in decoding from a board.
  struct NeededBoardInfo_t {
    std::string const name;
    nanoseconds bufferLength;
    nanoseconds preTriggerTime;
    nanoseconds PMTtriggerDelay;
    nanoseconds TTTresetDelay;
  };
  
  // --- BEGIN -- Configuration parameters -------------------------------------
  
  art::InputTag const fInputTag; ///< Data product with artDAQ data fragments.
  
  bool const fSurviveExceptions; ///< Whether to "ignore" errors.
  
  /// If true will spew endless messages to output.
  bool const fDiagnosticOutput;
  
  /// Whether info on all input boards is required.
  bool const fRequireKnownBoards;
  
  /// Whether setup info on all boards is required.
  bool const fRequireBoardConfig;
  
  /// Input tag of the PMT configuration.
  std::optional<art::InputTag> const fPMTconfigTag;
  
  art::InputTag const fTriggerTag; ///< Input tag of the global trigger.
  
  bool const fTTTresetEverySecond; ///< Whether V1730 TTT is reset every second.
  
  /// All board setup settings.
  std::vector<daq::details::BoardSetup_t> const fBoardSetup;
  
  std::string const fLogCategory; ///< Message facility category.
  
  // --- END ---- Configuration parameters -------------------------------------

  
  // --- BEGIN -- Services -----------------------------------------------------
  
  /// Interface to LArSoft configuration for detector timing.
  detinfo::DetectorTimings const fDetTimings;
  
  /// Fragment/channel mapping database.
  icarusDB::IICARUSChannelMap const& fChannelMap;

  // --- END ---- Services -----------------------------------------------------


  // --- BEGIN -- Cached values ------------------------------------------------
  
  /// Duration of the optical detector readout sampling tick (i.e. 2 ns; hush!).
  nanoseconds const fOpticalTick;
  
  /// Trigger time as reported by `DetectorClocks` service.
  detinfo::timescales::electronics_time const fNominalTriggerTime;
  
  // --- END ---- Cached values ------------------------------------------------
  
  
  // --- BEGIN -- Per-run data cache -------------------------------------------
  
  /// Find the information on a readout boards by fragment ID.
  std::optional<daq::details::BoardInfoLookup> fBoardInfoLookup;
  
  // --- END -- Per-run data cache ---------------------------------------------
  
  
  // --- BEGIN -- PMT readout configuration ------------------------------------
  
  /// Returns whether PMT configuration information is expected to be available.
  bool hasPMTconfiguration() const { return fPMTconfigTag.has_value(); }
  
  /// Updates the PMT configuration cache. How? Dunno. Placeholder.
  bool UpdatePMTConfiguration(sbn::PMTconfiguration const* PMTconfig);
  
  
  /**
   * @brief Returns a lookup object with board setup and configuration info.
   * @param PMTconfig the PMT configuration, if available
   * @return an object working like lookup table for all fragment information
   * 
   * This method merges the setup information from the tool configuration with
   * the PMT configuration specified in the argument, and returns an object
   * that can look up all the information as a single record, with the
   * fragment ID as key. In addition, a few intermediate quantities ("facts",
   * see `BoardFacts_t`) are computed and stored in this object.
   * 
   * If a fragment ID is missing, it means that no PMT configuration was
   * provided and that the setup information did not include a fragment ID.
   * If some information (configuration or setup) is missing, the "facts"
   * depending on the missing information will have default values.
   */
  daq::details::BoardInfoLookup matchBoardConfigurationAndSetup
    (sbn::PMTconfiguration const* PMTconfig) const;
  
  /// Puts together all the needed information for a board.
  NeededBoardInfo_t fetchNeededBoardInfo(
    daq::details::BoardInfoLookup::BoardInfo_t const* boardInfo,
    unsigned int fragmentID
    ) const;

  /// Extracts the Trigger Time Tag (31+1 bits) value from the fragment
  static unsigned int extractTriggerTimeTag(artdaq::Fragment const& fragment);


  // --- END -- PMT readout configuration --------------------------------------


  /**
   * @brief Returns the count of set bits for each set bit.
   * @tparam NBits the number of bits to test
   * @param value the bit mask to be analyzed
   * @return a pair: `first` an array with the `NBits` values,
   *         `second` the number of bits set
   * 
   * Much better to go with an example.
   * Let's process `setBitIndices<16U>(0x6701)`: `value` is `0110 0111 0000 0001`,
   * and we request the `16` least significant bits (`NBits`).
   * There are six bit that are set, in position 0, 8, 9, 10, 13 and 14.
   * The returned value includes the number of set bits (6) as the `second`
   * element of the pair, an as `first` an array of 16 indices, one for each bit
   * position.
   * Its values are `0` at `[0]`, `1` at `[8]`, `2` at `[9]`, `3` at `[10]`,
   * `4` at `[13]` and `5` at `[14]`. That is, each value represents the number
   * of that bit in the list of set bits. All the other indices, associated to
   * the bits which are not set, are assigned a value that is equal or larger than
   * the number of bits set (i.e. `6` or larger).
   */
  template <std::size_t NBits, typename T>
  static constexpr std::pair<std::array<std::size_t, NBits>, std::size_t>
    setBitIndices(T value) noexcept;

  
  // --- BEGIN -- Trees and their data ---------------------------------------
  
  /// Data structure for basic event information in simple ROOT trees.
  struct TreeData_EventID_t {
    unsigned int run;    ///< Run number.
    unsigned int subrun; ///< Subrun number.
    unsigned int event;  ///< Event number.
    SplitTimestamp_t timestamp; ///< Event timestamp (seconds from the epoch).
  }; // TreeData_EventID_t
  
  /// Structure collecting all data for a fragment ROOT tree.
  struct TreeFragment_t {
    struct Data_t: public TreeData_EventID_t {
      
      int fragmentID = 0; ///< ID of the fragment of this entry.
      
      ///< Trigger time tag from the fragment.
      unsigned long int TriggerTimeTag = 0;
      
      SplitTimestamp_t trigger; ///< Global trigger time.
      
      SplitTimestamp_t fragTime; ///< PMT fragment time stamp.
      
    }; // Data_t
    
    Data_t data;
    TTree* tree = nullptr;
  }; // TreeFragment_t
  
  
  std::unique_ptr<TreeData_EventID_t> fEventInfo; ///< Event ID for trees.
  
  ///< Tree with fragment information.
  std::unique_ptr<TreeFragment_t> fTreeFragment;
  
  // --- END ---- Trees and their data -----------------------------------------
  
  
  // --- BEGIN -- Timestamps ---------------------------------------------------
  
  /// Retrieves the global trigger time stamp from the event.
  SplitTimestamp_t fetchTriggerTimestamp(art::Event const& event) const;
  
  /// Returns the timestamp for the waveforms in the specified fragment.
  detinfo::timescales::electronics_time fragmentWaveformTimestamp(
    artdaq::Fragment const& artdaqFragment,
    NeededBoardInfo_t const& boardInfo,
    SplitTimestamp_t triggerTime
    ) const;
  
  /**
   * @brief Returns the timestamp for the waveforms in the specified fragment.
   * 
   * This method assumes that all PMT readout board triggers happened at the
   * same time as the global trigger.
   * Unsurprisingly, this does not work well with multiple PMT windows.
   */
  detinfo::timescales::electronics_time fragmentWaveformTimestampOnTrigger(
    artdaq::Fragment const& artdaqFragment,
    NeededBoardInfo_t const& boardInfo,
    SplitTimestamp_t triggerTime
    ) const;
  
  /**
   * @brief Returns the timestamp for the waveforms in the specified fragment.
   * 
   * This method assumes that the Trigger Time Tag counters are synchronised
   * with the global trigger and their value is reset on each new second of 
   * the global trigger timescale (International Atomic Time).
   * 
   * This assumptions enables timestamping of waveforms from the same readout
   * boards at different times ("multi-window PMT readout").
   * 
   * See `TTTresetEverySecond` configuration option.
   */
  detinfo::timescales::electronics_time fragmentWaveformTimestampFromTTT(
    artdaq::Fragment const& artdaqFragment,
    NeededBoardInfo_t const& boardInfo,
    SplitTimestamp_t triggerTime
    ) const;
  
  /// Returns the fragment ID to be used with databases.
  static constexpr std::size_t effectivePMTboardFragmentID
    (artdaq::Fragment::fragment_id_t id)
    { return id & 0x0fff; }
  
  // --- END ---- Timestamps ---------------------------------------------------
  
  
  /// Collection of useful information from fragment data.
  struct FragmentInfo_t {
    artdaq::Fragment::fragment_id_t fragmentID
      = std::numeric_limits<artdaq::Fragment::fragment_id_t>::max();
    artdaq::Fragment::timestamp_t fragmentTimestamp;
    std::uint32_t TTT = 0U;
    std::uint16_t enabledChannels = 0U;
    std::size_t nSamplesPerChannel = 0U;
    std::uint16_t const* data = nullptr;
  }; // FragmentInfo_t

  /// Extracts waveforms from the specified fragments from a board.
  std::vector<raw::OpDetWaveform> processBoardFragments(
    std::vector<artdaq::Fragment const*> const& artdaqFragment,
    SplitTimestamp_t triggerTime
    );
  
  /**
   * @brief Create waveforms and fills trees for the specified artDAQ fragment.
   * @param artdaqFragment the fragment to process
   * @param boardInfo board information needed, from configuration/setup
   * @param triggerTime absolute time of the trigger
   * @return collection of PMT waveforms from the fragment
   * 
   * This method fills the information for the PMT fragment tree
   * (`fillPMTfragmentTree()`) and creates PMT waveforms from the fragment data
   * (`createFragmentWaveforms()`).
   */
  std::vector<raw::OpDetWaveform> processFragment(
    artdaq::Fragment const& artdaqFragment,
    NeededBoardInfo_t const& boardInfo,
    SplitTimestamp_t triggerTime
    );

  /**
   * @brief Creates `raw::OpDetWaveform` objects from the fragment data.
   * @param fragInfo information extracted from the fragment
   * @param timeStamp timestamp of the waveforms in the fragment
   * @return collection of newly created `raw::OpDetWaveform`
   * 
   * All fragment information needed is enclosed in `fragInfo`
   * (`extractFragmentInfo()`). The timestamp can be obtained with a call to
   * `fragmentWaveformTimestamp()`.
   */
  std::vector<raw::OpDetWaveform> createFragmentWaveforms(
    FragmentInfo_t const& fragInfo,
    detinfo::timescales::electronics_time const timeStamp
    ) const;
  
  /// Extracts useful information from fragment data.
  FragmentInfo_t extractFragmentInfo
    (artdaq::Fragment const& artdaqFragment) const;
  
  /// Returns the board information for this fragment.
  NeededBoardInfo_t neededBoardInfo
    (artdaq::Fragment::fragment_id_t fragment_id) const;
  

  /// Sorts in place the specified waveforms in channel order, then in time.
  void sortWaveforms(std::vector<raw::OpDetWaveform>& waveforms) const;
  
  
  // --- BEGIN -- Tree-related methods -----------------------------------------
  
  /// Declares the use of event information.
  void usesEventInfo();
  
  /// Initializes all requested data trees.
  void initTrees(std::vector<std::string> const& treeNames);
  
  /// Initializes the event ID part of a tree.
  void initEventIDtree(TTree& tree, TreeData_EventID_t& data);
  
  /// Initializes the fragment data tree (`fTreeFragment`).
  void initFragmentsTree();

  /// Fills the base information of a tree data entry from an _art_ event.
  void fillTreeEventID
    (art::Event const& event, TreeData_EventID_t& treeData) const;

  /// Assigns the cached event information to the specified tree data.
  void assignEventInfo(TreeData_EventID_t& treeData) const;
  
  /// Fills the PMT fragment tree with the specified information
  /// (additional information needs to have been set already).
  void fillPMTfragmentTree
    (FragmentInfo_t const& fragInfo, SplitTimestamp_t triggerTime);
  

  /// Static initialization.
  static TreeNameList_t initTreeNames();
  
  // --- END ---- Tree-related methods -----------------------------------------
  
  
}; // icarus::DaqDecoderICARUSPMT


//------------------------------------------------------------------------------
// --- implementation
//------------------------------------------------------------------------------
namespace {
  
  /// Moves the contend of `src` into the end of `dest`.
  template <typename T>
  std::vector<T>& appendTo(std::vector<T>& dest, std::vector<T>&& src) {
    if (dest.empty()) dest = std::move(src);
    else {
      dest.reserve(dest.size() + src.size());
      std::move(src.begin(), src.end(), std::back_inserter(dest));
    }
    src.clear();
    return dest;
  } // appendTo()
  
} // local namespace


//------------------------------------------------------------------------------
namespace icarus {

  /// Special function `fhicl::TableAs` uses to convert BoardSetupConfig.
  daq::details::BoardSetup_t convert
    (DaqDecoderICARUSPMT::BoardSetupConfig const& config)
  {
    
    using ::util::fhicl::getOptionalValue;
    return {
        config.Name()                                          // name
      , getOptionalValue(config.FragmentID) 
          .value_or(daq::details::BoardSetup_t::NoFragmentID)  // fragmentID
      , config.TriggerDelay()                                  // triggerDelay
      , config.TTTresetDelay()                                 // TTTresetDelay
      };
  } // convert(BoardSetupConfig)

} // namespace icarus


//------------------------------------------------------------------------------
// --- icarus::DaqDecoderICARUSPMT::SplitTimestamp_t
//------------------------------------------------------------------------------
constexpr icarus::DaqDecoderICARUSPMT::SplitTimestamp_t::SplitTimestamp_t
  (int sec, unsigned int ns)
  : time { static_cast<long long int>(sec) * 1'000'000'000LL + ns }
  , split { sec, ns }
{}


//------------------------------------------------------------------------------
constexpr icarus::DaqDecoderICARUSPMT::SplitTimestamp_t::SplitTimestamp_t
  (long long int triggerTime)
  : time { triggerTime }
  , split {
      static_cast<int>(time / 1'000'000'000), // seconds
      static_cast<unsigned int>(time % 1'000'000'000) // nanoseconds
    }
{}


//------------------------------------------------------------------------------
namespace icarus {
  
  std::ostream& operator<<
    (std::ostream& out, DaqDecoderICARUSPMT::SplitTimestamp_t const& time)
  {
    out << time.split.seconds << '.'
      << std::setfill('0') << std::setw(9) << time.split.nanoseconds;
    return out;
  } // operator<< (std::ostream&, PMTDecoder::SplitTimestamp_t)
  
} // namespace icarus


//------------------------------------------------------------------------------
// --- icarus::DaqDecoderICARUSPMT
//------------------------------------------------------------------------------
icarus::DaqDecoderICARUSPMT::TreeNameList_t const
icarus::DaqDecoderICARUSPMT::TreeNames
  = icarus::DaqDecoderICARUSPMT::initTreeNames();

auto icarus::DaqDecoderICARUSPMT::initTreeNames() -> TreeNameList_t {
  TreeNameList_t names;
  names[static_cast<std::size_t>(DataTrees::Fragments)] = "PMTfragments";
  return names;
} // icarus::DaqDecoderICARUSPMT::initTreeNames()


//------------------------------------------------------------------------------
std::string icarus::DaqDecoderICARUSPMT::listTreeNames
  (std::string const& sep /* = " " */)
{
  std::string l;
  for (std::string const& name: TreeNames) {
    if (!l.empty()) l += sep;
    l += '\'';
    l += name;
    l += '\'';
  } // for
  return l;
} // icarus::DaqDecoderICARUSPMT::listTreeNames()


//------------------------------------------------------------------------------
icarus::DaqDecoderICARUSPMT::DaqDecoderICARUSPMT(Parameters const& params)
  : art::EDProducer(params)
  , fInputTag{ params().FragmentsLabel() }
  , fSurviveExceptions{ params().SurviveExceptions() }
  , fDiagnosticOutput{ params().DiagnosticOutput() }
  , fRequireKnownBoards{ params().RequireKnownBoards() }
  , fRequireBoardConfig{ params().RequireBoardConfig() }
  , fPMTconfigTag{ ::util::fhicl::getOptionalValue(params().PMTconfigTag) }
  , fTriggerTag{ params().TriggerTag() }
  , fTTTresetEverySecond{
    ::util::fhicl::getOptionalValue(params().TTTresetEverySecond)
      .value_or(!fTriggerTag.empty())
    }
  , fBoardSetup{ params().BoardSetup() }
  , fLogCategory{ params().LogCategory() }
  , fDetTimings
    { art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob() }
  , fChannelMap{ *(art::ServiceHandle<icarusDB::IICARUSChannelMap const>{}) }
  , fOpticalTick{ fDetTimings.OpticalClockPeriod() }
  , fNominalTriggerTime{ fDetTimings.TriggerTime() }
{
  //
  // consumed data products declaration
  //
  if (fPMTconfigTag) consumes<sbn::PMTconfiguration>(*fPMTconfigTag);
  consumes<std::vector<raw::ExternalTrigger>>(fTriggerTag);
  
  //
  // produced data products declaration
  //
  produces<std::vector<raw::OpDetWaveform>>();
  
  //
  // additional initialization
  //
  initTrees(params().DataTrees());
  
  
  //
  // configuration dump
  //
  mf::LogInfo log(fLogCategory);
  log << "Configuration:"
    << "\n * boards with setup: " << fBoardSetup.size();
  if (fPMTconfigTag)
    log << "\n * PMT configuration from '" << fPMTconfigTag->encode() << "'";
  else 
    log << "\n * PMT configuration not used (and some corrections will be skipped)";
  if (fRequireKnownBoards) {
    log << "\n * all readout boards in input must be known (from `"
      << params().BoardSetup.name() << "` or `"
      << params().PMTconfigTag.name() << "`)"
      ;
  }
  else {
    log << "\n * readout boards with no information (from neither `"
      << params().BoardSetup.name() << "` or `"
      << params().PMTconfigTag.name()
      << "`) are processed at the best we can (skipping corrections)"
      ;
  }
  if (fRequireBoardConfig) {
    log << "\n * all readout boards in `"
      << params().BoardSetup.name()
      << "` must appear in the PMT configuration from `"
      << params().PMTconfigTag.name() << "`"
      ;
  }
  else {
    log << "\n * all readout boards in `"
      << params().BoardSetup.name()
      << "` may lack a matching PMT configuration from `"
      << params().PMTconfigTag.name() << "`"
      ;
  }
  
} // icarus::DaqDecoderICARUSPMT::DaqDecoderICARUSPMT()


//------------------------------------------------------------------------------
void icarus::DaqDecoderICARUSPMT::beginRun(art::Run& run) {
  
  sbn::PMTconfiguration const* PMTconfig = fPMTconfigTag
    ? run.getPointerByLabel<sbn::PMTconfiguration>(*fPMTconfigTag): nullptr;
  
  UpdatePMTConfiguration(PMTconfig);
  
} // icarus::DaqDecoderICARUSPMT::beginRun()


//------------------------------------------------------------------------------
void icarus::DaqDecoderICARUSPMT::produce(art::Event& event) {
  
  // ---------------------------------------------------------------------------
  // preparation
  //
  
  //
  // global trigger
  //
  SplitTimestamp_t const triggerTimestamp = fetchTriggerTimestamp(event);
  mf::LogDebug(fLogCategory)
    << "Trigger time ('" << fTriggerTag.encode() << "'): "
    << triggerTimestamp << " s"
    ;
  
  //
  // event ID
  //
  
  // if needed, fill the record with the basic information of the event
  if (fEventInfo) fillTreeEventID(event, *fEventInfo);
  
  //
  // output data product initialization
  //
  std::vector<raw::OpDetWaveform> opDetWaveforms;
  
  
  // ---------------------------------------------------------------------------
  // pre-processing
  //
  
  
  // ---------------------------------------------------------------------------
  // processing
  //
  try {
    
    auto const& fragments = event.getByLabel<artdaq::Fragments>(fInputTag);
    
    for (artdaq::Fragment const& fragment: fragments) {
      std::vector<artdaq::Fragment const*> boardFragments { &fragment };
      appendTo(
        opDetWaveforms, processBoardFragments(boardFragments, triggerTimestamp)
        );
    } // for
    
  }
  catch (cet::exception const& e) {
    if (!fSurviveExceptions) throw;
    mf::LogError("DaqDecoderICARUSPMT")
      << "Error while attempting to decode PMT data:\n" << e.what() << '\n';
    opDetWaveforms.clear();
  }
  catch (...) {
    if (!fSurviveExceptions) throw;
    mf::LogError("DaqDecoderICARUSPMT")
      << "Error while attempting to decode PMT data.\n";
    opDetWaveforms.clear();
  }
  
  //
  // post-processing
  //
  sortWaveforms(opDetWaveforms);
  
  // ---------------------------------------------------------------------------
  // output
  //
  event.put(
    std::make_unique<std::vector<raw::OpDetWaveform>>(std::move(opDetWaveforms))
    );
  
} // icarus::DaqDecoderICARUSPMT::produce()


//------------------------------------------------------------------------------
bool icarus::DaqDecoderICARUSPMT::UpdatePMTConfiguration
  (sbn::PMTconfiguration const* PMTconfig)
{
  fBoardInfoLookup.emplace(matchBoardConfigurationAndSetup(PMTconfig));
  
  mf::LogDebug(fLogCategory)
    << "Board information as cached:\n" << *fBoardInfoLookup;
  
  return true;
} // icarus::DaqDecoderICARUSPMT::UpdatePMTConfiguration()


auto icarus::DaqDecoderICARUSPMT::matchBoardConfigurationAndSetup
  (sbn::PMTconfiguration const* PMTconfig) const
  -> daq::details::BoardInfoLookup
{
  /*
   * We need to support the case where no PMT configuration is known
   * (that is the standard situation in the online monitor).
   * The "strategy" is that in such cases we give up the correct time stamp
   * decoding; if the setup information contains a fragment ID, it may be
   * possible to do a little better, that is to use the setup information
   * (this is not possible without knowing the fragment ID that each bit of
   * setup information pertains).
   * 
   * So the cases for a board are:
   * * setup information is not present: encountering such a board will cause
   *   an exception to be thrown (implemented elsewhere)
   * * PMT configuration and setup present: full configuration
   *     * exception thrown if setup fragment ID is present and inconsistent
   * * PMT configuration not present: a general warning is printed;
   *     * boards with setup fragment ID information: add setup information
   *       to the "database" for the board: it will be used for partial
   *       timestamp reconstruction
   *     * boards without setup fragment ID information: board will not be
   *       added into the database; no specific correction will be performed;
   *       a warning is printed for each board
   * 
   */
  
  // dictionary of board configurations (if any)
  std::vector<std::pair<std::string, sbn::V1730Configuration const*>>
    configByName;
  if (PMTconfig) {
    if (!PMTconfig->boards.empty())
      configByName.reserve(PMTconfig->boards.size());
    for (sbn::V1730Configuration const& boardConfig: PMTconfig->boards)
      configByName.emplace_back(boardConfig.boardName, &boardConfig);
    std::sort(configByName.begin(), configByName.end()); // sorted by board name
  } // if we have configuration
  
  
  auto findPMTconfig = [this, &configByName]
      (std::string const& name) -> sbn::V1730Configuration const*
    {
      if (!hasPMTconfiguration()) return nullptr;
      auto const* ppBoardConfig
        = daq::details::binarySearch(configByName, name);
      if (!ppBoardConfig) {
        if (!fRequireBoardConfig) return nullptr;
        throw cet::exception("PMTDecoder")
          << "No DAQ configuration found for PMT readout board '"
          << name << "'\n"
          << "If this is expected, you may skip this check by setting "
          << "PMTDecoder tool configuration `RequireBoardConfig` to `false`.\n";
      }
      return ppBoardConfig->second;
    }; // findPMTconfig()
  
  // the filling is driven by boards configured in the tool
  // (which is how a setup entry is mandatory)
  daq::details::BoardInfoLookup::Database_t boardInfoByFragment;
  
  for (daq::details::BoardSetup_t const& boardSetup: fBoardSetup) {
    
    std::string const& boardName = boardSetup.name;
    
    sbn::V1730Configuration const* pBoardConfig = findPMTconfig(boardName);
    
    if (pBoardConfig) {
      // fragment ID from configuration and setup must match if both present
      if (boardSetup.hasFragmentID()
        && (boardSetup.fragmentID != pBoardConfig->fragmentID)
      ) {
        throw cet::exception("PMTDecoder")
          << "Board '" << boardName << "' has fragment ID "
          << std::hex << pBoardConfig->fragmentID << std::dec
          << " but it is set up as "
          << std::hex << boardSetup.fragmentID << std::dec
          << "!\n";
      } // if fragment ID in setup
    }
    else {
      if (boardSetup.hasFragmentID()) {
        mf::LogPrint(fLogCategory)
          << "Board '" << boardName
          << "' has no configuration information;"
            " some time stamp corrections will be skipped.";
        // to avoid this, make a PMT configuration available from input file
      }
      else {
        mf::LogPrint(fLogCategory)
          << "Board '" << boardName
          << "' can't be associated to a fragment ID;"
            " its time stamp corrections will be skipped.";
        // to avoid this, add a `BoardSetup.FragmentID` entry for it in the
        // configuration of this tool, or make a PMT configuration available
        continue; // no entry for this board at all
      }
    }
    
    unsigned int const fragmentID
      = pBoardConfig? pBoardConfig->fragmentID: boardSetup.fragmentID;
    assert(fragmentID != daq::details::BoardSetup_t::NoFragmentID);
    
    nanoseconds const preTriggerTime
      = pBoardConfig
      ? (pBoardConfig->bufferLength * (1.0f - pBoardConfig->postTriggerFrac))
        * fOpticalTick
      : nanoseconds{ 0.0 }
      ;
    
    daq::details::BoardFacts_t boardFacts {
      preTriggerTime  // ditto
      };
    
    boardInfoByFragment.push_back({
      fragmentID,               // fragmentID
      &boardSetup,              // setup
      pBoardConfig,             // config
      std::move(boardFacts)     // facts
      });
  } // for
  
  return daq::details::BoardInfoLookup{ std::move(boardInfoByFragment) };
  
} // icarus::DaqDecoderICARUSPMT::matchBoardConfigurationAndSetup()


//------------------------------------------------------------------------------
auto icarus::DaqDecoderICARUSPMT::fetchNeededBoardInfo(
  daq::details::BoardInfoLookup::BoardInfo_t const* boardInfo,
  unsigned int fragmentID
) const -> NeededBoardInfo_t {
  
  using util::quantities::intervals::microseconds;
  
  return NeededBoardInfo_t{
    // name
      ((boardInfo && boardInfo->config)
        ? boardInfo->config->boardName: ("<ID=" + std::to_string(fragmentID)))
    // bufferLength
    , ((boardInfo && boardInfo->config)
        ? boardInfo->config->bufferLength * fOpticalTick: nanoseconds{ 0.0 }
      )
    // preTriggerTime
    , (boardInfo? boardInfo->facts.preTriggerTime: nanoseconds{ 0.0 })
    // PMTtriggerDelay
    , ((boardInfo && boardInfo->setup)
        ? boardInfo->setup->triggerDelay: nanoseconds{ 0.0 })
    // TTTresetDelay
    , ((boardInfo && boardInfo->setup)
        ? boardInfo->setup->TTTresetDelay: nanoseconds{ 0.0 })
    };
      
} // icarus::DaqDecoderICARUSPMT::fetchNeededBoardInfo()


//------------------------------------------------------------------------------
auto icarus::DaqDecoderICARUSPMT::fetchTriggerTimestamp
  (art::Event const& event) const -> SplitTimestamp_t
{
  auto const& triggers
    = event.getByLabel<std::vector<raw::ExternalTrigger>>(fTriggerTag);
  if (triggers.size() != 1U) {
    // if this is hit, the decoder needs some development to correctly deal
    // with input with no trigger, or more than one
    throw cet::exception("PMTDecoder")
      << "Found " << triggers.size() << " triggers from '"
      << fTriggerTag.encode() << "', can deal only with 1.\n";
  }
  return { triggers.front().GetTrigTime() };
} // icarus::DaqDecoderICARUSPMT::fetchTriggerTimestamp()


//------------------------------------------------------------------------------
auto icarus::DaqDecoderICARUSPMT::processBoardFragments(
  std::vector<artdaq::Fragment const*> const& artdaqFragments,
  SplitTimestamp_t triggerTime
) -> std::vector<raw::OpDetWaveform> {
  
  if (artdaqFragments.empty()) return {};
  
  assert(artdaqFragments.front());
  
  NeededBoardInfo_t const boardInfo
    = neededBoardInfo(artdaqFragments.front()->fragmentID());
  
  std::vector<raw::OpDetWaveform> waveforms;
  for (artdaq::Fragment const* fragment: artdaqFragments)
    appendTo(waveforms, processFragment(*fragment, boardInfo, triggerTime));
  
  return waveforms;
  
} // icarus::DaqDecoderICARUSPMT::processBoardFragments()


//------------------------------------------------------------------------------
auto icarus::DaqDecoderICARUSPMT::processFragment(
  artdaq::Fragment const& artdaqFragment,
  NeededBoardInfo_t const& boardInfo,
  SplitTimestamp_t triggerTime
) -> std::vector<raw::OpDetWaveform> {
  
  FragmentInfo_t const fragInfo = extractFragmentInfo(artdaqFragment);
  
  if (fTreeFragment) fillPMTfragmentTree(fragInfo, triggerTime);
  
  auto const timeStamp
    = fragmentWaveformTimestamp(artdaqFragment, boardInfo, triggerTime);
  return (timeStamp != NoTimestamp)
    ? createFragmentWaveforms(fragInfo, timeStamp)
    : std::vector<raw::OpDetWaveform>{}
    ;
  
} // icarus::DaqDecoderICARUSPMT::processFragment()


//------------------------------------------------------------------------------
auto icarus::DaqDecoderICARUSPMT::createFragmentWaveforms(
  FragmentInfo_t const& fragInfo,
  detinfo::timescales::electronics_time const timeStamp
) const -> std::vector<raw::OpDetWaveform> {

  assert(timeStamp != NoTimestamp);
  
  auto const [ chDataMap, nEnabledChannels ]
    = setBitIndices<16U>(fragInfo.enabledChannels);

  std::vector<raw::OpDetWaveform> opDetWaveforms; // output collection
    
  std::optional<mf::LogVerbatim> diagOut;
  if (fDiagnosticOutput) diagOut.emplace(fLogCategory);
  
  icarusDB::DigitizerChannelChannelIDPairVec const& digitizerChannelVec
    = fChannelMap.getChannelIDPairVec
      (effectivePMTboardFragmentID(fragInfo.fragmentID))
    ;
  
  if (diagOut)
    (*diagOut) << "      " << digitizerChannelVec.size() << " channels:";
  
  // allocate the vector outside the loop since we'll reuse it over and over
  std::vector<std::uint16_t> wvfm(fragInfo.nSamplesPerChannel);

  // track what we do and what we want to
  std::uint16_t attemptedChannels = 0;
  
  // loop over the channels that we know might be in the fragment
  for(auto const [ digitizerChannel, channelID ]: digitizerChannelVec) {
    
    if (diagOut)
      (*diagOut) << " " << digitizerChannel << " [=> " << channelID << "];";
    
    attemptedChannels |= (1 << digitizerChannel);
    
    // find where this channel is in the data fragment
    std::size_t const channelPosInData = chDataMap[digitizerChannel];
    if (channelPosInData >= nEnabledChannels) {
      mf::LogTrace(fLogCategory)
        << "Digitizer channel " << digitizerChannel
        << " [=> " << channelID << "] skipped because not enabled.";
      continue;
    }
    
    std::size_t const ch_offset
      = channelPosInData * fragInfo.nSamplesPerChannel;
    
    std::copy_n(fragInfo.data + ch_offset, wvfm.size(), wvfm.begin());
    
    mf::LogTrace(fLogCategory)
      << "PMT channel " << channelID << " has " << wvfm.size()
      << " samples (read from entry #" << channelPosInData
      << " in fragment data) starting at electronics time " << timeStamp;
    opDetWaveforms.emplace_back(timeStamp.value(), channelID, wvfm);
    
  } // for channels
  
  if (diagOut) diagOut.reset(); // destroys and therefore prints out
  if (attemptedChannels != fragInfo.enabledChannels) {
    // this is mostly a warning; regularly, for example,
    // we effectively have 15 channels per board; but all 16 are enabled,
    // so one channel is not decoded at all
    mf::LogTrace log(fLogCategory);
    log << "Not all data read:";
    for (int const bit: util::counter(16U)) {
      std::uint16_t const mask = (1 << bit);
      bool const attempted = bool(attemptedChannels & mask);
      bool const enabled = bool(fragInfo.enabledChannels & mask);
      if (attempted == enabled) continue;
      if (!enabled) // and attempted
        log << "\n  requested channel " << bit << " was not enabled";
      if (!attempted) // and enabled
        log << "\n  data for enabled channel " << bit << " was ignored";
    } // for bits
  } // if request and availability did not match

  if (fDiagnosticOutput) {
    mf::LogVerbatim(fLogCategory)
      << "      - number of waveforms decoded: " << opDetWaveforms.size();
  }
  
  return opDetWaveforms;
  
} // icarus::DaqDecoderICARUSPMT::createFragmentWaveforms()


//------------------------------------------------------------------------------
void icarus::DaqDecoderICARUSPMT::fillPMTfragmentTree
  (FragmentInfo_t const& fragInfo, SplitTimestamp_t triggerTime)
{
  
  if (!fTreeFragment) return;
  
  fTreeFragment->data.fragmentID = fragInfo.fragmentID;
  fTreeFragment->data.TriggerTimeTag = fragInfo.TTT;
  fTreeFragment->data.trigger = triggerTime;
  fTreeFragment->data.fragTime
    = { static_cast<long long int>(fragInfo.fragmentTimestamp) };
  assignEventInfo(fTreeFragment->data);
  fTreeFragment->tree->Fill();
  
} // icarus::DaqDecoderICARUSPMT::fillPMTfragmentTree()


//------------------------------------------------------------------------------
auto icarus::DaqDecoderICARUSPMT::extractFragmentInfo
  (artdaq::Fragment const& artdaqFragment) const -> FragmentInfo_t
{
  artdaq::Fragment::fragment_id_t const fragment_id
    = artdaqFragment.fragmentID();

  sbndaq::CAENV1730Fragment const fragment { artdaqFragment };
  sbndaq::CAENV1730FragmentMetadata const& metafrag = *(fragment.Metadata());
  sbndaq::CAENV1730EventHeader const& header = fragment.Event()->Header;

  // chDataMap tells where in the buffer each digitizer channel is;
  // if nowhere, then the answer is a number no smaller than nEnabledChannels
  std::uint16_t const enabledChannels = header.ChannelMask();

  artdaq::Fragment::timestamp_t const fragmentTimestamp
    = artdaqFragment.timestamp();
    
  unsigned int const TTT =  header.triggerTimeTag;
  
  std::size_t const nChannelsPerBoard = metafrag.nChannels;
    
  std::uint32_t const ev_size_quad_bytes = header.eventSize;
  constexpr std::uint32_t evt_header_size_quad_bytes
    = sizeof(sbndaq::CAENV1730EventHeader)/sizeof(std::uint32_t);
  std::uint32_t const data_size_double_bytes
    = 2*(ev_size_quad_bytes - evt_header_size_quad_bytes);
  std::size_t const nSamplesPerChannel
    = data_size_double_bytes/nChannelsPerBoard;
  
  if (fDiagnosticOutput) {
    
    mf::LogVerbatim(fLogCategory)
      << "----> PMT Fragment ID: " << std::hex << fragment_id << std::dec
        << ", nChannelsPerBoard: " << nChannelsPerBoard
        << ", nSamplesPerChannel: " << nSamplesPerChannel
        << ", enabled: " << icarus::ns::util::bin(enabledChannels)
      << "\n      size: " << ev_size_quad_bytes
        << ", data size: " << data_size_double_bytes
        << ", samples/channel: " << nSamplesPerChannel
        << ", trigger time tag: " << TTT
        << ", time stamp: " << (fragmentTimestamp / 1'000'000'000UL)
          << "." << (fragmentTimestamp % 1'000'000'000UL) << " s"
      ;
  } // if diagnostics

  std::uint16_t const* data_begin = reinterpret_cast<std::uint16_t const*>
    (artdaqFragment.dataBeginBytes() + sizeof(sbndaq::CAENV1730EventHeader));

  return { // C++20: write the member names explicitly
    fragment_id,
    fragmentTimestamp,
    TTT,
    enabledChannels,
    nSamplesPerChannel,
    data_begin
    };
  
} // icarus::DaqDecoderICARUSPMT::extractFragmentInfo()


//------------------------------------------------------------------------------
auto icarus::DaqDecoderICARUSPMT::fragmentWaveformTimestamp(
  artdaq::Fragment const& artdaqFragment,
  NeededBoardInfo_t const& boardInfo,
  SplitTimestamp_t triggerTime
) const -> detinfo::timescales::electronics_time {
  
  // check availability of mapping for this board, otherwise can't do anything
  std::size_t const fragment_id = artdaqFragment.fragmentID();
  std::size_t const eff_fragment_id = effectivePMTboardFragmentID(fragment_id);
  
  if (!fChannelMap.hasPMTDigitizerID(eff_fragment_id)) {
    mf::LogError(fLogCategory)
      << "*** PMT could not find channel information for fragment: "
      << artdaqFragment.fragmentID()
      ;
    return NoTimestamp;
  }
  
  if (fTTTresetEverySecond) {
    return
      fragmentWaveformTimestampFromTTT(artdaqFragment, boardInfo, triggerTime);
  }
  else {
    return fragmentWaveformTimestampOnTrigger
      (artdaqFragment, boardInfo, triggerTime);
  }
  
} // icarus::DaqDecoderICARUSPMT::fragmentWaveformTimestamp()


//------------------------------------------------------------------------------
auto icarus::DaqDecoderICARUSPMT::fragmentWaveformTimestampOnTrigger(
  artdaq::Fragment const& artdaqFragment,
  NeededBoardInfo_t const& boardInfo,
  SplitTimestamp_t /* triggerTime */
) const -> detinfo::timescales::electronics_time {
  
  nanoseconds const preTriggerTime = boardInfo.preTriggerTime;
  nanoseconds const PMTtriggerDelay = boardInfo.PMTtriggerDelay;
  
  auto const timestamp = fNominalTriggerTime - PMTtriggerDelay - preTriggerTime;
  mf::LogTrace(fLogCategory) << "V1730 board '" << boardInfo.name
    << "' has data starting at electronics time " << timestamp
    << " = " << fNominalTriggerTime
    << " - " << PMTtriggerDelay << " - " << preTriggerTime
    ;
  return timestamp;
  
} // icarus::DaqDecoderICARUSPMT::fragmentWaveformTimestampOnTrigger()


//------------------------------------------------------------------------------
auto icarus::DaqDecoderICARUSPMT::fragmentWaveformTimestampFromTTT(
  artdaq::Fragment const& artdaqFragment,
  NeededBoardInfo_t const& boardInfo,
  SplitTimestamp_t triggerTime
) const -> detinfo::timescales::electronics_time {
  
  /*
   * 1. goal is a timestamp in electronics time
   * 2. we have the global trigger time in electronics time
   *    (from raw::Trigger data product or from DetectorClocks service)
   * 3. board TTT and global trigger time are on the same timescale:
   *    their difference directly describes the time of the board trigger
   *    relative to the global trigger
   * 4. TTT tags the last (or after the last?) sample of the collected waveform;
   *    the time of the first sample precedes that tag by the full buffer length
   * 5. the PMT trigger itself is subject to a sequence of delays compared to
   *    the (local or global) trigger from SPEXi; here we quantify these delays
   *    from calibration offsets collectively passed via job configuration.
   */
  
  using namespace util::quantities::time_literals;
  
  //
  // 2. global trigger time
  //
  detinfo::timescales::electronics_time waveformTime = fNominalTriggerTime;
  
  //
  // 3. PMT readout board trigger relative to global trigger time
  //
  unsigned int const triggerTimeNS = triggerTime.split.nanoseconds;
  
  // converted into nanoseconds (each TTT tick is 8 nanoseconds):
  unsigned int const TTT = extractTriggerTimeTag(artdaqFragment) * 8;
  
  /*
   * The trigger time tag (TTT) on the PMT readout board is incremented every 8
   * nanoseconds, and the board is sent a reset signal every second, matching
   * the time of the change of second of the global trigger time scale
   * (International Atomic Time from White Rabbit). If the global trigger and
   * the trigger of the PMT readout happen at the same instant (which would be
   * ideally true for one fragment for each board and each event), the TTT will
   * represent exactly the number of nanoseconds of the global trigger passed
   * since the last crossing of a second boundary.
   * (in practice, this is biassed by the signal creation, propagation and
   * interpretation delays, and smeared by the clocks of the SPEXi board sending
   * the signal and the CAEN V1730 readout board, which have a period of 25 ns
   * and 8 ns, respectively).
   *
   * Multiple fragments can be collected from one PMT readout board for the same
   * event by sending the board multiple "local" triggers; these triggers are
   * all supposed to be in the neighbourhood of the global trigger time;
   * in fact, they should not be more than 1 ms away from it, since PMT readout
   * enable window is +/- 1 ms around the expected beam arrival time.
   *
   * In the most common case when global trigger and fragment tag are not
   * separated by a boundary of the second
   * (e.g. global trigger at 1620'284'028 seconds + 799'800'000 nanoseconds
   * and fragment tag 0.5 milliseconds earlier, at 1620'284'028 seconds
   * + 799'300'000 nanoseconds, i.e. with a trigger tag around 799'300'000);
   * in this example, the fragment is 0.5 ms earlier:
   * -500'000 nanoseconds = 799'300'000 (TTT) - 799'800'000 (glob. trigger)
   */
  int fragmentRelTime = static_cast<int>(TTT) - triggerTimeNS;
  if (TTT > triggerTimeNS + 500'000'000U) { // 0.5 seconds
    /*
     * case when global trigger arrives just after the boundary of the second
     * (e.g. global trigger at 1620'284'029 seconds + 300'000 nanoseconds)
     * and this fragment was tagged before that crossing (e.g. 0.5 milliseconds
     * earlier, at 1620'284'028 seconds + 999'800'000 nanoseconds, i.e. with a
     * trigger tag around 999'800'000);
     * in this example, the fragment is 0.5 ms earlier: -500'000 nanoseconds =
     * 999'800'000 (TTT) - 1'000'000'000 (second step) - 300'000 (glob. trigger)
     * and the plain difference,
     * 999'800'000 (TTT) - 300'000 (glob. trigger) = 999'500'000,
     * must be corrected by removing a whole second:
     */
    fragmentRelTime -= 1'000'000'000;
  }
  else if (TTT + 500'000'000U < triggerTimeNS) { // 0.5 seconds
    /*
     * case when global trigger arrives just before the boundary of the second
     * (e.g. global trigger at 1620'284'028 seconds + 999'800'000 nanoseconds)
     * and this fragment was tagged after that crossing (e.g. 0.5 milliseconds
     * later, at 1620'284'029 seconds + 300'000 nanoseconds, i.e. with a
     * trigger tag around 300'000 because of the pulse-per-second reset);
     * in this example, the fragment is 0.5 ms late: 500'000 nanoseconds =
     * 300'000 (TTT) + 1'000'000'000 (second step) - 999'800'000 (glob. trigger)
     * and the plain difference,
     * 300'000 (TTT) - 999'800'000 (glob. trigger) = -999'500'000,
     * must be corrected by adding a whole second:
     */
    fragmentRelTime += 1'000'000'000;
  }
  
  waveformTime += nanoseconds::castFrom(fragmentRelTime);
  
  //
  // 4. correction for relative position of sample tagged by TTT in the buffer
  //
  waveformTime -= boardInfo.bufferLength;
  
  //
  // 5. correction for calibrated delays
  //
  /*
   * Waveform time has been expressed based on the "absolute" trigger time plus
   * an offset based on the Trigger Time Tag, which is synchronous with the
   * global trigger and reset every second.
   * We are missing a possible delay between the time of the trigger time scale
   * stepping into a new second and the the time TTT reset is effective.
   * 
   * 
   */
  
  waveformTime += boardInfo.TTTresetDelay;
  
  mf::LogTrace(fLogCategory) << "V1730 board '" << boardInfo.name
    << "' has data starting at electronics time " << waveformTime
    << " = " << fNominalTriggerTime << " (global trigger)"
    << " + " << nanoseconds(fragmentRelTime) << " (TTT - global trigger)"
    << " - " << boardInfo.bufferLength << " (buffer size)"
    << " + " << boardInfo.TTTresetDelay << " (reset delay)"
    ;
  
  return waveformTime;
  
  
} // icarus::DaqDecoderICARUSPMT::fragmentWaveformTimestampFromTTT()


//------------------------------------------------------------------------------
auto icarus::DaqDecoderICARUSPMT::neededBoardInfo
  (artdaq::Fragment::fragment_id_t fragment_id) const -> NeededBoardInfo_t
{
  
  assert(fBoardInfoLookup);

  /*
   * The trigger time is always the nominal one, because that is the
   * reference time of the whole DAQ (PMT, TPC...).
   * We only need to know how sooner than the trigger the V1730 buffer
   * starts. Oh, and the delay from the global trigger time to when
   * the readout board receives and processes the trigger signal.
   */
  daq::details::BoardInfoLookup::BoardInfo_t const* boardInfo
    = fBoardInfoLookup->findBoardInfo(fragment_id);
  if (!boardInfo) {
    if (fRequireKnownBoards) {
      cet::exception e("PMTDecoder");
      e << "Input fragment has ID " << fragment_id
        << " which has no associated board information (`BoardSetup`";
      if (!hasPMTconfiguration()) e << " + `.FragmentID`";
      throw e << ").\n";
    }
  }
  else {
    assert(boardInfo->fragmentID == fragment_id);
    assert(boardInfo->setup);
  }
  
  return fetchNeededBoardInfo(boardInfo, fragment_id);
} // icarus::DaqDecoderICARUSPMT::neededBoardInfo()


//------------------------------------------------------------------------------
unsigned int icarus::DaqDecoderICARUSPMT::extractTriggerTimeTag
  (artdaq::Fragment const& fragment)
{
  sbndaq::CAENV1730Fragment const V1730fragment { fragment };
  sbndaq::CAENV1730EventHeader const header = V1730fragment.Event()->Header;
  
  return { header.triggerTimeTag }; // prevent narrowing
  
} // icarus::DaqDecoderICARUSPMT::extractTriggerTimeTag()


//------------------------------------------------------------------------------
void icarus::DaqDecoderICARUSPMT::sortWaveforms
  (std::vector<raw::OpDetWaveform>& waveforms) const
{
  auto byChannelThenTime = []
    (raw::OpDetWaveform const& left, raw::OpDetWaveform const& right)
    {
      return (left.ChannelNumber() != right.ChannelNumber())
        ? left.ChannelNumber() < right.ChannelNumber()
        : left.TimeStamp() < right.TimeStamp();
    };
  
  std::sort(waveforms.begin(), waveforms.end(), byChannelThenTime);

} // icarus::DaqDecoderICARUSPMT::sortWaveforms()


//------------------------------------------------------------------------------
template <std::size_t NBits, typename T>
constexpr std::pair<std::array<std::size_t, NBits>, std::size_t>
icarus::DaqDecoderICARUSPMT::setBitIndices(T value) noexcept {
  
  std::pair<std::array<std::size_t, NBits>, std::size_t> res;
  auto& [ indices, nSetBits ] = res;
  for (std::size_t& index: indices) {
    index = (value & 1)? nSetBits++: NBits;
    value >>= 1;
  } // for
  return res;
  
} // icarus::DaqDecoderICARUSPMT::setBitIndices()


//------------------------------------------------------------------------------
void icarus::DaqDecoderICARUSPMT::initTrees
  (std::vector<std::string> const& treeNames)
{
  
  auto findTree = [](std::string const& name)
    {
      return static_cast<DataTrees>(
          std::distance(TreeNames.begin(),
          std::find(TreeNames.begin(), TreeNames.end(), name))
        );
    };
  
  for (std::string const& name: treeNames) {
    switch (findTree(name)) {
      case DataTrees::Fragments: initFragmentsTree(); break;
      case DataTrees::N:
      default:
        throw cet::exception("DaqDecoderICARUSPMT")
          << "initTrees(): no data tree supported with name '" << name
          << "'.\n";
    } // switch
  } // for names
  
} // icarus::DaqDecoderICARUSPMT::initTrees()


//------------------------------------------------------------------------------
void icarus::DaqDecoderICARUSPMT::initEventIDtree
  (TTree& tree, TreeData_EventID_t& data)
{
  
  usesEventInfo(); // this tree includes event information
  
  tree.Branch("run", &data.run);
  tree.Branch("subrun", &data.subrun);
  tree.Branch("event", &data.event);
  tree.Branch("timestamp", &data.timestamp.time, "timestamp/L");
  
} // icarus::DaqDecoderICARUSPMT::initEventIDtree()


//------------------------------------------------------------------------------
void icarus::DaqDecoderICARUSPMT::initFragmentsTree() {
  
  if (fTreeFragment) return;
  
  TTree* tree = art::ServiceHandle<art::TFileService>()
    ->make<TTree>("PMTfragments", "PMT fragment data");
  
  fTreeFragment = std::make_unique<TreeFragment_t>();
  fTreeFragment->tree = tree;
  auto& data = fTreeFragment->data;
  
  initEventIDtree(*tree, data);
  
  tree->Branch("fragmentID", &data.fragmentID);
  tree->Branch("fragTime", &data.fragTime.time, "fragTime/L"); // ROOT 6.24 can't detect 64-bit
  tree->Branch("fragTimeSec", &data.fragTime.split.seconds);
  tree->Branch("TTT", &data.TriggerTimeTag, "TTT/l"); // ROOT 6.24 can't detect 64-bit
  tree->Branch("trigger", &data.trigger.time, "trigger/L"); // ROOT 6.24 can't detect 64-bit
  tree->Branch("triggerSec", &data.trigger.split.seconds);
  tree->Branch("triggerNS", &data.trigger.split.nanoseconds);
  
} // icarus::DaqDecoderICARUSPMT::initFragmentsTree()


//------------------------------------------------------------------------------
void icarus::DaqDecoderICARUSPMT::usesEventInfo() {
  
  // the allocation of fEventInfo is the flag for event information usage
  if (!fEventInfo) fEventInfo = std::make_unique<TreeData_EventID_t>();
  
} // icarus::DaqDecoderICARUSPMT::usesEventInfo()


void icarus::DaqDecoderICARUSPMT::assignEventInfo
  (TreeData_EventID_t& treeData) const
{
  
  assert(fEventInfo);
  treeData = *fEventInfo; // nice slicing
  
} // icarus::DaqDecoderICARUSPMT::assignEventInfo()


//------------------------------------------------------------------------------
void icarus::DaqDecoderICARUSPMT::fillTreeEventID
  (art::Event const& event, TreeData_EventID_t& treeData) const
{
  art::EventID const& id = event.id();
  treeData.run    = id.run();
  treeData.subrun = id.subRun();
  treeData.event  = id.event();
  
  art::Timestamp const& timestamp = event.time();
  treeData.timestamp
    = { static_cast<int>(timestamp.timeHigh()), timestamp.timeLow() };
  
} // icarus::DaqDecoderICARUSPMT::fillTreeEventID()


//------------------------------------------------------------------------------
DEFINE_ART_MODULE(icarus::DaqDecoderICARUSPMT)


//------------------------------------------------------------------------------
