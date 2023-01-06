////////////////////////////////////////////////////////////////////////
// Class:       FilterCRTPMTMatching
// Plugin Type: analyzer (Unknown Unknown)
// File:        FilterCRTPMTMatching_module.cc
//
// Generated at Thu Feb 17 13:31:00 2022 by Biswaranjan Behera using cetskelgen
// from  version .
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
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"

// Data product includes
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "icaruscode/Decode/DataProducts/TriggerConfiguration.h"
// C++ includes
#include <vector>
#include <map>
#include <numeric>

// ROOT includes
#include "TTree.h"
#include "TVector3.h"

namespace icarus
{
    namespace crt
    {
        class FilterCRTPMTMatching;
        struct CRTPMT
        {
            double tof;
            sbn::crt::CRTHit const* CRTHit=nullptr;
        };
        struct MatchedCRT
        {
            std::vector<CRTPMT> entering;
            std::vector<CRTPMT> exiting;
	    /**
            * Type of the matching.
            * 
            * Secret decoder ring:
            *  * `0`: no matched CRT
            *  * `1`: particle entering from Top CRT
            *  * `2`: particle entering from Side CRT
            *  * `3`: particle entering from Top and exiting from Side CRT
            *  * `4`: particle exiting/coincidence with Top CRT
            *  * `5`: particle exiting/coincidence with Side CRT
            *  * `6`: multiple particles entering from Top and Side
            *  * `7`: multiple particles entering from Top and Side and exiting from side
            *  * `8`: all other cases
            * 
            */
	    int matchType;
        };
        struct TrajPoint
        {
            double X;
            double Y;
            double Z;
            double T;
        };
    }
}

using namespace icarus::crt;

bool flashInTime(double const &flashTime, int gateType, double gateDiff, double gateWidth)
{
    //   {1, 2.2},  // BNB
    //   {2, 10.1}, // NuMI
    //   {3, 2.2},  // BNB Offbeam
    //   {4, 10.1}, // NuMI Offbeam
   

    //As a reminder, I will leave here the commented part of the vetoOffset, in case something changes in the future
    /*double vetoOffset = 3.5;*/
    double activeGate = gateWidth /*- vetoOffset*/;

    double relFlashTime = flashTime + gateDiff / 1000. /*- vetoOffset*/;
    mf::LogInfo("FilterCRTPMTMatching::flashInTime") << "Gate Diff " << gateDiff / 1000 << " Ftime+gateDiff " << flashTime + gateDiff / 1000. << " " << activeGate;

    return ((relFlashTime > 0) && (relFlashTime < activeGate));
}

icarus::crt::MatchedCRT CRTHitmatched(double flashTime, std::vector<sbn::crt::CRTHit>const crtHits, double interval)
{

    std::vector<icarus::crt::CRTPMT> enteringCRTHits;
    std::vector<icarus::crt::CRTPMT> exitingCRTHits;

    int hasCRTHit = 0;
    int topen = 0, topex = 0, sideen = 0, sideex = 0;
    for (auto const& crtHit : crtHits)
    {
        double tof = crtHit.ts1_ns / 1e3 - flashTime;
        // std::cout<<"TOF "<<tof<<" "<<crtHit->ts1_ns<<" "<<crtHit->plane<<std::endl;
        if (tof < 0 && abs(tof) < interval)
        {
	// crtHit->Plane >36 are Side CRT, crtHit->Plane <36 are Top CRT. I will update this with an existing function: isSideCRT
            if (crtHit.plane > 36)
                sideen++;
            else
                topen++;
	    CRTPMT m_match = {tof, &crtHit};
            enteringCRTHits.emplace_back(m_match);
        }
        else if (tof >= 0 && abs(tof) < interval)
        {
            if (crtHit.plane > 36)
                sideex++;
            else
                topex++;
            CRTPMT m_match = {tof, &crtHit};
            exitingCRTHits.emplace_back(m_match);
        }
    }

    // Secret decoder ring, do not remove!
    // Coup attempted. Unfortunatly it failed. Secrecy will live on, until this is approved, then I will move it to a readme/online documentation
    // hasCRTHit = 0, no matched CRT
    // hasCRTHit = 1, Muon entering from Top CRT
    // hasCRTHit = 2, Muon entering from Side CRT
    // hasCRTHit = 3, Muon entering from Top and exiting from Side CRT
    // hasCRTHit = 4, Muon exiting/coincidence with Top CRT
    // hasCRTHit = 5, Muon exiting/coincidence with Side CRT
    // hasCRTHit = 6, Multiple muons entering from Top and Side
    // hasCRTHit = 7, Multiple muons entering from Top and Side and exiting from side
    // hasCRTHit = 8, all other cases
    if (topen == 0 && sideen == 0 && topex == 0 && sideex == 0)
        hasCRTHit = 0;
    else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 0)
        hasCRTHit = 1;
    else if (topen == 0 && sideen == 1 && topex == 0 && sideex == 0)
        hasCRTHit = 2;
    else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 1)
        hasCRTHit = 3;
    else if (topen == 0 && sideen == 0 && topex == 1 && sideex == 0)
        hasCRTHit = 4;
    else if (topen == 0 && sideen == 0 && topex == 0 && sideex == 1)
        hasCRTHit = 5;
    else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex == 0)
        hasCRTHit = 6;
    else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex >= 1)
        hasCRTHit = 7;
    else
        hasCRTHit = 8;

    icarus::crt::MatchedCRT matches{enteringCRTHits, exitingCRTHits, hasCRTHit};

    return matches;
}

class icarus::crt::FilterCRTPMTMatching : public art::EDFilter
{
public:
    using CRTHit = sbn::crt::CRTHit;

    explicit FilterCRTPMTMatching(fhicl::ParameterSet const &p);

    // Required functions.
    bool filter(art::Event &) override;

private:
    // Declare member data here.

    static bool HitCompare(const art::Ptr<CRTHit> &h1, const art::Ptr<CRTHit> &h2);
    void ClearVecs();

    std::vector<art::InputTag> fOpFlashModuleLabelVec;
    art::InputTag              fCrtHitModuleLabel;
    art::InputTag              fTriggerLabel;
    art::InputTag	       fTriggerConfiguration;

    //double                     fFlashTimeCut;

    int                        fEvent;  ///< number of the event being processed.
    int                        fRun;    ///< number of the run being processed.
    int                        fSubRun; ///< number of the sub-run being processed.

    // add trigger data product vars
    unsigned int               m_gate_type;
    std::string                m_gate_name;
    uint64_t                   m_trigger_timestamp;
    uint64_t                   m_gate_start_timestamp;
    uint64_t                   m_trigger_gate_diff;
    uint64_t		       m_gate_width;

    std::string		       fFilterLevel;		// Filter level, default is loose 
    double		       fOpHitAmplitude;    	// ADC threshold for the PMT
    int                        fnOpHitToTrigger;        // Number of OpHit above threshold to mimic the trigger
    double 		       fTimeOfFlightInterval;   // Time of Flight interval to find the match
    bool                       fOutputTree;             // Output tree or not
    map<int, art::InputTag>    fFlashLabels;

    TTree*                     fMatchTree;

    // matchTree vars
    vector<double>             fOpFlashPos;             ///< Position of the optical Flash Barycenter.
    double                     fOpFlashTime;            ///< Time of the optical Flash w.r.t. Global Trigger.
    bool                       fInTime;                 ///< Was the OpFlash in beam spill time?
    double                     fOpFlashPE;              ///< Total PEs of the optical Flash.
    vector<double>             fOpHitX;                 ///< X position of the OP hit.
    vector<double>             fOpHitY;                 ///< Y position of the OP hit.
    vector<double>             fOpHitZ;                 ///< Z position of the OP hit.
    vector<double>             fOpHitT;                 ///< T of the OP hit.
    vector<double>             fOpHitA;                 ///< Amplitude (ADC) of the OP hit.
    int                        fNCRTmatch;
    vector<vector<double>>     fMatchedCRTpos;          ///< Position of the matched CRTs.
    vector<double>             fMatchedCRTtime;         ///< Time of the matched CRT Hits w.r.t. Global Trigger.
    vector<int>                fMatchedCRTregion;       ///< Region of the matched CRT Hits.
    vector<vector<int>>        fMatchedCRTmodID;        ///< ModID of the CRT module.
    vector<int>                fMatchedCRTsys;          ///< Subsystem of the matched CRT Hit.
    vector<double>             fMatchedCRTamplitude;    ///< Amplitude of the matched CRT Hit.
    vector<double>             fTofOpFlash;             ///< Time of Flight between matched CRT and Optical Flash.
    vector<double>             fCRTGateDiff;            ///< Difference between CRTHit and BeamGate opening.
    int                        fEventType;              ///< Was triggered the event?

    geo::GeometryCore const*   fGeometryService; ///< pointer to Geometry provider.
};

icarus::crt::FilterCRTPMTMatching::FilterCRTPMTMatching(fhicl::ParameterSet const &p) : EDFilter{p} 
    ,fOpFlashModuleLabelVec(p.get<std::vector<art::InputTag>>("OpFlashModuleLabelVec", {"opflashE","opflashW"}))
    ,fCrtHitModuleLabel(p.get<art::InputTag>("CrtHitModuleLabel", "crthit"))
    ,fTriggerLabel(p.get<art::InputTag>("TriggerLabel", "daqTrigger")) 
    ,fTriggerConfiguration(p.get<art::InputTag>("TriggerConfiguration", "triggerconfig"))
    ,fFilterLevel(p.get<std::string>("FilterLevel","FilterOption"))
    ,fOpHitAmplitude(p.get<double>("OpHitAmplitude", 0.))
    ,fnOpHitToTrigger(p.get<int>("nOpHitToTrigger", 0))
    ,fTimeOfFlightInterval(p.get<double>("TimeOfFlightInterval", 0.))
    ,fOutputTree(p.get<bool>("OutputTree","treeOption"))
    ,fGeometryService(lar::providerFrom<geo::Geometry>())
{

    art::ServiceHandle<art::TFileService> tfs;
    if (fOutputTree==true) {

    fMatchTree = tfs->make<TTree>("matchTree", "CRTHit - OpHit/Flash matching analysis");
    fMatchTree->Branch("event", &fEvent, "event/I");
    fMatchTree->Branch("run", &fRun, "run/I");
    fMatchTree->Branch("subrun", &fSubRun, "subrun/I");
    fMatchTree->Branch("opFlash_pos", &fOpFlashPos);
    fMatchTree->Branch("opFlash_time", &fOpFlashTime);
    fMatchTree->Branch("opHit_x", &fOpHitX);
    fMatchTree->Branch("opHit_y", &fOpHitY);
    fMatchTree->Branch("opHit_z", &fOpHitZ);
    fMatchTree->Branch("opHit_t", &fOpHitT);
    fMatchTree->Branch("opHit_amplitude", &fOpHitA);
    fMatchTree->Branch("inTime", &fInTime);
    fMatchTree->Branch("CRT_matches", &fNCRTmatch);
    fMatchTree->Branch("CRT_pos", &fMatchedCRTpos);
    fMatchTree->Branch("CRT_Time", &fMatchedCRTtime);
    fMatchTree->Branch("CRT_region", &fMatchedCRTregion);
    fMatchTree->Branch("CRT_FEB", &fMatchedCRTmodID);
    fMatchTree->Branch("CRT_subsys", &fMatchedCRTsys);
    fMatchTree->Branch("CRT_pes", &fMatchedCRTamplitude);
    fMatchTree->Branch("TOF_opflash", &fTofOpFlash);
    fMatchTree->Branch("CRT_gate_diff", &fCRTGateDiff);
    fMatchTree->Branch("eventType", &fEventType);
    fMatchTree->Branch("gate_type", &m_gate_type, "gate_type/b");
    fMatchTree->Branch("gate_name", &m_gate_name);
    fMatchTree->Branch("trigger_timestamp", &m_trigger_timestamp, "trigger_timestamp/l");
    fMatchTree->Branch("gate_start_timestamp", &m_gate_start_timestamp, "gate_start_timestamp/l");
    fMatchTree->Branch("trigger_gate_diff", &m_trigger_gate_diff, "trigger_gate_diff/l");
    }
}

bool icarus::crt::FilterCRTPMTMatching::filter(art::Event &e)
{
    mf::LogDebug("FilterCRTPMTMatching: ") << "beginning analyis";

    // Start by fetching some basic event information for our n-tuple.
    fEvent = e.id().event();
    fRun = e.run();
    fSubRun = e.subRun();

    ClearVecs();

    // add trigger info
    auto const& trigger_handle = e.getProduct<sbn::ExtraTriggerInfo>(fTriggerLabel);
    sbn::triggerSource bit = trigger_handle.sourceType;
    m_gate_type = value(bit);
    m_gate_name = bitName(bit);
    m_trigger_timestamp = trigger_handle.triggerTimestamp;
    m_gate_start_timestamp = trigger_handle.beamGateTimestamp;
    m_trigger_gate_diff = trigger_handle.triggerTimestamp - trigger_handle.beamGateTimestamp;

    //Read Beam Gate Size
    art::Run const& run = e.getRun();
    auto const& trigger_configuration = run.getProduct<icarus::TriggerConfiguration>(fTriggerConfiguration);
    m_gate_width=trigger_configuration.getGateWidth(m_gate_type);

    // CRTHits

    auto const& crtHitList = e.getProduct<std::vector<CRTHit>>(fCrtHitModuleLabel);

    bool Cosmic = false;
    int  Type   = 9;

    if (fFilterLevel != "loose" || fFilterLevel != "tight" || fFilterLevel != "medium") fFilterLevel="loose";

    for (auto const& flashLabel : fOpFlashModuleLabelVec)
    {
        art::Handle<std::vector<recob::OpFlash>> flashHandle;

        if (!e.getByLabel(flashLabel, flashHandle))
        {
            mf::LogError("FilterCRTPMTMatching") << "Did not find opflash object with label: " << flashLabel;
            continue;
        }

        std::vector<art::Ptr<recob::OpFlash>> opFlashList;

        art::fill_ptr_vector(opFlashList, flashHandle);

        // This is returning only ophit associations in the first TPC? Is that what is wanted?
//        art::FindManyP<recob::OpHit> findManyHits(flashHandles[flashList.first], e, fFlashLabels[flashList.first]);
        art::FindManyP<recob::OpHit> findManyHits(flashHandle, e, flashLabel);

        for (size_t iflash = 0; iflash < opFlashList.size(); iflash++)
        {

            int eventType = 0;
            auto const &flash = opFlashList[iflash];

            double tflash = flash->Time();

            vector<art::Ptr<recob::OpHit>> hits = findManyHits.at(iflash);
            int nOpHitsTriggering = 0;
            double firstTime = 999999;
            double flash_pos[3] = {0, 0, 0};
            double ampsum = 0, t_m = 0;
            for (auto const &hit : hits)
            {
                if (hit->Amplitude() > fOpHitAmplitude)
                    nOpHitsTriggering++;
                if (firstTime > hit->PeakTime())
                    firstTime = hit->PeakTime();
                geo::Point_t pos = fGeometryService->OpDetGeoFromOpChannel(hit->OpChannel()).GetCenter();
                double amp = hit->Amplitude();
                ampsum += amp;
                fOpHitX.push_back(pos.X());
                fOpHitY.push_back(pos.Y());
                fOpHitZ.push_back(pos.Z());
                fOpHitT.push_back(hit->PeakTime());
                fOpHitA.push_back(amp);
                flash_pos[0] = flash_pos[0] + pos.X() * amp;
                flash_pos[1] = flash_pos[1] + pos.Y() * amp;
                flash_pos[2] = flash_pos[2] + pos.Z() * amp;
                t_m = t_m + hit->PeakTime();
            }
            flash_pos[0] = flash_pos[0] / ampsum;
            flash_pos[1] = flash_pos[1] / ampsum;
            flash_pos[2] = flash_pos[2] / ampsum;
            t_m = t_m / nOpHitsTriggering;
            if (nOpHitsTriggering < fnOpHitToTrigger)
            {
                continue;
            }
            bool   inTime   = flashInTime(firstTime, m_gate_type, m_trigger_gate_diff, m_gate_width);

            mf::LogTrace("FilterCRTPMTMatching") << "\nFlash Time " << tflash << "  First Op Hit " << firstTime << "  " << inTime << "\n" <<
                        "Average Pos X " << flash_pos[0] << "  Y " << flash_pos[1] << "  Z " << flash_pos[2] << " nPMT " << nOpHitsTriggering << " " << ampsum << "\n" <<
                        "Flash X " << flash->XCenter() << "  " << flash->YCenter() << "  " << flash->ZCenter();
            // Now get the CRT. Search the CRT Hits within -100 from the flash time
            // NB currently the selection uses only top CRT.
            //  for the future a differentiation between Top and Side and some considerations based
            //  on the proximity of the light are necessary
            // std::vector< art::Ptr<CRTHit> > selCRTHits;
            icarus::crt::MatchedCRT CRTmatches = CRTHitmatched(firstTime, crtHitList, fTimeOfFlightInterval);
            auto nCRTHits = CRTmatches.entering.size() + CRTmatches.exiting.size();

	    {	        
                // The following meant to try to preserve the original formatting
                mf::LogDebug log("FilterCRTPMTMatching");
                log << "Matched CRT " << nCRTHits << "  entering: \n";

                for (auto &entering : CRTmatches.entering)
                {
                    log << "TOF " << entering.tof <<  " Region " << entering.CRTHit->plane << "\n";
                }

                log << "Exiting: \n";

                for (auto &exiting : CRTmatches.exiting)
                {
                    log << "TOF " << exiting.tof << " Region " << exiting.CRTHit->plane << "\n";
                }
	    }
            bool   matched = false;
            double peflash = flash->TotalPE();

            if (nCRTHits > 0) matched = true;

            if (matched == true)
            {
                eventType = CRTmatches.matchType;

                for (auto &entering : CRTmatches.entering)
                {
                    vector<double> CRTpos     = {entering.CRTHit->x_pos, entering.CRTHit->y_pos, entering.CRTHit->z_pos};
                    double         CRTtime    = entering.CRTHit->ts1_ns / 1e3;
                    ULong64_t      CRTAbsTime = entering.CRTHit->ts0_s;
                    int            CRTRegion  = entering.CRTHit->plane;
                    int            CRTSys     = 0;

                    if (CRTRegion >= 36) CRTSys = 1;

                    double CRTPe = entering.CRTHit->peshit;
                    double CRTTof_opflash = CRTtime - tflash;
                    double CRTGateDiff    = CRTAbsTime - m_gate_start_timestamp;

                    std::vector<int> HitFebs;
                    for (auto &crts : entering.CRTHit->feb_id)
                    {
                        HitFebs.emplace_back((int)crts);
                    }

                    fMatchedCRTmodID.emplace_back(HitFebs);
                    fMatchedCRTpos.emplace_back(CRTpos);
                    fMatchedCRTtime.emplace_back(CRTtime);
                    fTofOpFlash.emplace_back(CRTTof_opflash);
                    fMatchedCRTregion.emplace_back(CRTRegion);
                    fMatchedCRTsys.emplace_back(CRTSys);
                    fMatchedCRTamplitude.emplace_back(CRTPe);
                    fCRTGateDiff.emplace_back(CRTGateDiff);
                    HitFebs.clear();
                }

                for (auto &exiting : CRTmatches.exiting)
                {
                    vector<double> CRTpos     = {exiting.CRTHit->x_pos, exiting.CRTHit->y_pos, exiting.CRTHit->z_pos};
                    double         CRTtime    = exiting.CRTHit->ts1_ns / 1e3;
                    ULong64_t      CRTAbsTime = exiting.CRTHit->ts0_s;
                    int            CRTRegion  = exiting.CRTHit->plane;
                    int            CRTSys     = 0;

                    if (CRTRegion >= 36) CRTSys = 1;

                    double CRTPe          = exiting.CRTHit->peshit;
                    double CRTTof_opflash = CRTtime - tflash;
                    double CRTGateDiff    = CRTAbsTime - m_gate_start_timestamp;

                    std::vector<int> HitFebs;
                    for (auto &crts : exiting.CRTHit->feb_id)
                    {
                        HitFebs.emplace_back((int)crts);
                    }

                    fMatchedCRTmodID.emplace_back(HitFebs);
                    fMatchedCRTpos.emplace_back(CRTpos);
                    fMatchedCRTtime.emplace_back(CRTtime);
                    fMatchedCRTregion.emplace_back(CRTRegion);
                    fMatchedCRTsys.emplace_back(CRTSys);
                    fMatchedCRTamplitude.emplace_back(CRTPe);
                    fTofOpFlash.emplace_back(CRTTof_opflash);
                    fCRTGateDiff.emplace_back(CRTGateDiff);
                    HitFebs.clear();
                    // Filla i vettori
                }
                if ((eventType == 1 || eventType == 3 || eventType == 6 || eventType == 7) && inTime == true) Cosmic = true;

                if (inTime == true) Type = eventType;
            } // if matched

            fOpFlashPE = peflash;
            fOpFlashPos = {flash_pos[0], flash_pos[1], flash_pos[2]};
            fOpFlashTime = firstTime;
            fInTime = inTime;
            fEventType = eventType;
            fNCRTmatch = nCRTHits;
 	    if (fOutputTree==true) fMatchTree->Fill();
            fOpHitX.clear();
            fOpHitY.clear();
            fOpHitZ.clear();
            fOpHitT.clear();
            fOpHitA.clear();
            fOpFlashPos.clear();
            fMatchedCRTmodID.clear();
            fMatchedCRTpos.clear();
            fMatchedCRTtime.clear();
            fMatchedCRTregion.clear();
            fMatchedCRTsys.clear();
            fMatchedCRTamplitude.clear();
            fTofOpFlash.clear();
            fCRTGateDiff.clear();
        } // for Flash
    }
    fEventType   = Type;

    if (fFilterLevel == "loose") Cosmic = false;  // By default, with loose Filtering, everything is "intersting", nothing tagged as clear cosmic
    else if (fFilterLevel == "medium") {
	if (fEventType == 1 || fEventType == 3 || fEventType == 6 || fEventType == 7) Cosmic = true;
	// With Medium filter, everything (inTime) which is associated with Top CRT Hit before the Flash is filtered as clear cosmic
	else Cosmic = false;
    }
    else if (fFilterLevel == "thight"){
	if (fEventType != 0 || fEventType != 4 || fEventType != 5) Cosmic = true;
	// With Thight filter, everything (inTime) associated with a CRT Hit before the Flash is filtered as clear cosmic
	else Cosmic = false;
    }
    return Cosmic;
}

void icarus::crt::FilterCRTPMTMatching::ClearVecs()
{
    // matchTree
    fOpHitX.clear();
    fOpHitY.clear();
    fOpHitZ.clear();
    fOpHitT.clear();
    fOpHitA.clear();
    fOpFlashPos.clear();
    fMatchedCRTmodID.clear();
    fMatchedCRTpos.clear();
    fMatchedCRTtime.clear();
    fMatchedCRTregion.clear();
    fMatchedCRTsys.clear();
    fMatchedCRTamplitude.clear();
    fTofOpFlash.clear();
    fCRTGateDiff.clear();
}

bool icarus::crt::FilterCRTPMTMatching::HitCompare(const art::Ptr<CRTHit> &hit1, const art::Ptr<CRTHit> &hit2)
{

    if (hit1->ts1_ns != hit2->ts1_ns)
        return false;
    if (hit1->plane != hit2->plane)
        return false;
    if (hit1->x_pos != hit2->x_pos)
        return false;
    if (hit1->y_pos != hit2->y_pos)
        return false;
    if (hit1->z_pos != hit2->z_pos)
        return false;
    if (hit1->x_err != hit2->x_err)
        return false;
    if (hit1->y_err != hit2->y_err)
        return false;
    if (hit1->z_err != hit2->z_err)
        return false;
    if (hit1->tagger != hit2->tagger)
        return false;

    return true;
}

DEFINE_ART_MODULE(icarus::crt::FilterCRTPMTMatching)
