#include "LTrackerCalibrationManager.hh"
#include "LTrackerTools.hh"



#include <iostream>

LTrackerCalibrationManager::LTrackerCalibrationManager()
{
    calRunFile = 0;
    InitialTargetRun = -1;
    FinalTargetRun = -1;
    verboseFLAG = true;
}


LTrackerCalibrationManager& LTrackerCalibrationManager::GetInstance()
{
    static LTrackerCalibrationManager instance; // Guaranteed to be destroyed.
    // Instantiated on first use.
    return instance;
}


int LTrackerCalibrationManager::LoadRun (const char *fileInp)
{
    // Reset previous loaded runs
    if (calRunFile) {
        if (calRunFile->IsOpen() ) calRunFile->Close();
        calRunFile = 0;
        // Sure we want also to reset the target runs? Today we reply yes... check
        InitialTargetRun = -1;
        FinalTargetRun = -1;
    }
    calRunFile = new LEvRec0File (fileInp);
    if (!calRunFile || ! (calRunFile->IsOpen() ) ) {
        std::cerr << "Error! Attempt to load a calibration run, but the file cannot be opened."
                  << std::endl;
        calRunFile = 0;
        return -999;
    }
    return 0;
}

void LTrackerCalibrationManager::SetTargetRuns (const int InitialRun, const int FinalRun)
{
    InitialTargetRun = InitialRun;
    FinalTargetRun = FinalRun;
    return;
}

LTrackerCalibration* LTrackerCalibrationManager::Calibrate (const int nEvents, const int skipEvents)
{
    if (calRunFile == 0 || ! (calRunFile->IsOpen() ) ) {
        std::cerr << "Error! Attempt to call the \"Calibrate\" method, but no calibration run loaded."
                  << std::endl;
        return 0;
    }
    int nEntries = calRunFile->GetEntries();
    if (nEntries < skipEvents + nEvents) return 0;
    std::vector<int> pivot = CalculateCalibrationSlots (nEvents, skipEvents, nEntries);
    LTrackerCalibration *result = CreateTrackerCalibration();
    for (uint is = 0; is < pivot.size() - 1; ++is) result->Add (CalibrateSlot (pivot[is], pivot[is + 1]) );
    return result;
}

std::vector<int> LTrackerCalibrationManager::CalculateCalibrationSlots (const int nEvents, const int skipEvents, const int nEntries)
{
    int nEv = nEvents;
    int skipEv = skipEvents;
    if (nEvents == -1) nEv = nEntries;
    if (skipEvents == -1) skipEv = 0;
    const int nSlots = (nEntries - skipEv) / nEv;
    std::vector<int> pivot (nSlots + 1);
    for (int is = 0; is < nSlots + 1; ++is) pivot[is] = skipEv + is * nEv;
    return pivot;
}

LTrackerCalibrationSlot* LTrackerCalibrationManager::CalibrateSlot (const int StartEntry, const int StopEntry)
{
    if (calRunFile == 0 || ! (calRunFile->IsOpen() ) ) {
        std::cerr << "Error! Attempt to call the \"CalibrateSlot\" method, but no calibration run loaded."
                  << std::endl;
        return 0;
    }
    LTrackerSlotCalibrator slotcal(calRunFile);
    slotcal.SetEntries(StartEntry, StopEntry);
    std::vector<statq> statraw   = slotcal.RawMeanSigma();
    std::vector<statq> statclean = slotcal.CleanedMeanSigma ();
    LTrackerMask CN_mask         = slotcal.ComputeCNMask ();
    std::vector<statq> statCNcorr= slotcal.CNCorrectedSigma ();
    std::vector<double> ngindex  = slotcal.GaussianityIndex ();

    // Start and stop events
    LEvRec0 cev;
    calRunFile->SetTheEventPointer (cev);
    calRunFile->GetEntry (StartEntry);
    int StartEvent = static_cast<int> (cev.event_index);
    calRunFile->GetEntry (StopEntry - 1);
    int StopEvent = static_cast<int> (cev.event_index);
    // Result

    unsigned int ncor=statCNcorr.size();
    double ped[ncor];
    double sig[ncor];
    for (int i=0; i<ncor; i++) {
      ped[i]=statCNcorr[i].GetMean();
      sig[i]=statCNcorr[i].GetStdDev();
    }

    bool cnm[NCHAN];
    for (int ichan=0; ichan<NCHAN; ichan++) cnm[ichan]=CN_mask[ichan]; // vector<bool> does not have c_array functionality :(
    LTrackerCalibrationSlot *result = new LTrackerCalibrationSlot
      (StartEvent, StopEvent, ped, sig, &ngindex[0], cnm);
    return result;
}







LTrackerCalibrationManager::~LTrackerCalibrationManager()
{
    // do not care about singleton destructor
}


LTrackerCalibration* LTrackerCalibrationManager::CreateTrackerCalibration()
{
    if (calRunFile == 0 || ! (calRunFile->IsOpen() ) ) {
        std::cerr << "Error! Attempt to create a tracker calibration but no calibration run loaded."
                  << std::endl;
        return 0;
    }
    int RunId = calRunFile->GetRunId();
    if (InitialTargetRun == -1 || FinalTargetRun == -1) {
        std::cerr << "Warning! Target run interval of current calibration not defined."
                  << std::endl;
    }
    LTrackerCalibration *result = new LTrackerCalibration (RunId, InitialTargetRun, FinalTargetRun);
    return result;
}
