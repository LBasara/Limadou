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

LTrackerCalibrationSlot* LTrackerCalibrationManager::CalibrateSlot (const int startEntry, const int stopEntry)
{
    if (calRunFile == 0 || ! (calRunFile->IsOpen() ) ) {
        std::cerr << "Error! Attempt to call the \"CalibrateSlot\" method, but no calibration run loaded."
                  << std::endl;
        return 0;
    }
    this->StartEntry=startEntry;
    this->StopEntry =stopEntry;
    // RawMeanSigma
    std::vector<statq> statraw = RawMeanSigma ();
    // First cleaning
    std::vector<statq> statclean = CleanedMeanSigma (statraw);
    // Compute CN mask
    std::vector<bool> CN_mask = ComputeCNMask (statclean);
    // CNCorrectedSigma
    std::vector<statq> statCNcorr =CNCorrectedSigma (statclean, CN_mask);
    // Gaussianity
    std::vector<double> ngindex=GaussianityIndex (statCNcorr, CN_mask);
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




std::vector<bool> LTrackerCalibrationManager::ComputeCNMask (std::vector<statq> cleanstat)
{


    std::vector<bool> CN_mask (NCHAN);
    // Create "histograms" of sigmas per VA and set them to zero
    int hSigma1[N_VA][NSIGMA1BIN + 1]={{0}}; // overflow

    double xBin[NSIGMA1BIN + 1];
    for (int iBin = 0; iBin < NSIGMA1BIN + 1; ++iBin)
      xBin[iBin] = MINSIGMA1 + (MAXSIGMA1 - MINSIGMA1) * iBin / static_cast<double> (NSIGMA1BIN);
    // Fill histograms
    for (int iChan = 0; iChan < NCHAN; ++iChan) {
        int iVA = ChanToVA (iChan);
        for (int iBin = NSIGMA1BIN; iBin > -1; --iBin) {
            if (cleanstat[iChan].GetStdDev() > xBin[iBin]) {
                ++hSigma1[iVA][iBin];
                break;
            }
        }
    }
    // Find histos' maxima
    double hMaxima[N_VA];
    for (int iVA = 0; iVA < N_VA; ++iVA) {
        int tmpMax = -99999;
        hMaxima[iVA] = -99999.;
        for (int iBin = 0; iBin < NSIGMA1BIN; ++iBin) {
            if (hSigma1[iVA][iBin] > tmpMax) {
                hMaxima[iVA] = xBin[iBin] + 0.5 * (MAXSIGMA1 - MINSIGMA1) / NSIGMA1BIN;
                tmpMax = hSigma1[iVA][iBin];
            }
        }
    }
    // Compute CN mask
    for (int iChan = 0; iChan < NCHAN; ++iChan) {
        int iVA = ChanToVA (iChan);
        CN_mask[iChan] = std::fabs (cleanstat[iChan].GetStdDev() - hMaxima[iVA]) <= HALFSIGMA1WIDTH;
    }
    if (verboseFLAG) {
        std::cout << "CNmask computed" << std::endl;
    }
    return CN_mask;
}





std::vector<statq> LTrackerCalibrationManager::RawMeanSigma ()
{
    std::vector<statq> rawstat;
    LEvRec0 cev;
    calRunFile->SetTheEventPointer (cev);
    std::vector<std::vector<float>> matchan (NCHAN);
    for (int iEntry = StartEntry; iEntry < StopEntry; iEntry++) {
        calRunFile->GetEntry (iEntry);
        for (int iChan = 0; iChan < NCHAN; ++iChan) {
            float x = static_cast<float> (cev.strip[iChan]);
                  matchan[iChan].push_back(x);
        }
    }
    for (auto v : matchan)
        rawstat.push_back (statq (v) );
    if (verboseFLAG) std::cout << "RawMeanSigma computed" << std::endl;
    return rawstat;
}






std::vector<statq> LTrackerCalibrationManager::CleanedMeanSigma (std::vector<statq> rawstat)
{
    std::vector<statq> cleanstat;
    LEvRec0 cev;
    calRunFile->SetTheEventPointer (cev);
    // Average counts and squares
    std::vector<std::vector<float>> matchan (NCHAN);
    for (int iEntry = StartEntry; iEntry < StopEntry; ++iEntry) {
        calRunFile->GetEntry (iEntry);
        for (int iChan = 0; iChan < NCHAN; ++iChan) {
            float x = static_cast<float> (cev.strip[iChan]);
            float absdiff = std::fabs (x - rawstat[iChan].GetMean() );
            float cleanthres = CHANCLEANINGTHRESHOLD * rawstat[iChan].GetStdDev();
            if ( absdiff > cleanthres) continue;
            matchan[iChan].push_back (x);
        }
    }
    for (auto v : matchan) cleanstat.push_back (statq (v) );
    if (verboseFLAG) std::cout << "CleanedMeanSigma computed" << std::endl;
    return cleanstat;
}






std::vector<statq>  LTrackerCalibrationManager::CNCorrectedSigma (std::vector<statq> statclean, const std::vector<bool> CN_mask)
{

  std::vector<statq> statCNcorr;

    LEvRec0 cev;
    calRunFile->SetTheEventPointer (cev);


    // Average counts and squares

    std::vector<std::vector<float>> matchan (NCHAN);

    for (int iEntry = StartEntry; iEntry < StopEntry; ++iEntry) {
        calRunFile->GetEntry (iEntry);
        std::vector<double> CN=ComputeCN (cev.strip, statclean, CN_mask);
        for (int iChan = 0; iChan < NCHAN; ++iChan) {
            double x = static_cast<double> (cev.strip[iChan]);
            double diff = (x - statclean[iChan].GetMean());
            if (std::fabs (diff) > CHANCLEANINGTHRESHOLD * statclean[iChan].GetStdDev()) continue;
            double y = (x - CN[ChanToVA (iChan)]); // CN corrected!
            matchan[iChan].push_back (y);
        }
    }

    for (auto v : matchan)
        statCNcorr.push_back (statq (v) );

    if (verboseFLAG) std::cout << "CNCorrectedSigma computed" << std::endl;
    return statCNcorr;
}






std::vector<double> LTrackerCalibrationManager::GaussianityIndex (std::vector<statq> statCNcorr, const std::vector<bool> CN_mask)
{

  std::vector<double> ngindex(NCHAN);

    LEvRec0 cev;
    calRunFile->SetTheEventPointer (cev);
    // Gaussianity index
    int ngcounter[NCHAN];
    for (int iChan = 0; iChan < NCHAN; ++iChan) {
        ngindex[iChan] = 0.;
        ngcounter[iChan] = 0;
    }
    for (int iEntry = StartEntry; iEntry < StopEntry; ++iEntry) {
        calRunFile->GetEntry (iEntry);
        std::vector<double> CN=ComputeCN (cev.strip, statCNcorr, CN_mask);
        for (int iChan = 0; iChan < NCHAN; ++iChan) {
            double x = (static_cast<double> (cev.strip[iChan]) - statCNcorr[iChan].GetMean() - CN[ChanToVA (iChan)]);
            if (std::fabs (x) > GAUSSIANITYSIGMATHRESHOLD * statCNcorr[iChan].GetStdDev()) ++ngindex[iChan];
            ++ngcounter[iChan];
        }
    }
    for (int iChan = 0; iChan < NCHAN; ++iChan) {
        double outliers_expected = GAUSSIANITYEVRACTHRESHOLD * ngcounter[iChan];
        double delta = ngindex[iChan] - outliers_expected;
        double significance = sqrt (outliers_expected + ngindex[iChan] - 2 * GAUSSIANITYEVRACTHRESHOLD * ngindex[iChan]);
        ngindex[iChan] = (delta / significance);
    }
    if (verboseFLAG) std::cout << "GaussianityIndex computed" << std::endl;
    return ngindex;
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
