#include "LEvRec0File.hh"


LEvRec0File::LEvRec0File(const char *inpFile) {
  inputCalib=TFile::Open(inpFile, "READ");
  //calibration run processing...
  treeCalib=(TTree*)inputCalib->Get("T");
  TTree *Tmd = (TTree*)inputCalib->Get("Tmd");
  unsigned short run_id;
  Tmd->SetBranchAddress("run_id",&run_id);
  Tmd->GetEntry(0);
  RunId = static_cast<int>(run_id);
}


int LEvRec0File::SetTheEventPointer(LEvRec0 &ev) {
  treeCalib->SetBranchAddress("strip[4608]",&ev.strip);
  treeCalib->SetBranchAddress("trigger_index", &ev.trigger_index);
  treeCalib->SetBranchAddress("hepd_time", &ev.hepd_time);
  treeCalib->SetBranchAddress("event_index", &ev.event_index);
  treeCalib->SetBranchAddress("event_length", &ev.event_length);
  treeCalib->SetBranchAddress("pmt_high[64]", &ev.pmt_high);
  treeCalib->SetBranchAddress("pmt_low[64]", &ev.pmt_low);
  treeCalib->SetBranchAddress("rate_meter[9]", &ev.rate_meter);
  treeCalib->SetBranchAddress("trigger_flag[64]", &ev.trigger_flag);
  treeCalib->SetBranchAddress("alive_time", &ev.alive_time);
  treeCalib->SetBranchAddress("dead_time", &ev.dead_time);

  treeCalib->SetBranchStatus("*",kFALSE);
  treeCalib->SetBranchStatus("strip[4608]",kTRUE);
  treeCalib->SetBranchStatus("event_index",kTRUE);

  return 0;
}



int LEvRec0File::GetEntry(int iEntry) {
  treeCalib->GetEntry(iEntry);
  return 0;
}



int LEvRec0File::GetEntries() {
  return treeCalib->GetEntries();
}

void LEvRec0File::Close() {
  if(inputCalib) {
    treeCalib=0;
    inputCalib->Close();
    inputCalib=0;
  }
  return;
}



LEvRec0File::~LEvRec0File() {
  if(inputCalib) {
    treeCalib=0;
    inputCalib->Close();
  }
}
