#include <iostream>
#include "silicon_calibration.hh"
#include <stdlib.h>

int main(int argc,char *argv[]){
  if(argc<2){//we can add as mutch options as we need
 std::cout<<"Please type: "<<std::endl<<
      "./event 1 \t Eta study and track reconstruction (Ester)"<<std::endl<<
      "./event 2 \t  (Francesco)"<<std::endl;
 return 0;

  }
    int who=atoi(argv[1]);
  if(who==1){
  /*  
>>>>>>> 45d8c217da05a58b9729daf7acf168b5539405e6
  analysis_ondata("/home/ester/limadou/protons_workdirectory/Servo-OFF/Run-20161112/20161112-112421-Run_3C_37MeV_SERVO_EASIROC2.root","../first_test/30MeV_ondata");
  analysis_ondata("/home/ester/limadou/protons_workdirectory/Servo-OFF/Run-20161112/RUN_3C_51MeV_SERVO_EASIROC2_HOT.root","../first_test/50MeV_ondata");
  analysis_ondata("/home/ester/limadou/protons_workdirectory/Servo-OFF/Run-20161111/RUN_POS4XC_70MeV_EASIROC2_TrigMask0_SERVO_HOT.root","../first_test/70MeV_ondata");
 analysis_ondata("/home/ester/limadou/protons_workdirectory/Servo-OFF/Run-20161112/RUN_3C_100MeV_SERVO_EASIROC2_HOT.root","../first_test/100MeV_ondata");
 analysis_ondata("/home/ester/limadou/protons_workdirectory/Servo-OFF/Run-20161112/RUN_3C_125MeV_SERVO_EASIROC2_HOT.root","../first_test/125MeV_ondata");
  analysis_ondata("/home/ester/limadou/protons_workdirectory/Servo-OFF/Run-20161112/RUN_3C_154MeV_SERVO_EASIROC2_HOT.root","../first_test/154MeV_ondata");
  analysis_ondata("/home/ester/limadou/protons_workdirectory/Servo-OFF/Run-20161112/RUN_3C_174MeV_SERVO_EASIROC2_HOT.root","../first_test/174MeV_ondata");
   analysis_ondata("/home/ester/limadou/protons_workdirectory/Servo-OFF/Run-20161112/RUN_3C_228MeV_SERVO_EASIROC2_HOT.root","../first_test/228MeV_ondata");

  analysis_ondata("/home/ester/limadou/protons_workdirectory/Servo-OFF/Run-20161112/RUN_ROTAZIONE_70MeV_SERVO_EASIROC2_HOT.root","../first_test/Rotation_ondata");
    
    analysis_ondata("/home/ester/limadou/backup_geant4_20_12_2016/muons/20161024-152620-Run_48059_4660.root ","muons");
    */
    
    //std::cout<<"It works!!!"<<std::endl;
*/
  }
  else if(who==2){


  }
  else
    std::cout<<"Wrong option. Please type:"<<std::endl<<
      "./event 1 \t Eta study and track reconstruction (Ester)"<<std::endl<<
      "./event 2 \t  (Francesco)"<<std::endl;
      
    


  return 1;



}
