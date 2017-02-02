#include <iostream>
#include "silicon_calibration.hh"
#include <stdlib.h>

int main(int argc,char *argv[]){
  if(argc==4) analysis(argv[1],argv[2],argv[3]);
  else if(argc==3) analysis_ondata(argv[1],argv[2]);
  else{
    std::cout<<"Please type: "<<std::endl<<
      "./event data_file_path calibration_file output_name"<<std::endl;
      return 0; 
  }

  return 1;



}
