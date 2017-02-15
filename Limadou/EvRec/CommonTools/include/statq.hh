/*
 * statq.hh
 *
 *
 *
 */


#ifndef STATQ_HH
#define STATQ_HH

#include <math.h>
#include <vector>
#include <numeric>
#include <algorithm>


class statq
{
   public:
      statq(std::vector<float> v);
      float GetSum   (){return sum;}
      float GetMean  (){return mean;}
      float GetSqSum (){return sq_sum;}
      float GetStdDev(){return stdev;}
      float GetSize() {return vsize;}


   private:
      std::vector<float> vraw;
      unsigned int vsize;
      float sum;
      float mean;
      float sq_sum;
      float stdev;
      void ComputeVars();

};

#endif /* STATQ_HH */
