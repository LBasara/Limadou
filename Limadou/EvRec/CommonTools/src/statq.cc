/*
 * statq.cc
 *
 *
 *
 */


#include "statq.hh"


statq::statq (std::vector<float> v)
{
    vraw = v;
    vsize = v.size();
    ComputeVars();
}

void statq::ComputeVars()
{
    sum = std::accumulate (std::begin (vraw), std::end (vraw), 0.0);
    mean = sum / vsize;
    sq_sum = 0.0;
    std::for_each (std::begin (vraw), std::end (vraw), [&] (const double d) {
        sq_sum += (d - mean) * (d - mean);
    });
    stdev = sqrt (sq_sum / vsize );
}

