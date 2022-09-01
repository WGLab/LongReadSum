/*
BasicStatistics.cpp:
Utility functions for basic statistics

*/

#include <vector>
#include <numeric>  // std::accumulate
#include<algorithm>  // std::foreach
#include <math.h>  // sqrt

#include "BasicStatistics.h"

// Compute basic statistics
double computeMean(std::vector<int> values)
{
    // Mean
    int size = values.size();
    double sum = std::accumulate(std::begin(values), std::end(values), 0.0);
    double mean =  sum / size;

    return mean
}

double computeStd(std::vector<int> values)
{
    // Standard deviation
    // Σ(value - mean)²
    double accum = 0.0;
    std::for_each (std::begin(values), std::end(values), [&](const double d) {
        accum += (d - mean) * (d - mean);
    });

    // sqrt( Σ(value - mean)² / N-1 )
    double std = sqrt(accum / (values.size()-1));

    return std;
}

double computeMedian(std::vector<int> values)
{
    // Median
    double median;
    std::sort(values.begin(), values.end());
    if (values.size() % 2 != 0) {
        // Median is the middle value
        median = (double)values[size/2];
    } else {
        // Median is the mean of the two middle values
        median = (double)(values[(size-1)/2] + values[size/2])/2.0;
    }

    return median
}
