/*
BasicStatistics.cpp:
Utility functions for basic statistics

*/

#include <vector>
#include <numeric>  // std::accumulate
#include<algorithm>  // std::foreach
#include <math.h>  // sqrt, pow

#include "basic_statistics.h"

// Compute basic statistics
double computeMean(std::vector<int> values)
{
    // Mean
    int size = values.size();
    double sum = std::accumulate(std::begin(values), std::end(values), 0.0);
    double mean = sum / size;

    return mean;
}

double computeStd(std::vector<int> values)
{
    // Standard deviation
    // Σ(value - mean)²
    double mean = computeMean(values);
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
    int size = values.size();
    if (size % 2 != 0) {
        // Median is the middle value
        median = (double)values[size/2];
    } else {
        // Median is the mean of the two middle values
        median = (double)(values[(size-1)/2] + values[size/2])/2.0;
    }

    return median;
};

double computeMoment(int moment_value, double mean, std::vector<int> values)
{
    // Moment
    double accum = 0.0;
    std::for_each (std::begin(values), std::end(values), [&](const double d) {
        accum += pow((d - mean), moment_value);
    });
    double moment = accum / values.size();

    return moment;
}

double computePearsonsSkewnessCoeff(std::vector<int> values)
{
    // Pearson's Skewness Coefficient
    double mean = computeMean(values);
    double m3 = computeMoment(3, mean, values);
    double m2 = computeMoment(2, mean, values);
    double coeff = m3 / pow(sqrt(m2), 3);

    return coeff;
};

double computeKurtosis(std::vector<int> values)
{
    // Kurtosis
    double mean = computeMean(values);
    double m2 = computeMoment(2, mean, values);
    double m4 = computeMoment(4, mean, values);
    double kurtosis = m4 / pow(m2, 2.);

    return kurtosis;
};
