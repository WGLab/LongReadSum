#ifndef BASICSTATISTICS_MODULE_H_
#define BASICSTATISTICS_MODULE_H_


// Compute basic statistics
double computeMean(std::vector<int> values);

double computeStd(std::vector<int> values);

double computeMedian(std::vector<int> values);

double computeMoment(int moment_value, double mean, std::vector<int> values);

double computePearsonsSkewnessCoeff(std::vector<int> values);

double computeKurtosis(std::vector<int> values);

#endif
