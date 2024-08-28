#ifndef TIN_STATS_H
#define TIN_STATS_H

// Define a structure containing TIN count, mean, median, and standard
// deviation
struct TINStats
{
    int num_transcripts;
    double mean;
    double median;
    double stddev;

    TINStats() : num_transcripts(0), mean(0), median(0), stddev(0) {}
};

#endif
