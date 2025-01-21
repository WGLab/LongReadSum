#include "utils.h"

/// @cond
#include <iostream>
#include <iomanip>
#include <string>
#include <mutex>
#include <stdio.h>
#include <sys/resource.h>  // getrusage
/// @endcond


// Define print mutex
std::mutex print_mtx;

// Thread-safe print message function
void printMessage(std::string message)
{
    std::lock_guard<std::mutex> lock(print_mtx);
    std::cout << message << std::endl;
}

// Thread-safe print error function
void printError(std::string message)
{
    std::lock_guard<std::mutex> lock(print_mtx);
    std::cerr << message << std::endl;
}

void printMemoryUsage(const std::string& functionName) {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

    // Convert from KB to GB
    double mem_usage_gb = (double)usage.ru_maxrss / 1024.0 / 1024.0;
    std::lock_guard<std::mutex> lock(print_mtx);
    std::cout << functionName << " memory usage: "
              << std::fixed << std::setprecision(2) << mem_usage_gb << " GB" << std::endl;
}
