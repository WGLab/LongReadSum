#include "utils.h"

/// @cond
#include <stdio.h>
#include <string>
#include <iostream>
#include <mutex>
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
