// Utility functions
#ifndef UTILS_H
#define UTILS_H

/// @cond
#include <string>
/// @endcond

// Print a message to stdout in a thread-safe manner
void printMessage(std::string message);

// Print an error message to stderr in a thread-safe manner
void printError(std::string message);

#endif // UTILS_H
