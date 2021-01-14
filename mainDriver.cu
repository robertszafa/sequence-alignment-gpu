#include "SequenceAlignment.hpp"

#include <cuda.h>

#include <iostream>
#include <fstream>


int main(int argc, const char *argv[])
{
    SequenceAlignment::Request request = {};
    SequenceAlignment::Response response;

    // Fill the request with user specified arguments.
    if (parseArguments(argc, argv, &request) == -1) return -1;

    if (request.deviceType == SequenceAlignment::programArgs::CPU)
        SequenceAlignment::alignSequenceCPU(request, &response);

    std::cout << std::string(response.alignedTextBytes, (response.alignedTextBytes + response.numAlignmentBytes)) << "\n";
    std::cout << std::string(response.alignedPatternBytes, (response.alignedPatternBytes + response.numAlignmentBytes)) << "\n";
    std::cout << "Score: " << response.score << "\n";

    return 0;
}
