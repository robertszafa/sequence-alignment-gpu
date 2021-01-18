#include "SequenceAlignment.hpp"

#include <cuda.h>

#include <iostream>
#include <fstream>


int main(int argc, const char *argv[])
{
    SequenceAlignment::Request request;
    SequenceAlignment::Response response;

    // Fill the request with user specified arguments.
    if (parseArguments(argc, argv, &request) == -1) return -1;

    // Runtime dispatch based on device and alignment algorithm.
    if (request.deviceType == SequenceAlignment::programArgs::CPU &&
        request.alignmentType == SequenceAlignment::programArgs::GLOBAL)
        SequenceAlignment::alignSequenceGlobalCPU(request, &response);

    prettyAlignmentPrint(response.alignedTextBytes, response.alignedPatternBytes,
                         response.numAlignmentBytes, std::cout);
    std::cout << "Score: " << response.score << "\n";

    return 0;
}
