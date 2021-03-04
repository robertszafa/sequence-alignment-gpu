#include "SequenceAlignment.hpp"

#include <iostream>

int main(int argc, const char *argv[])
{
    SequenceAlignment::Request request;
    SequenceAlignment::Response response;

    // Fill the request with user specified arguments.
    if (parseArguments(argc, argv, &request) == -1) return -1;

    // Runtime dispatch based on device and alignment algorithm.
    if (request.deviceType == SequenceAlignment::programArgs::CPU)
        SequenceAlignment::alignSequenceCPU(request, &response);
    else if (request.deviceType == SequenceAlignment::programArgs::GPU)
        SequenceAlignment::alignSequenceGPU(request, &response);

    prettyAlignmentPrint(response, std::cout);

    return 0;
}
