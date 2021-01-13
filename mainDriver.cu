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

    if (SequenceAlignment::deviceType == SequenceAlignment::programArgs::CPU)
        SequenceAlignment::alignSequenceCPU(request, &response);

    return 0;
}
