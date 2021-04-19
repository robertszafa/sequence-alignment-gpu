#include "SequenceAlignment.hpp"


int main(int argc, const char *argv[])
{
    int ERR = 0;

    SequenceAlignment::Request request;
    SequenceAlignment::Response response;

    // Fill the request with user specified arguments.
    ERR = parseArguments(argc, argv, &request);

    if (ERR) return 1;

    // Runtime dispatch based on device and alignment algorithm.
    if (request.deviceType == SequenceAlignment::programArgs::CPU)
        ERR = SequenceAlignment::alignSequenceCPU(request, &response);
    else if (request.deviceType == SequenceAlignment::programArgs::GPU)
        ERR = SequenceAlignment::alignSequenceGPU(request, &response);

    if (ERR) return 1;

    prettyAlignmentPrint(response, std::cout);

    return 0;
}
