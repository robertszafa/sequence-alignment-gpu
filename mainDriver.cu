#include "SequenceAlignment.hpp"

#include <cuda.h>

#include <iostream>
#include <fstream>


int main(int argc, const char *argv[])
{
    if (parseArguments(argc, argv) == -1) return -1;

    if (SequenceAlignment::deviceType == SequenceAlignment::programArgs::CPU)
        SequenceAlignment::alignSequenceCPU();

    return 0;
}
