#include "alignSequence.hpp"

#include <cuda.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>

void parseArguments(int argc, char *argv[])
{
    if (argc < 2)
        std::cerr << "Usage: " << SequenceAlignment::USAGE << std::endl;

    for (int i = 1; i < argc; ++i)
    {
        if (argv[i][0] == '-' && strlen(argv[i]) > 1)
        {
            if (SequenceAlignment::argumentMap.count(argv[i][1]) > 0)
            {
                auto setArg = SequenceAlignment::argumentMap.at(argv[i][1]);
                SequenceAlignment::deviceType = (setArg == SequenceAlignment::programArgs::CPU) ||
                                                (setArg == SequenceAlignment::programArgs::GPU)
                                                ? setArg
                                                : SequenceAlignment::deviceType;
                SequenceAlignment::sequenceType = (setArg == SequenceAlignment::programArgs::DNA) ||
                                                  (setArg == SequenceAlignment::programArgs::PROTEIN)
                                                  ? setArg
                                                  : SequenceAlignment::sequenceType;
            }
            else
            {
                std::cerr << "Ignoring \"" << argv[i] << "\"" << std::endl;
            }
        }
        else
        {
            std::ifstream f(argv[i]);
            if (f.good())
            {
                // :TODO: Read 2 files: text and pattern.
                //        What format are the sequences in? Just one line of characters?
            }
            else
            {
                std::cout << argv[i] << " file does not exist" << std::endl;
            }

            f.close();
        }
    }
}

int main(int argc, char *argv[])
{
    parseArguments(argc, argv);
}
