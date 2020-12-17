#include "SequenceAlignment.hpp"

#include <cuda.h>

#include <iostream>
#include <fstream>

void parseArguments(int argc, char *argv[])
{
    // We need at least the text and pattern file.
    if (argc < 3) std::cerr << "Usage: " << SequenceAlignment::USAGE << std::endl;

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
                // Use string's range constructor to copy over entire file to memory.
                std::string fileString((std::istreambuf_iterator<char>(f)),
                                        std::istreambuf_iterator<char>());

                if (SequenceAlignment::textNumBytes == 0)
                {
                    SequenceAlignment::textBytes = fileString.c_str();
                    SequenceAlignment::textNumBytes = fileString.length();
                }
                else
                {
                    SequenceAlignment::patternBytes = fileString.c_str();
                    SequenceAlignment::patternNumBytes = fileString.length();
                }

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

    std::cout << SequenceAlignment::textNumBytes << " " << SequenceAlignment::patternNumBytes << "\n";
}
