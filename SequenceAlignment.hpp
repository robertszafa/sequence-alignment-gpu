#pragma once

#include <string>
#include <map>


namespace SequenceAlignment
{
    const std::string USAGE = "alignSequence [-p -d -cpu] [textSequenceFile] [patternSequenceFile]";

    enum programArgs
    {
        /// Which device should be used to for the algorithm?
        CPU = 1, GPU = 2,
        /// What is the type of the input sequence?
        DNA = 4, PROTEIN = 8
    };
    const std::map<char, int> argumentMap = {
        { 'c', programArgs::CPU},
        { 'g', programArgs::GPU},
        { 'd', programArgs::DNA},
        { 'p', programArgs::PROTEIN},
    };
    /// Default is CPU.
    int deviceType = programArgs::CPU;
    /// Default is DNA.
    int sequenceType = programArgs::DNA;


    const char* textBytes;
    int textNumBytes = 0;

    const char* patternBytes;
    int patternNumBytes = 0;

} // namespace SequenceAlignment


void parseArguments(int, const char**);
