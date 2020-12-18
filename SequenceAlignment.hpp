#pragma once

#include <string>
#include <map>


namespace SequenceAlignment
{
    const std::string USAGE = "Usage: alignSequence [-p -d -cpu] [textSequenceFile] [patternSequenceFile]\n";

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

    const unsigned int MAX_SEQUENCE_LEN = 524288;

    char textBytes[MAX_SEQUENCE_LEN];
    int textNumBytes = 0;

    char patternBytes[MAX_SEQUENCE_LEN];
    int patternNumBytes = 0;

} // namespace SequenceAlignment
