#pragma once

#include <string>
#include <unordered_map>


namespace SequenceAlignment
{
    const std::string USAGE = "Usage: alignSequence [-p -d -c -g] [-m scoreMatrixFile] textSequenceFile patternSequenceFile \
                              \n\t-p, -d             -- protein or dna sequence (dna by default)\
                              \n\t-c, -g             -- cpu or gpu device (cpu by default)\
                              \n\t-m scoreMatrixFile -- file with custom scoring matrix (blast for dna, XYZ for protein by default)\n";
    const std::string SEQ_NOT_READ_ERROR = "Error:\ttextSequence or patternSequence not read\n";
    const std::string SCORE_MATRIX_NOT_READ_WARNING = "Warining:\ttextSequence or patternSequence not read\n";

    enum programArgs
    {
        /// Which device should be used to for the algorithm?
        CPU = 1, GPU = 2,
        /// What is the type of the input sequence?
        DNA = 4, PROTEIN = 8,
        /// The next argument is a file with the custom score matrix.
        SCORE_MATRIX = 16
    };
    const std::unordered_map<char, int> argumentMap = {
        { 'c', programArgs::CPU},
        { 'g', programArgs::GPU},
        { 'd', programArgs::DNA},
        { 'p', programArgs::PROTEIN},
        { 'm', programArgs::SCORE_MATRIX},
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

    /// The empty character is also counted.
    const unsigned int NUM_DNA_CHARS = 5;
    const unsigned int NUM_PROTEIN_CHARS = 25;

    /// Each combination of 2 sequence characters is mapped to a score value.
    /// A unique 16-bit integer key can be created from 2 characters: (char1 << 8) & char2
    std::unordered_map<short, short> scoreMap;

    /// The order of scores in a file with a scoring matrix is fixed ('*' represents a gap).
    const char dnaScoreMatrixCharOrder[] =  {'A', 'T', 'C', 'G', '*'};
    const char proteinScoreMatrixCharOrder[] =  {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                                                 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
                                                 'B', 'J', 'Z', 'X', '*'};

    const std::string defaultDnaScoreMatrixFile = "scoreMatrices/dna/blast.txt";
    const std::string defaultProteinScoreMatrixFile = "scoreMatrices/protein/blosum50.txt";

} // namespace SequenceAlignment
