#pragma once

#include <string>
#include <unordered_map>

namespace SequenceAlignment
{
    enum programArgs
    {
        /// Which device should be used to for the algorithm?
        CPU, GPU,
        /// What is the type of the input sequence?
        DNA, PROTEIN,
        /// Which type of alignment algorithm should be used.
        GLOBAL, LOCAL, SEMI_GLOBAL,
        /// The next argument is a file with the custom score matrix.
        SCORE_MATRIX,
        /// The next argument is the gap penalty.
        GAP_PENALTY,
    };
    const std::unordered_map<std::string, programArgs> argumentMap = {
        { "--cpu", programArgs::CPU}, { "-c", programArgs::CPU},
        { "--gpu", programArgs::GPU}, { "-g", programArgs::GPU},
        { "--dna", programArgs::DNA}, { "-d", programArgs::DNA},
        { "--protein", programArgs::PROTEIN}, { "-p", programArgs::PROTEIN},
        { "--global", programArgs::GLOBAL},
        { "--local", programArgs::LOCAL},
        { "--semi-global", programArgs::SEMI_GLOBAL},
        { "--score-matrix", programArgs::SCORE_MATRIX}, { "-s", programArgs::SCORE_MATRIX},
        { "--gap-penalty", programArgs::GAP_PENALTY},
    };

    /// User messages for stdout.
    const std::string USAGE = "\
Usage: alignSequence [-p -d -c -g] [-m scoreMatrixFile] textSequenceFile patternSequenceFile\n\
       -d, --dna             - align dna sequences (default)\n\
       -p, --protein         - align protein sequence\n\
       -c, --cpu             - use cpu device (default)\n\
       -g, --gpu             - use gpu device\n\
       --global              - use global alignment (default)\n\
       --semi-global         - use semi global alignment\n\
       --local               - use local alignment\n\
       -s, --score-matrix    - next argument is a score matrix file\n\
       --gap-penalty         - [int] next argument is a gap open penalty (default 5)\n";
    const std::string SEQ_NOT_READ_ERROR = "error: text sequence or pattern sequence not read\n";
    const std::string MEM_ERROR = "error: sequence is too long, not enough memory\n";
    const std::string TEXT_SHORTER_THAN_PATTERN_ERROR =
        "error: text sequence cannot be shorter than the pattern sequence.\n";
    const std::string SCORE_MATRIX_NOT_READ_ERROR =
        "error: matrix scores not read. Only integer scores accepted (-32,768 to 32,767)\n";
    const std::string GAP_PENALTY_NOT_READ_ERROR =
        "error: gap penalty not read. Only integer scores accepted (-32,768 to 32,767)\n";

    const unsigned int NUM_DNA_CHARS = 4;
    const unsigned int NUM_PROTEIN_CHARS = 23;
    /// The order of scores in a file with a scoring matrix is fixed ('*' represents a gap).
    /// Characters are transformed into ordered integers, e.g. for DNA A->0, T->1, C->3, ...
    const char DNA_ALPHABET[] =  {'A', 'T', 'C', 'G', '-'};
    const char PROTEIN_ALPHABET[] =  {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
                                      'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '-'};

    const unsigned int MAX_SEQUENCE_LEN = 8192;

    /// Default program arguments.
    const programArgs DEFAULT_DEVICE = programArgs::CPU;
    const programArgs DEFAULT_SEQUENCE = programArgs::DNA;
    const programArgs DEFAULT_ALIGNMENT_TYPE = programArgs::GLOBAL;
    static const char *DEFAULT_ALPHABET = DNA_ALPHABET;
    const int DEFAULT_ALPHABET_SIZE = NUM_DNA_CHARS;
    const short DEFAULT_GAP_PENALTY = 5;
    const std::string DEFAULT_DNA_SCORE_MATRIX_FILE = "scoreMatrices/dna/blast.txt";
    const std::string DEFAULT_PROTEIN_SCORE_MATRIX_FILE = "scoreMatrices/protein/blosum50.txt";


    struct Request
    {
        /// Which device should be used for the request.
        programArgs deviceType;
        /// What is the type of the sequences.
        programArgs sequenceType;
        /// What is the type of the aignment algorithm.
        programArgs alignmentType;
        /// Buffer holding the text sequence.
        char *textBytes = nullptr;
        int textNumBytes;
        /// Buffer holding the pattern sequence.
        char *patternBytes = nullptr;
        int patternNumBytes;
        /// Alphabet of the sequence.
        const char *alphabet; int alphabetSize;
        /// Substitution matrix stored in row major order.
        short scoreMatrix[NUM_PROTEIN_CHARS * NUM_PROTEIN_CHARS];
        /// Penalty for a gap.
        short gapPenalty;

        ~Request()
        {
            if (textBytes) delete[] textBytes;
            if (patternBytes) delete[] patternBytes;
            textBytes = nullptr;
            patternBytes = nullptr;
        }
    };

    struct Response
    {
        /// Buffer holding the aligned text sequence.
        char *alignedTextBytes = nullptr;
        /// Buffer holding the aligned pattern sequence.
        char *alignedPatternBytes = nullptr;
        int numAlignmentBytes;
        int score;

        ~Response()
        {
            if (alignedTextBytes) delete[] alignedTextBytes;
            if (alignedPatternBytes) delete[] alignedPatternBytes;
            alignedTextBytes = nullptr;
            alignedPatternBytes = nullptr;
        }
    };

    enum DIR { LEFT, DIAG, TOP};


    void alignSequenceGlobalCPU(const Request&, Response*);

    void alignSequenceGlobalGPU(const Request&, Response*);

    void traceBack(const char*, const unsigned int, const unsigned int, const Request&, Response*);

} // namespace SequenceAlignment


/** Implementation files */

#include "utilities.cpp"
#include "alignSequenceCPU.cpp"
#include "alignSequenceGPU.cu"
