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
        /// The next argument is a file with the custom score matrix.
        SCORE_MATRIX,
        /// The next argument is the gap open penalty.
        GAP_OPEN,
        /// The next argument is the gap extend penalty.
        GAP_EXTEND,
    };
    const std::unordered_map<std::string, programArgs> argumentMap = {
        { "--cpu", programArgs::CPU}, { "-c", programArgs::CPU},
        { "--gpu", programArgs::GPU}, { "-g", programArgs::GPU},
        { "--dna", programArgs::DNA}, { "-d", programArgs::DNA},
        { "--protein", programArgs::PROTEIN}, { "-p", programArgs::PROTEIN},
        { "--score-matrix", programArgs::SCORE_MATRIX}, { "-s", programArgs::SCORE_MATRIX},
        { "--gap-open", programArgs::GAP_OPEN},
        { "--gap-extend", programArgs::GAP_EXTEND},
    };

    /// User messages for stdout.
    const std::string USAGE = "\
Usage: alignSequence [-p -d -c -g] [-m scoreMatrixFile] textSequenceFile patternSequenceFile\n\
       -d, --dna             - align dna sequences (default)\n\
       -p, --protein         - align protein sequence\n\
       -c, --cpu             - use cpu device (default)\n\
       -g, --gpu             - use gpu device\n\
       -s, --score-matrix    - next argument is a score matrix file\n\
       --gap-open            - next argument is a gap open penalty (default 10)\n\
       --gap-extend          - next argument is a gap extend penalty (default 1)\n";
    const std::string SEQ_NOT_READ_ERROR =
        "error: text sequence or pattern sequence not read.\n";
    const std::string SCORE_MATRIX_NOT_READ_ERROR =
        "error: matrix scores not read. Only integer scores accepted (-32,768 to 32,767)\n";
    const std::string GAP_OPEN_NOT_READ_ERROR =
        "error: gap penalty not read. Only integer scores accepted (-32,768 to 32,767)\n";
    const std::string GAP_EXTEND_NOT_READ_ERROR =
        "error: gap extend penalty not read. Only integer scores accepted (-32,768 to 32,767)\n";

    const unsigned int NUM_DNA_CHARS = 4;
    const unsigned int NUM_PROTEIN_CHARS = 24;
    /// The order of scores in a file with a scoring matrix is fixed ('*' represents a gap).
    /// Characters are transformed into ordered integers, e.g. for DNA A->0, T->1, C->3, ...
    const char DNA_ALPHABET[] =  {'A', 'T', 'C', 'G'};
    const char PROTEIN_ALPHABET[] =  {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
                                      'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X'};

    const unsigned int MAX_SEQUENCE_LEN = 4096;

    /// Default program arguments.
    const programArgs DEFAULT_DEVICE = programArgs::CPU;
    const programArgs DEFAULT_SEQUENCE = programArgs::DNA;
    static const char *DEFAULT_ALPHABET = DNA_ALPHABET;
    const int DEFAULT_ALPHABET_SIZE = NUM_DNA_CHARS;
    const short DEFAULT_GAP_OPEN_SCORE = -5;
    const short DEFAULT_GAP_EXTEND_SCORE = -1;
    const std::string DEFAULT_DNA_SCORE_MATRIX_FILE = "scoreMatrices/dna/blast.txt";
    const std::string DEFAULT_PROTEIN_SCORE_MATRIX_FILE = "scoreMatrices/protein/blosum50.txt";


    struct Request
    {
        /// Which device should be used for the request.
        programArgs deviceType;
        /// What is the type of the sequences.
        programArgs sequenceType;
        /// Buffer holding the text sequence.
        char textBytes[MAX_SEQUENCE_LEN]; int textNumBytes;
        /// Buffer holding the pattern sequence.
        char patternBytes[MAX_SEQUENCE_LEN]; int patternNumBytes;
        /// Alphabet of the sequence.
        const char *alphabet; int alphabetSize;
        /// Substitution matrix stored in row major order.
        short scoreMatrix[NUM_PROTEIN_CHARS * NUM_PROTEIN_CHARS];
        /// Penalties for opening and extending gaps.
        short gapOpenScore; short gapExtendScore;
    };

    struct Response
    {
        /// Buffer holding the aligned text sequence.
        char alignedTextBytes[MAX_SEQUENCE_LEN]; int alignedTextNumBytes;
        /// Buffer holding the aligned pattern sequence.
        char alignedPatternBytes[MAX_SEQUENCE_LEN]; int alignedPatternNumBytes;
    };


    void alignSequenceCPU(const Request&, Response*);

} // namespace SequenceAlignment


/** Implementation files */

#include "utilities.cpp"
#include "alignSequenceCPU.cpp"
