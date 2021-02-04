#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"


#include "../SequenceAlignment.hpp"


TEST_CASE("indexOfLetter")
{
    REQUIRE(indexOfLetter('A', SequenceAlignment::DNA_ALPHABET, SequenceAlignment::NUM_DNA_CHARS) == 0);
    REQUIRE(indexOfLetter('H', SequenceAlignment::DNA_ALPHABET, SequenceAlignment::NUM_DNA_CHARS) == -1);
    REQUIRE(indexOfLetter('H', SequenceAlignment::PROTEIN_ALPHABET, SequenceAlignment::NUM_PROTEIN_CHARS) == 8);
}

TEST_CASE("parseScoreMatrixFile")
{
    SequenceAlignment::Request request;
    request.alphabet = SequenceAlignment::DNA_ALPHABET;
    request.alphabetSize = SequenceAlignment::NUM_DNA_CHARS;
    parseScoreMatrixFile("scoreMatrices/dna/blast.txt", request.alphabetSize, request.scoreMatrix);

    REQUIRE(getScore('A', 'A', request.alphabet, request.alphabetSize, request.scoreMatrix) == 5);
    REQUIRE(getScore('G', 'T', request.alphabet, request.alphabetSize, request.scoreMatrix) == -4);
}

TEST_CASE("readSequenceBytes")
{
    const int argc = 3;
    const char *argv[argc] = { "./alignSequence", "data/dna/dna_01.txt", "data/dna/dna_02.txt"};
    SequenceAlignment::Request request;
    parseArguments(argc, argv, &request);

    const char expectedText[] = {0, 2, 0, 2};
    const char expectedPattern[] = {2, 2, 1, 0};

    REQUIRE(std::equal(request.textBytes, request.textBytes + request.textNumBytes, expectedText));
    REQUIRE(std::equal(request.patternBytes, request.patternBytes + request.patternNumBytes, expectedPattern));
}

TEST_CASE("parseArguments")
{
    // Catch stderr to string to test error messages.
    std::stringstream buffer;
    std::streambuf *old = std::cerr.rdbuf(buffer.rdbuf());

    SECTION("usage")
    {
        const int argc = 1;
        const char *argv[argc] = { "./alignSequence"};
        SequenceAlignment::Request request;
        parseArguments(argc, argv, &request);

        std::string stderrString = buffer.str();
        REQUIRE(stderrString == SequenceAlignment::USAGE);
    }

    SECTION("no or empty sequence files")
    {
        const int argc = 3;
        const char *argv[argc] = { "./alignSequence", "-p", "-c"};
        SequenceAlignment::Request request;
        parseArguments(argc, argv, &request);

        std::string expectedMsg = SequenceAlignment::SEQ_NOT_READ_ERROR + SequenceAlignment::USAGE;
        std::string stderrString = buffer.str();
        REQUIRE(stderrString == expectedMsg);

        REQUIRE(request.deviceType == SequenceAlignment::programArgs::CPU);
        REQUIRE(request.sequenceType == SequenceAlignment::programArgs::PROTEIN);
    }

    SECTION("incorrect score matrix")
    {
        const int argc = 5;
        const char *argv[argc] = { "./alignSequence", "--score-matrix", "test/corruptScoreMatrix.txt",
                                   "data/dna/dna_01.txt", "data/dna/dna_02.txt"};
        SequenceAlignment::Request request;
        parseArguments(argc, argv, &request);

        std::string expectedMsg = SequenceAlignment::SCORE_MATRIX_NOT_READ_ERROR;
        std::string stderrString = buffer.str();
        REQUIRE(stderrString == expectedMsg);
    }

    // Restore old cerr.
    std::cerr.rdbuf(old);
}

TEST_CASE("alignSequenceGlobalCPU")
{

    SECTION("DNA_01")
    {
        const int argc = 6;
        const char *argv[argc] = { "./alignSequence",  "--gap-penalty", "5", "--global",
                                "data/dna/dna_01.txt", "data/dna/dna_02.txt"};
        SequenceAlignment::Request request = {};
        SequenceAlignment::Response response;
        parseArguments(argc, argv, &request);

        SequenceAlignment::alignSequenceGlobalCPU(request, &response);

        // Don't check the aligned sequences since there can be multiple alignments with the same score.
        const int expectedScore = -4;
        REQUIRE(expectedScore == response.score);
    }

    SECTION("DNA_02")
    {
        SequenceAlignment::Request request;
        SequenceAlignment::Response response;

        std::string text = "GCCT";
        std::string pattern = "GGTC";
        request.deviceType = SequenceAlignment::programArgs::CPU;
        request.sequenceType = SequenceAlignment::programArgs::DNA;
        request.alignmentType = SequenceAlignment::programArgs::GLOBAL;
        request.alphabet = SequenceAlignment::DNA_ALPHABET;
        request.alphabetSize = SequenceAlignment::NUM_DNA_CHARS;
        request.gapPenalty = 5;
        request.textNumBytes = validateAndTransform(text, request.alphabet, request.alphabetSize);
        request.textBytes = new char[request.textNumBytes];
        std::copy_n(text.begin(), request.textNumBytes, request.textBytes);
        request.patternNumBytes = validateAndTransform(pattern, request.alphabet, request.alphabetSize);
        request.patternBytes = new char[request.patternNumBytes];
        std::copy_n(pattern.begin(), request.patternNumBytes, request.patternBytes);
        parseScoreMatrixFile(SequenceAlignment::DEFAULT_DNA_SCORE_MATRIX_FILE, request.alphabetSize, request.scoreMatrix);

        SequenceAlignment::alignSequenceGlobalCPU(request, &response);

        // Don't check the aligned sequences since there can be multiple alignments with the same score.
        const int expectedScore = -4;
        REQUIRE(expectedScore == response.score);
    }

    SECTION("DNA_03")
    {
        SequenceAlignment::Request request;
        SequenceAlignment::Response response;

        std::string text = "TTCGCCT";
        std::string pattern = "CTCGGTC";
        request.deviceType = SequenceAlignment::programArgs::CPU;
        request.sequenceType = SequenceAlignment::programArgs::DNA;
        request.alignmentType = SequenceAlignment::programArgs::GLOBAL;
        request.alphabet = SequenceAlignment::DNA_ALPHABET;
        request.alphabetSize = SequenceAlignment::NUM_DNA_CHARS;
        request.gapPenalty = 5;
        request.textNumBytes = validateAndTransform(text, request.alphabet, request.alphabetSize);
        request.textBytes = new char[request.textNumBytes];
        std::copy_n(text.begin(), request.textNumBytes, request.textBytes);
        request.patternNumBytes = validateAndTransform(pattern, request.alphabet, request.alphabetSize);
        request.patternBytes = new char[request.patternNumBytes];
        std::copy_n(pattern.begin(), request.patternNumBytes, request.patternBytes);
        parseScoreMatrixFile(SequenceAlignment::DEFAULT_DNA_SCORE_MATRIX_FILE, request.alphabetSize, request.scoreMatrix);

        SequenceAlignment::alignSequenceGlobalCPU(request, &response);

        // Don't check the aligned sequences since there can be multiple alignments with the same score.
        const int expectedScore = 2;
        REQUIRE(expectedScore == response.score);
    }

    SECTION("DNA_04")
    {
        SequenceAlignment::Request request;
        SequenceAlignment::Response response;

        std::string text =
            "CATAAAACTCTCGGTCGGGCTTAGTACCAGGACCGGCGCACCAGAGTGTCAATCACGACCCTTCACACTTTGTGC";
        std::string pattern =
            "ATGAAGTTGTTCGCCTTACTTTTAATTCTACTCTCTCCTCGAGATTCGTCCGCTGAAAAATCTCTCAGCG";
        // There may be multiple alignments with the same score. These expected alignments
        // are for regression testing.
        const std::string expectedAlignedText =
            "CATAAAACTCTCGGTCGGGCTTAGTACCAGGAC--CGGCGCACCA-GAG-TGTCAATCACGACCCTTCACACTTTGT--GC-";
        const std::string expectedAlignedPattern =
            "-ATGAAG-T-T-GTTCGC-CTTACTTTTAATTCTACT-CTCTCCTCGAGAT-TCG-TC-CG-C--TGAAAAATCTCTCAGCG";
        const int expectedScore = 22;

        request.deviceType = SequenceAlignment::programArgs::CPU;
        request.sequenceType = SequenceAlignment::programArgs::DNA;
        request.alignmentType = SequenceAlignment::programArgs::GLOBAL;
        request.alphabet = SequenceAlignment::DNA_ALPHABET;
        request.alphabetSize = SequenceAlignment::NUM_DNA_CHARS;
        request.gapPenalty = 5;
        request.textNumBytes = validateAndTransform(text, request.alphabet, request.alphabetSize);
        request.textBytes = new char[request.textNumBytes];
        std::copy_n(text.begin(), request.textNumBytes, request.textBytes);
        request.patternNumBytes = validateAndTransform(pattern, request.alphabet, request.alphabetSize);
        request.patternBytes = new char[request.patternNumBytes];
        std::copy_n(pattern.begin(), request.patternNumBytes, request.patternBytes);
        parseScoreMatrixFile(SequenceAlignment::DEFAULT_DNA_SCORE_MATRIX_FILE, request.alphabetSize, request.scoreMatrix);

        SequenceAlignment::alignSequenceGlobalCPU(request, &response);

        auto gotAlignedText = std::string(response.alignedTextBytes,
                                          (response.alignedTextBytes + response.numAlignmentBytes));
        auto gotAlignedPattern = std::string(response.alignedPatternBytes,
                                             (response.alignedPatternBytes + response.numAlignmentBytes));

        REQUIRE(expectedScore == response.score);
        REQUIRE(expectedAlignedText == gotAlignedText);
        REQUIRE(expectedAlignedPattern == gotAlignedPattern);
    }

    SECTION("DNA_05")
    {
        const int argc = 6;
        const char *argv[argc] = {"./alignSequence", "--gap-penalty", "5", "--global",
                                  "data/dna/NC_018874.txt", "data/dna/mutated_NC_018874.txt"};
        SequenceAlignment::Request request = {};
        SequenceAlignment::Response response;
        parseArguments(argc, argv, &request);

        SequenceAlignment::alignSequenceGlobalCPU(request, &response);

        // Don't check the aligned sequences since there can be multiple alignments with the same score.
        const int expectedScore = 2673;
        REQUIRE(expectedScore == response.score);
    }


    SECTION("PROTEIN_01")
    {
        SequenceAlignment::Request request;
        SequenceAlignment::Response response;

        std::string text =
            "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR";
        std::string pattern =
            "MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHGKKVADALASAAGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLASHHPADFTPAVHASLDKFLASVSTVLTSKYR";
        // There may be multiple alignments with the same score. These expected alignments
        // are for regression testing.
        const std::string expectedAlignedText =
            "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR";
        const std::string expectedAlignedPattern =
            "MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHGKKVADALASAAGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLASHHPADFTPAVHASLDKFLASVSTVLTSKYR";
        const int expectedScore = 821;

        request.deviceType = SequenceAlignment::programArgs::CPU;
        request.sequenceType = SequenceAlignment::programArgs::PROTEIN;
        request.alignmentType = SequenceAlignment::programArgs::GLOBAL;
        request.alphabet = SequenceAlignment::PROTEIN_ALPHABET;
        request.alphabetSize = SequenceAlignment::NUM_PROTEIN_CHARS;
        request.gapPenalty = 5;
        request.textNumBytes = validateAndTransform(text, request.alphabet, request.alphabetSize);
        request.textBytes = new char[request.textNumBytes];
        std::copy_n(text.begin(), request.textNumBytes, request.textBytes);
        request.patternNumBytes = validateAndTransform(pattern, request.alphabet, request.alphabetSize);
        request.patternBytes = new char[request.patternNumBytes];
        std::copy_n(pattern.begin(), request.patternNumBytes, request.patternBytes);
        parseScoreMatrixFile(SequenceAlignment::DEFAULT_PROTEIN_SCORE_MATRIX_FILE, request.alphabetSize, request.scoreMatrix);

        SequenceAlignment::alignSequenceGlobalCPU(request, &response);

        auto gotAlignedText = std::string(response.alignedTextBytes,
                                          (response.alignedTextBytes + response.numAlignmentBytes));
        auto gotAlignedPattern = std::string(response.alignedPatternBytes,
                                             (response.alignedPatternBytes + response.numAlignmentBytes));

        REQUIRE(expectedScore == response.score);
        REQUIRE(expectedAlignedText == gotAlignedText);
        REQUIRE(expectedAlignedPattern == gotAlignedPattern);
    }

    SECTION("PROTEIN_02")
    {
        const int argc = 7;
        const char *argv[argc] = {"./alignSequence", "--protein", "--gap-penalty", "5", "--global",
                                  "data/protein/P0C6B8.txt", "data/protein/mutated_P0C6B8.txt"};
        SequenceAlignment::Request request = {};
        SequenceAlignment::Response response;
        parseArguments(argc, argv, &request);

        SequenceAlignment::alignSequenceGlobalCPU(request, &response);

        // Don't check the aligned sequences since there can be multiple alignments with the same score.
        const int expectedScore = 32095;
        REQUIRE(expectedScore == response.score);
    }

    SECTION("PROTEIN_03")
    {
        const int argc = 8;
        const char *argv[argc] = {"./alignSequence", "--protein", "--cpu", "--gap-penalty", "5", "--global",
                                  "data/protein/p0.txt", "data/protein/mutated_P0.txt"};
        SequenceAlignment::Request request = {};
        SequenceAlignment::Response response;
        parseArguments(argc, argv, &request);

        SequenceAlignment::alignSequenceGlobalCPU(request, &response);

        // Don't check the aligned sequences since there can be multiple alignments with the same score.
        const int expectedScore = 1328;
        REQUIRE(expectedScore == response.score);
    }

}


TEST_CASE("alignSequenceGlobalGPU")
{
    SECTION("PROTEIN_01")
    {
        const int argc = 8;
        const char *argv[argc] = {"./alignSequence", "--protein", "--gpu", "--gap-penalty", "5", "--global",
                                  "data/protein/p0.txt", "data/protein/mutated_P0.txt"};
        SequenceAlignment::Request request = {};
        parseArguments(argc, argv, &request);

        SequenceAlignment::Response responseGPU;
        SequenceAlignment::Response responseCPU;
        SequenceAlignment::alignSequenceGlobalGPU(request, &responseGPU);
        SequenceAlignment::alignSequenceGlobalCPU(request, &responseCPU);

        REQUIRE(responseCPU.score == responseGPU.score);
    }

    SECTION("PROTEIN_02")
    {
        const int argc = 8;
        const char *argv[argc] = {"./alignSequence", "--protein", "--gpu", "--gap-penalty", "5", "--global",
                                  "data/protein/P0C6B8.txt", "data/protein/mutated_P0C6B8.txt"};
        SequenceAlignment::Request request = {};
        SequenceAlignment::Response response;
        parseArguments(argc, argv, &request);

        SequenceAlignment::alignSequenceGlobalGPU(request, &response);

        // Don't check the aligned sequences since there can be multiple alignments with the same score.
        const int expectedScore = 32095;
        REQUIRE(expectedScore == response.score);
    }

}
