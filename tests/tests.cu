#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"

// For filesystem.
#include <stdio.h>
#include <string.h>
#include <dirent.h>

#include "../SequenceAlignment.hpp"


/// Given a directory name, return the full path of all files in the directory.
std::vector<std::string> getFilenamesInDirectory (const std::string dir)
{
    std::vector<std::string> allFiles = {};

    struct dirent *entry = nullptr;
    DIR *dp = opendir(&dir[0]);

    if (dp != nullptr)
    {
        while ((entry = readdir(dp)))
        {
            auto f = std::string(entry->d_name);
            if (f != "." && f != "..")
                allFiles.push_back(dir + "/" + f);
        }
    }

    closedir(dp);
    return allFiles;
}

TEST_CASE("indexOfLetter")
{
    CHECK(indexOfLetter('A', SequenceAlignment::DNA_ALPHABET, SequenceAlignment::NUM_DNA_CHARS) == 0);
    CHECK(indexOfLetter('H', SequenceAlignment::DNA_ALPHABET, SequenceAlignment::NUM_DNA_CHARS) == -1);
    CHECK(indexOfLetter('H', SequenceAlignment::PROTEIN_ALPHABET, SequenceAlignment::NUM_PROTEIN_CHARS) == 8);
}

TEST_CASE("parseScoreMatrixFile")
{
    SequenceAlignment::Request request;
    request.alphabet = SequenceAlignment::DNA_ALPHABET;
    request.alphabetSize = SequenceAlignment::NUM_DNA_CHARS;
    parseScoreMatrixFile("scoreMatrices/dna/blast.txt", request.alphabetSize, request.scoreMatrix);

    CHECK(getScore('A', 'A', request.alphabet, request.alphabetSize, request.scoreMatrix) == 5);
    CHECK(getScore('G', 'T', request.alphabet, request.alphabetSize, request.scoreMatrix) == -4);
}

TEST_CASE("readSequenceBytes")
{
    const int argc = 3;
    const char *argv[argc] = { "./alignSequence", "data/dna/dna_01.txt", "data/dna/dna_02.txt"};
    SequenceAlignment::Request request;
    parseArguments(argc, argv, &request);

    const char expectedText[] = {0, 2, 0, 2};
    const char expectedPattern[] = {2, 2, 1, 0};

    CHECK(std::equal(request.textBytes, request.textBytes + request.textNumBytes, expectedText));
    CHECK(std::equal(request.patternBytes, request.patternBytes + request.patternNumBytes, expectedPattern));
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
        CHECK(stderrString == SequenceAlignment::USAGE);
    }

    SECTION("no or empty sequence files")
    {
        const int argc = 3;
        const char *argv[argc] = { "./alignSequence", "-p", "-c"};
        SequenceAlignment::Request request;
        parseArguments(argc, argv, &request);

        std::string expectedMsg = SequenceAlignment::SEQ_NOT_READ_ERROR + SequenceAlignment::USAGE;
        std::string stderrString = buffer.str();
        CHECK(stderrString == expectedMsg);

        CHECK(request.deviceType == SequenceAlignment::programArgs::CPU);
        CHECK(request.sequenceType == SequenceAlignment::programArgs::PROTEIN);
    }

    SECTION("incorrect score matrix")
    {
        const int argc = 5;
        const char *argv[argc] = { "./alignSequence", "--score-matrix", "tests/corruptScoreMatrix.txt",
                                   "data/dna/dna_01.txt", "data/dna/dna_02.txt"};
        SequenceAlignment::Request request;
        parseArguments(argc, argv, &request);

        std::string expectedMsg = SequenceAlignment::SCORE_MATRIX_NOT_READ_ERROR;
        std::string stderrString = buffer.str();
        CHECK(stderrString == expectedMsg);
    }

    // Restore old cerr.
    std::cerr.rdbuf(old);
}

TEST_CASE("alignSequenceCPU - Global")
{

    SECTION("DNA_01")
    {
        const int argc = 6;
        const char *argv[argc] = { "./alignSequence",  "--gap-penalty", "5", "--global",
                                "data/dna/dna_01.txt", "data/dna/dna_02.txt"};
        SequenceAlignment::Request request;
        SequenceAlignment::Response response;
        parseArguments(argc, argv, &request);

        SequenceAlignment::alignSequenceCPU(request, &response);

        // Don't check the aligned sequences since there can be multiple alignments with the same score.
        const int expectedScore = -4;
        CHECK(expectedScore == response.score);
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

        SequenceAlignment::alignSequenceCPU(request, &response);

        // Don't check the aligned sequences since there can be multiple alignments with the same score.
        const int expectedScore = -4;
        CHECK(expectedScore == response.score);
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

        SequenceAlignment::alignSequenceCPU(request, &response);

        // Don't check the aligned sequences since there can be multiple alignments with the same score.
        const int expectedScore = 2;
        CHECK(expectedScore == response.score);
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

        SequenceAlignment::alignSequenceCPU(request, &response);

        auto gotAlignedText = std::string(response.alignedTextBytes,
                                          (response.alignedTextBytes + response.numAlignmentBytes));
        auto gotAlignedPattern = std::string(response.alignedPatternBytes,
                                             (response.alignedPatternBytes + response.numAlignmentBytes));

        CHECK(expectedScore == response.score);
        CHECK(expectedAlignedText == gotAlignedText);
        CHECK(expectedAlignedPattern == gotAlignedPattern);
    }

    SECTION("DNA_05")
    {
        const int argc = 6;
        const char *argv[argc] = {"./alignSequence", "--gap-penalty", "5", "--global",
                                  "data/dna/NC_018874.txt", "data/dna/GCA_003231495.txt"};
        SequenceAlignment::Request request;
        SequenceAlignment::Response response;
        parseArguments(argc, argv, &request);

        SequenceAlignment::alignSequenceCPU(request, &response);

        // Don't check the aligned sequences since there can be multiple alignments with the same score.
        const int expectedScore = -5991;
        CHECK(expectedScore == response.score);
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

        SequenceAlignment::alignSequenceCPU(request, &response);

        auto gotAlignedText = std::string(response.alignedTextBytes,
                                          (response.alignedTextBytes + response.numAlignmentBytes));
        auto gotAlignedPattern = std::string(response.alignedPatternBytes,
                                             (response.alignedPatternBytes + response.numAlignmentBytes));

        CHECK(expectedScore == response.score);
        CHECK(expectedAlignedText == gotAlignedText);
        CHECK(expectedAlignedPattern == gotAlignedPattern);
    }

    SECTION("PROTEIN_02")
    {
        const int argc = 7;
        const char *argv[argc] = {"./alignSequence", "--protein", "--gap-penalty", "5", "--global",
                                  "data/protein/P02232.fasta", "data/protein/P03989.fasta"};
        SequenceAlignment::Request request;
        SequenceAlignment::Response response;
        parseArguments(argc, argv, &request);

        SequenceAlignment::alignSequenceCPU(request, &response);

        // Don't check the aligned sequences since there can be multiple alignments with the same score.
        const int expectedScore = -597;
        CHECK(expectedScore == response.score);
    }

    SECTION("PROTEIN_03")
    {
        const int argc = 8;
        const char *argv[argc] = {"./alignSequence", "--protein", "--cpu", "--gap-penalty", "5", "--global",
                                  "data/protein/P05013.fasta", "data/protein/P07327.fasta"};
        SequenceAlignment::Request request;
        SequenceAlignment::Response response;
        parseArguments(argc, argv, &request);

        SequenceAlignment::alignSequenceCPU(request, &response);

        // Don't check the aligned sequences since there can be multiple alignments with the same score.
        const int expectedScore = -423;
        CHECK(expectedScore == response.score);
    }

}

TEST_CASE("alignSequenceCPU - Local")
{
    SECTION("DNA_01")
    {
        const int argc = 6;
        const char *argv[argc] = { "./alignSequence", "--gap-penalty", "5", "--local",
                                "data/dna/GCA_003231495.txt", "data/dna/dna_01.txt"};
        SequenceAlignment::Request request;
        SequenceAlignment::Response response;
        parseArguments(argc, argv, &request);

        SequenceAlignment::alignSequenceCPU(request, &response);

        auto gotAlignedText = std::string(response.alignedTextBytes,
                                          (response.alignedTextBytes + response.numAlignmentBytes));
        auto gotAlignedPattern = std::string(response.alignedPatternBytes,
                                             (response.alignedPatternBytes + response.numAlignmentBytes));
        CHECK(20 == response.score);
        CHECK("ACAC" == gotAlignedText);
        CHECK("ACAC" == gotAlignedPattern);
        CHECK(248 == response.startInAlignedText);
        CHECK(0 == response.startInAlignedPattern);
    }

    SECTION("PROTEIN_01")
    {
        const int argc = 7;
        const char *argv[argc] = {"./alignSequence", "--protein", "--gap-penalty", "10", "--local",
                                  "data/protein/P08519.fasta", "data/protein/P10635.fasta"};
        SequenceAlignment::Request request;
        SequenceAlignment::Response response;
        parseArguments(argc, argv, &request);

        SequenceAlignment::alignSequenceCPU(request, &response);

        CHECK(57 == response.score);
        CHECK(4204 == response.startInAlignedText);
        CHECK(95 == response.startInAlignedPattern);
    }

}

TEST_CASE("alignSequenceGPU - Global")
{
    SECTION("PROTEIN_01")
    {
        const int argc = 8;
        const char *argv[argc] = {"./alignSequence", "--protein", "--gpu", "--gap-penalty", "11", "--global",
                                  "data/protein/P10635.fasta", "data/protein/P02232.fasta"};
        SequenceAlignment::Request request;
        parseArguments(argc, argv, &request);

        SequenceAlignment::Response responseGPU;
        SequenceAlignment::Response responseCPU;
        SequenceAlignment::alignSequenceGPU(request, &responseGPU);
        SequenceAlignment::alignSequenceCPU(request, &responseCPU);

        CHECK(responseCPU.score == responseGPU.score);
        CHECK(std::string(responseCPU.alignedTextBytes, (responseCPU.alignedTextBytes + responseCPU.numAlignmentBytes)) ==
                std::string(responseGPU.alignedTextBytes, (responseGPU.alignedTextBytes + responseGPU.numAlignmentBytes)));
        CHECK(std::string(responseCPU.alignedPatternBytes, (responseCPU.alignedPatternBytes + responseCPU.numAlignmentBytes)) ==
                std::string(responseGPU.alignedPatternBytes, (responseGPU.alignedPatternBytes + responseGPU.numAlignmentBytes)));
    }

    SECTION("PROTEIN_02")
    {
        const int argc = 8;
        const char *argv[argc] = {"./alignSequence", "--protein", "--gpu", "--gap-penalty", "5", "--global",
                                  "data/protein/P27895.fasta", "data/protein/P27895.fasta"};
        SequenceAlignment::Request request = {};
        parseArguments(argc, argv, &request);

        SequenceAlignment::Response responseGPU;
        SequenceAlignment::Response responseCPU;
        SequenceAlignment::alignSequenceGPU(request, &responseGPU);
        SequenceAlignment::alignSequenceCPU(request, &responseCPU);

        CHECK(responseCPU.score == responseGPU.score);
        CHECK(std::string(responseCPU.alignedTextBytes, (responseCPU.alignedTextBytes + responseCPU.numAlignmentBytes)) ==
                std::string(responseGPU.alignedTextBytes, (responseGPU.alignedTextBytes + responseGPU.numAlignmentBytes)));
        CHECK(std::string(responseCPU.alignedPatternBytes, (responseCPU.alignedPatternBytes + responseCPU.numAlignmentBytes)) ==
                std::string(responseGPU.alignedPatternBytes, (responseGPU.alignedPatternBytes + responseGPU.numAlignmentBytes)));
    }

}


TEST_CASE("alignSequenceGPU - Local")
{
    SECTION("DNA_01")
    {
        const int argc = 6;
        const char *argv[argc] = {"./alignSequence", "--gap-penalty", "5", "--local",
                                  "data/dna/GCA_003231495.txt", "data/dna/dna_01.txt"};
        SequenceAlignment::Request request;
        SequenceAlignment::Response responseCPU;
        SequenceAlignment::Response responseGPU;
        parseArguments(argc, argv, &request);

        SequenceAlignment::alignSequenceCPU(request, &responseCPU);
        SequenceAlignment::alignSequenceGPU(request, &responseGPU);

        CHECK(responseGPU.score == responseCPU.score);
        CHECK(std::string(responseCPU.alignedTextBytes, (responseCPU.alignedTextBytes + responseCPU.numAlignmentBytes)) ==
                std::string(responseGPU.alignedTextBytes, (responseGPU.alignedTextBytes + responseGPU.numAlignmentBytes)));
        CHECK(std::string(responseCPU.alignedPatternBytes, (responseCPU.alignedPatternBytes + responseCPU.numAlignmentBytes)) ==
                std::string(responseGPU.alignedPatternBytes, (responseGPU.alignedPatternBytes + responseGPU.numAlignmentBytes)));
        CHECK(responseCPU.startInAlignedText == responseGPU.startInAlignedText);
        CHECK(responseCPU.startInAlignedPattern == responseGPU.startInAlignedPattern);
    }

    SECTION("PROTEIN_01")
    {
        const int argc = 7;
        const char *argv[argc] = {"./alignSequence", "--protein", "--gap-penalty", "5", "--local",
                                  "data/protein/P33450.fasta", "data/protein/P07327.fasta"};
        SequenceAlignment::Request request;
        SequenceAlignment::Response responseCPU;
        SequenceAlignment::Response responseGPU;
        parseArguments(argc, argv, &request);

        SequenceAlignment::alignSequenceCPU(request, &responseCPU);
        SequenceAlignment::alignSequenceGPU(request, &responseGPU);

        CHECK(responseGPU.score == responseCPU.score);
        CHECK(std::string(responseCPU.alignedTextBytes, (responseCPU.alignedTextBytes + responseCPU.numAlignmentBytes)) ==
                std::string(responseGPU.alignedTextBytes, (responseGPU.alignedTextBytes + responseGPU.numAlignmentBytes)));
        CHECK(std::string(responseCPU.alignedPatternBytes, (responseCPU.alignedPatternBytes + responseCPU.numAlignmentBytes)) ==
                std::string(responseGPU.alignedPatternBytes, (responseGPU.alignedPatternBytes + responseGPU.numAlignmentBytes)));
        CHECK(responseCPU.startInAlignedText == responseGPU.startInAlignedText);
        CHECK(responseCPU.startInAlignedPattern == responseGPU.startInAlignedPattern);
    }
}

TEST_CASE("Batch DNA alignment")
{
    std::vector<std::string> allDnaFiles = getFilenamesInDirectory("data/dna");

    std::cout << "\nTesting all possible combinations of 2 sequences in data/dna ...\n";

    for (auto &alignType : {"--global", "--local"})
    {
        std::cout << "Alignment type: " << alignType << "\n";
        for (int i=0; i < allDnaFiles.size(); ++i)
        {
            for (int j=i+1; j < allDnaFiles.size(); ++j)
            {
                const int argc = 7;
                const char *argv[argc] = {"./alignSequence", "--dna", "--gap-penalty", "11",
                                          &alignType[0], &allDnaFiles[i][0], &allDnaFiles[j][0]};

                SequenceAlignment::Request request;
                SequenceAlignment::Response responseCPU;
                SequenceAlignment::Response responseGPU;
                parseArguments(argc, argv, &request);

                // Avoid spending time on very long sequences.
                if (request.textNumBytes > 10000)
                    continue;

                SequenceAlignment::alignSequenceCPU(request, &responseCPU);
                SequenceAlignment::alignSequenceGPU(request, &responseGPU);
                CHECK(responseGPU.score == responseCPU.score);
                // Don't check text since there can be multiple optimal alignments.
            }
        }
    }
}

TEST_CASE("Batch Protein alignment")
{
    std::vector<std::string> allProteinFiles = getFilenamesInDirectory("data/protein");

    std::cout << "\nTesting all possible combinations of 2 sequences in data/protein ...\n";

    for (auto &alignType : {"--global", "--local"})
    {
        std::cout << "Alignment type: " << alignType << "\n";
        for (int i=0; i < allProteinFiles.size(); ++i)
        {
            for (int j=i+1; j < allProteinFiles.size(); ++j)
            {
                const int argc = 7;
                const char *argv[argc] = {"./alignSequence", "--protein", "--gap-penalty", "5",
                                          &alignType[0], &allProteinFiles[i][0], &allProteinFiles[j][0]};

                SequenceAlignment::Request request;
                SequenceAlignment::Response responseCPU;
                SequenceAlignment::Response responseGPU;
                parseArguments(argc, argv, &request);

                // Avoid spending time on very long sequences.
                if (request.textNumBytes > 10000)
                    continue;

                SequenceAlignment::alignSequenceCPU(request, &responseCPU);
                SequenceAlignment::alignSequenceGPU(request, &responseGPU);
                CHECK(responseGPU.score == responseCPU.score);
                // Don't check text since there can be multiple optimal alignments.
            }
        }
    }
}


TEST_CASE("Very long (>200k) DNA alignment")
{
    // // Comment out to no spend too much time for tests.
    // const int argc = 8;
    // const char *argv[argc] = {"./alignSequence", "--dna", "--gpu", "--gap-penalty", "5", "--global",
    //                           "data/dna/AbHV_ORF111.txt", "data/dna/mutated_AbHV_ORF111.txt"};
    // SequenceAlignment::Request request = {};
    // parseArguments(argc, argv, &request);

    // SequenceAlignment::Response responseGPU;
    // SequenceAlignment::Response responseCPU;
    // SequenceAlignment::alignSequenceGPU(request, &responseGPU);
    // SequenceAlignment::alignSequenceCPU(request, &responseCPU);

    // CHECK(responseCPU.score == responseGPU.score);
    // CHECK(std::string(responseCPU.alignedTextBytes, (responseCPU.alignedTextBytes + responseCPU.numAlignmentBytes)) ==
    //         std::string(responseGPU.alignedTextBytes, (responseGPU.alignedTextBytes + responseGPU.numAlignmentBytes)));
    // CHECK(std::string(responseCPU.alignedPatternBytes, (responseCPU.alignedPatternBytes + responseCPU.numAlignmentBytes)) ==
    //         std::string(responseGPU.alignedPatternBytes, (responseGPU.alignedPatternBytes + responseGPU.numAlignmentBytes)));
}

TEST_CASE("Very long (~70k) Protein alignment")
{
    //// Comment out to no spend too much time for tests.
    // const int argc = 7;
    // const char *argv[argc] = {"./alignSequence", "--protein", "--gap-penalty", "7", "--global",
    //                             "data/protein/qbpln50.txt", "data/protein/mutated_qbpln50.txt"};
    // SequenceAlignment::Request request = {};
    // parseArguments(argc, argv, &request);

    // SequenceAlignment::Response responseGPU;
    // SequenceAlignment::Response responseCPU;
    // SequenceAlignment::alignSequenceGPU(request, &responseGPU);
    // SequenceAlignment::alignSequenceCPU(request, &responseCPU);

    // CHECK(responseCPU.score == responseGPU.score);
    // CHECK(std::string(responseCPU.alignedTextBytes, (responseCPU.alignedTextBytes + responseCPU.numAlignmentBytes)) ==
    //         std::string(responseGPU.alignedTextBytes, (responseGPU.alignedTextBytes + responseGPU.numAlignmentBytes)));
    // CHECK(std::string(responseCPU.alignedPatternBytes, (responseCPU.alignedPatternBytes + responseCPU.numAlignmentBytes)) ==
    //         std::string(responseGPU.alignedPatternBytes, (responseGPU.alignedPatternBytes + responseGPU.numAlignmentBytes)));
}

