#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"

#include "../SequenceAlignment.hpp"
#include "../utils.hpp"


TEST_CASE( "parseArguments")
{
    // Catch stderr to string to test error messages.
    std::stringstream buffer;
    std::streambuf *old = std::cerr.rdbuf(buffer.rdbuf());

    SECTION("usage")
    {
        const int argc = 1;
        const char *argv[argc] = { "./alignSequence"};
        parseArguments(argc, argv);

        std::string stderrString = buffer.str();
        REQUIRE(stderrString == SequenceAlignment::USAGE);
    }

    SECTION("no or empty sequence files")
    {
        const int argc = 3;
        const char *argv[argc] = { "./alignSequence", "-p", "-c"};
        parseArguments(argc, argv);

        std::string expectedMsg = "textSequence or patternSequence not read.\n" + SequenceAlignment::USAGE;
        std::string stderrString = buffer.str();
        REQUIRE(stderrString == expectedMsg);

        REQUIRE(SequenceAlignment::deviceType == SequenceAlignment::programArgs::CPU);
        REQUIRE(SequenceAlignment::sequenceType == SequenceAlignment::programArgs::PROTEIN);
    }

    // Restore original buffer.
    std::cout.rdbuf(old);
}
