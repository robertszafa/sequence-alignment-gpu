#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include "../SequenceAlignment.hpp"

TEST_CASE( "Parsing command line arguments", "[parseArguments]" ) {
    char* argv[10];

    argv[0] = "alignSequence";
    parseArguments(1, (const char**) argv);
    // REQUIRE(parseArguments(1, argv) == SequenceAlignment::USAGE);
}
