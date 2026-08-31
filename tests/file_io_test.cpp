#include <stdexcept>

#include <catch2/catch_test_macros.hpp>

#include "include/FILE_IO.h"

// Regression test for a real crash: FILE_IO used to assert(fileout.is_open())
// on a bad output path, which aborted the whole process instead of failing
// gracefully -- and our CMAKE_CXX_FLAGS_RELEASE override ("-O2" only, no
// -DNDEBUG) means that assert never gets compiled out even in Release
// builds. Now it throws a catchable, actionable error instead.
TEST_CASE("FILE_IO::WriteTimeStep throws instead of aborting when the output directory doesn't exist", "[file_io]") {
    FILE_IO FileIO;
    vector<double> timeStep = {0.1, 0.2};

    REQUIRE_THROWS_AS(
        FileIO.WriteTimeStep("/this/directory/does/not/exist", "test", timeStep),
        runtime_error);
}
