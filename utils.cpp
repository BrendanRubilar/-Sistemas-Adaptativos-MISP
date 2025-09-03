/** Author: LELE */

#ifndef UTILS
#define UTILS

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

inline void validate_input(int argc, char *argv[], std::int64_t& runs,
    std::int64_t& lower, std::int64_t& upper, std::int64_t& step)
{
    if (argc != 6) {
        std::cerr << "Usage: <filename> <RUNS> <LOWER> <UPPER> <STEP>" << std::endl;
        std::cerr << "<filename> is the name of the file where performance data will be written." << std::endl;
        std::cerr << "It is recommended for <filename> to have .csv extension and it should not previously exist." << std::endl;
        std::cerr << "<RUNS>: numbers of runs per test case: should be >= 32." << std::endl;
        std::cerr << "<LOWER> <UPPER> <STEP>: range of test cases." << std::endl;
        std::cerr << "These should all be positive." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Read command line arguments
    try {
        runs = std::stoll(argv[2]);
        lower = std::stoll(argv[3]);
        upper = std::stoll(argv[4]);
        step = std::stoll(argv[5]);
    } catch (std::invalid_argument const& ex) {
        std::cerr << "std::invalid_argument::what(): " << ex.what() << std::endl;
        std::exit(EXIT_FAILURE);
    } catch (std::out_of_range const& ex) {
        std::cerr << "std::out_of_range::what(): " << ex.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Validate arguments
    if (runs < 4) {
        std::cerr << "<RUNS> must be at least 4." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (step <= 0 or lower <= 0 or upper <= 0) {
        std::cerr << "<STEP>, <LOWER> and <UPPER> have to be positive." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (lower > upper) {
        std::cerr << "<LOWER> must be at most equal to <UPPER>." << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

inline void display_progress(std::int64_t u, std::int64_t v)
{
    const double progress = u / double(v);
    const std::int64_t width = 70;
    const std::int64_t p = width * progress;
    std::int64_t i;

    std::cerr << "\033[1m[";
    for (i = 0; i < width; i++) {
        if (i < p)
            std::cerr << "=";
        else if (i == p)
            std::cerr << ">";
        else
            std::cerr << " ";
    }
    std::cerr << "] " << std::int64_t(progress * 100.0) << "%\r\033[0m";
    std::cerr.flush();
}

#endif
