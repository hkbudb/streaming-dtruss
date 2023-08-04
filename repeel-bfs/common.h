#ifndef COMMON_H_
#define COMMON_H_

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <cstdint>
#include <iostream>

#include <string>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>

#define ASSERT(truth) \
    if (!(truth)) { \
      std::cerr << "\x1b[1;31mASSERT\x1b[0m: " \
                << "LINE " << __LINE__ \
                << ", " << __FILE__ \
                << std::endl; \
      std::exit(EXIT_FAILURE); \
    } else

#define ASSERT_MSG(truth, msg) \
    if (!(truth)) { \
      std::cerr << "\x1b[1;31mASSERT\x1b[0m: " \
                << "LINE " << __LINE__ \
                << ", " << __FILE__ << '\n' \
                << "\x1b[1;32mINFO\x1b[0m: " << msg \
                << std::endl; \
      std::exit(EXIT_FAILURE); \
    } else

#define COMMON_MAX(valA, valB) ((valA > valB)?valA:valB)
#define COMMON_MIN(valA, valB) ((valA < valB)?valA:valB)
using namespace std;

#endif
