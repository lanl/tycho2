#ifndef __TYPEDEF_HH__
#define __TYPEDEF_HH__

#include <cstdint>
#include <cinttypes>
#include <limits>

#define UNUSED_VARIABLE(x) (void)(x)


typedef uint64_t UINT;
const UINT UINT_MAX = std::numeric_limits<uint64_t>::max();

#endif
