# lazy
Lazy containers for reduced memory footprint and optional initialization

### Vector

lazy/vector.h contains a pseudo-STL implementation. While it doesn't provide the full std::vector interface, most functionality is covered.
The primary advantages are optional initialization and manual control of resizing. This avoids rounding size to the nearest power of 2, but adds some work to the user.
Using push_back instead of emplace_back causes an exponential resizing set by LAZY_PUSH_BACK_RESIZING_FACTOR. If unset, it defaults to 1.25, a rather conservative value.

This currently requires support for C++17.
