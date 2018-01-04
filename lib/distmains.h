#ifndef DIST_MAINS_H__
#define DIST_MAINS_H__
#include <functional>
#include <fstream>
#include <omp.h>
#include "lib/util.h"
#include "lib/database.h"
#include "lib/bitmap.h"
#include "lib/setcmp.h"
#include <sstream>

namespace emp {
using MainFnPtr = int (*) (int, char **);
int sketch_main(int, char **);
int dist_main(int, char **);
int setdist_main(int, char **);
}


#endif // #ifndef DIST_MAINS_H__
