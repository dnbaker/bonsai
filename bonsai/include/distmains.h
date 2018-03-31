#ifndef DIST_MAINS_H__
#define DIST_MAINS_H__
#include <fstream>
#include <omp.h>
#include "util.h"
#include "database.h"
#include "bitmap.h"
#include "setcmp.h"
#include "klib/kthread.h"
#include <sstream>

namespace bns {
using MainFnPtr = int (*) (int, char **);
int sketch_main(int, char **);
int dist_main(int, char **);
int setdist_main(int, char **);
}


#endif // #ifndef DIST_MAINS_H__
