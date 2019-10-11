
CXX?=g++
SUF=`$(PYCONF) --extension-suffix`
PYCONF?=python3-config
UNDEFSTR=-undefined dynamic_lookup
WARNINGS=-Wall -Wextra -Wno-char-subscripts \
         -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
         -Wformat -Wcast-align -Wno-unused-function -Wno-unused-parameter \
         -pedantic -Wunused-variable -Wno-attributes -Wno-ignored-attributes -Wpedantic

FLAGS=-O1 -funroll-loops -pipe -march=native -I. -Ivec/blaze -Ivec -Ipybind11/include -Iinclude -fpic -Wall $(WARNINGS) \
     -fno-strict-aliasing \
      -DXXH_INLINE_ALL  \
      -Wno-attributes -Wno-pragmas -Wno-ignored-qualifiers

CXXFLAGS=$(FLAGS) -Wreorder
INCLUDES=-I. -I.. -I../hll -I../bonsai -I../hll/vec -I../hll/libpopcnt -I../circularqueue

%.cpython.so: %.cpp
	$(CXX) $(UNDEFSTR) $(INCLUDES) -fopenmp -O3 -Wall $(CXXFLAGS) -shared $(STD) -fPIC `python3 -m pybind11 --includes` $< -o $*$(SUF) -lz && \
    ln -fs $*$(SUF) $@
