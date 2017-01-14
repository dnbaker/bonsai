.PHONY=all tests clean obj unit_tests
CXX=g++-mp-6
CC=gcc-mp-6
WARNINGS=-Werror -Wall -pedantic -Wextra -Wno-char-subscripts \
         -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
         -Wformat -Wcast-align -Wno-unused-function -Wno-unused-parameter
DBG=# -DNDEBUG # -fno-inline
OPT= $(DBG) -O3 -funroll-loops -fno-asynchronous-unwind-tables -ffast-math \
			-pipe -fno-strict-aliasing -march=native -Wa,-q #-flto#  -mavx512f # -msse2
XXFLAGS=-fno-rtti
CXXFLAGS=$(OPT) $(XXFLAGS) -std=c++14 $(WARNINGS)
CCFLAGS=$(OPT) -std=c11 $(WARNINGS)
LIB=-lz -pthread -lhll
LD=-L.

OBJS=$(patsubst %.c,%.o,$(wildcard lib/*.c) klib/kthread.o) $(patsubst %.cpp,%.o,$(wildcard lib/*.cpp))

TEST_OBJS=$(patsubst %.cpp,%.o,$(wildcard test/*.cpp))

EXEC_OBJS=$(patsubst %.cpp,%.o,$(wildcard src/*.cpp))

EX=$(patsubst src/%.o,%,$(EXEC_OBJS))

HEADERS=lib/encoder.h lib/kmerutil.h lib/spacer.h lib/misc.h \
        lib/cms.h lib/kseq_declare.h lib/feature_min.h hll/hll.h lib/hash.h lib/db.h


INCLUDE=-I. -Ilib



all: $(OBJS) $(EX) unit_tests

libhll.a:
	cd hll && make && cp libhll.a ..

obj: $(OBJS) $(EXEC_OBJS) libhll.a

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@ $(LIB)

#%.o: %.h %.cpp
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

%: src/%.o $(OBJS) libhll.a
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LD) $(OBJS) $< -o $@ $(LIB)


tests: clean unit_tests

unit_tests: $(OBJS) $(TEST_OBJS) libhll.a
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TEST_OBJS) $(LD) $(OBJS) -o $@ $(LIB)


clean:
	rm -f $(EXEC_OBJS) $(OBJS) $(EX) $(TEST_OBJS) unit_tests lib/*o src/*o libhll.a \
	&& cd hll && make clean && cd ..

mostlyclean: clean
