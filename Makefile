.PHONY=all tests clean obj
CXX=g++
CC=gcc
WARNINGS=-Wall -Wextra -Wno-char-subscripts \
         -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
         -Wformat -Wcast-align -Wno-unused-function -Wno-unused-parameter \
         # -pedantic
DBG= -DNDEBUG # -fno-inline
OPT:= -O3 -funroll-loops -fno-asynchronous-unwind-tables -ffast-math \
      -pipe -fno-strict-aliasing -march=native -mpclmul # -Wa,-q #-flto#  -mavx512f # -msse2
OS:=$(shell uname)
ifeq ($(OS),Darwin)
	OPT := $(OPT) -Wa,-q
endif
XXFLAGS=-fno-rtti
CXXFLAGS=$(OPT) $(XXFLAGS) -std=c++14 $(WARNINGS)
CCFLAGS=$(OPT) -std=c11 $(WARNINGS)
LIB=-lz -pthread -lhll -lcrypto
LD=-L.

OBJS=$(patsubst %.c,%.o,$(wildcard lib/*.c) klib/kthread.o) $(patsubst %.cpp,%.o,$(wildcard lib/*.cpp)) klib/kstring.o clhash.o

TEST_OBJS=$(patsubst %.cpp,%.o,$(wildcard test/*.cpp))

EXEC_OBJS=$(patsubst %.cpp,%.o,$(wildcard src/*.cpp))

EX=$(patsubst src/%.o,%,$(EXEC_OBJS))

HEADERS=lib/encoder.h lib/kmerutil.h lib/spacer.h lib/misc.h \
        lib/cms.h lib/kseq_declare.h lib/feature_min.h hll/hll.h lib/hash.h lib/db.h

INCLUDE=-I. -Ilib



all: $(OBJS) $(EX) unit

libhll.a:
	cd hll && make && cp libhll.a ..

obj: $(OBJS) $(EXEC_OBJS) libhll.a

clhash.o: clhash/src/clhash.c
	cd clhash && git checkout mac && make && cd .. && mv clhash/clhash.o .

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@ $(LIB)

test/%.o: test/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

%: src/%.o $(OBJS) libhll.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(OBJS) $< -o $@ $(LIB)


tests: clean unit

unit: $(OBJS) $(TEST_OBJS) libhll.a
	#$(CXX) $(CXXFLAGS) $(INCLUDE) $(TEST_OBJS) $(LD) $(OBJS) -o $@ $(LIB)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TEST_OBJS) $(LD) $(OBJS) -o $@ $(LIB)


clean:
	rm -f $(EXEC_OBJS) $(OBJS) $(EX) $(TEST_OBJS) unit lib/*o src/*o libhll.a \
	&& cd hll && make clean && cd ..

mostlyclean: clean
