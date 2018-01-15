.PHONY=all tests clean obj
CXX=g++
CC=gcc

GMATCH=$(findstring g++,$(CXX))

CLHASH_CHECKOUT = "&& git checkout master"
WARNINGS=-Wall -Wextra -Wno-char-subscripts \
		 -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
		 -Wformat -Wcast-align -Wno-unused-function -Wno-unused-parameter \
		 -pedantic -DUSE_PDQSORT -Wunused-variable \
		-Wduplicated-branches -Wdangling-else  # -Wconversion
ifndef EXTRA
	EXTRA:=
endif
ifndef INCPLUS
	INCPLUS:=
endif
DBG:=
OS:=$(shell uname)
FLAGS=

OPT:= -O3 -funroll-loops \
	  -fopenmp \
	  -pipe -fno-strict-aliasing -march=native -mpclmul $(FLAGS) $(EXTRA) -DHLL_HEADER_ONLY
XXFLAGS=-fno-rtti
CXXFLAGS=$(OPT) $(XXFLAGS) -std=c++17 $(WARNINGS)
CCFLAGS=$(OPT) $(CFLAGS) -std=c11 $(WARNINGS)
LIB=-lz #-lhll
LD=-L.

OBJS=$(patsubst %.c,%.o,$(wildcard lib/*.c) klib/kthread.o) $(patsubst %.cpp,%.o,$(wildcard lib/*.cpp)) klib/kstring.o clhash.o

TEST_OBJS=$(patsubst %.cpp,%.o,$(wildcard test/*.cpp))

EXEC_OBJS=$(patsubst %.cpp,%.o,$(wildcard src/*.cpp))

EX=$(patsubst src/%.cpp,%,$(wildcard src/*.cpp))

HEADERS=lib/encoder.h lib/kmerutil.h lib/spacer.h lib/misc.h \
		lib/kseq_declare.h lib/feature_min.h hll/hll.h lib/hash.h lib/db.h

INCLUDE=-I. -Ilib $(INCPLUS)

all: $(OBJS) $(EX) unit

#libhll.a:
#	cd hll && make CXX=$(CXX) && cp libhll.a ..

obj: $(OBJS)

clhash.o: clhash/src/clhash.c
	ls $@ || ln -s clhash/clhash.o . || (cd clhash $(CLHASH_CHECKOUT) && make && cd .. && ln -s clhash/clhash.o .)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@ $(LIB)

test/%.o: test/%.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

%: src/%.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(OBJS) $< -o $@ $(LIB)

fahist: src/fahist.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) klib/kthread.o $< -o $@ -lz

tests: clean unit

unit: $(OBJS) $(TEST_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TEST_OBJS) $(LD) $(OBJS) -o $@ $(LIB)
	# $(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) test/test_hll.o test/test_main.o $(LD) $(OBJS) -o $@ $(LIB)


clean:
	rm -f $(OBJS) $(EX) $(TEST_OBJS) unit lib/*o src/*o \
	&& cd hll && make clean && cd ..

mostlyclean:
	rm -f $(OBJS) $(EX) $(TEST_OBJS) unit lib/*o src/*o
