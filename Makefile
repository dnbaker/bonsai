.PHONY=all tests clean obj
CXX=g++
CC=gcc

GMATCH=$(findstring g++,$(CXX))

CLHASH_CHECKOUT = "&& git checkout master"
WARNINGS=-Wall -Wextra -Wno-char-subscripts \
		 -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
		 -Wformat -Wcast-align -Wno-unused-function -Wno-unused-parameter
		 # -pedantic
DBG:= -D_GLIBCXX_DEBUG # -DNDEBUG # -fno-inline
OS:=$(shell uname)
FLAGS=

ifeq (,$(findstring g++,$(CXX)))
	ifeq ($(shell uname),Darwin)
		ifeq (,$(findstring clang,$(CXX)))
			FLAGS := $(FLAGS) -Wa,-q
			CLHASH_CHECKOUT := "&& git checkout mac"
		endif
	endif
endif

OPT:= -O3 -funroll-loops -ffast-math \
	  -fopenmp \
	  -pipe -fno-strict-aliasing -march=native -mpclmul $(FLAGS) # -DUSE_PAR_HELPERS
XXFLAGS=-fno-rtti
CXXFLAGS=$(OPT) $(XXFLAGS) -std=c++1z $(WARNINGS)
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
	cd hll && make CXX=$(CXX) && cp libhll.a ..

obj: $(OBJS) $(EXEC_OBJS) libhll.a

clhash.o: clhash/src/clhash.c
	cp clhash/clhash.o . || (cd clhash $(CLHASH_CHECKOUT) && make && cd .. && mv clhash/clhash.o .)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@ $(LIB)

test/%.o: test/%.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

%.o: %.cpp
	echo "match: " $(GMATCH) && \
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

%: src/%.o $(OBJS) libhll.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(OBJS) $< -o $@ $(LIB)


tests: clean unit

unit: $(OBJS) $(TEST_OBJS) libhll.a
	#$(CXX) $(CXXFLAGS) $(INCLUDE) $(TEST_OBJS) $(LD) $(OBJS) -o $@ $(LIB)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) test/test_hll.o test/test_main.o $(LD) $(OBJS) -o $@ $(LIB)


clean:
	rm -f $(EXEC_OBJS) $(OBJS) $(EX) $(TEST_OBJS) unit lib/*o src/*o libhll.a \
	&& cd hll && make clean && cd ..

mostlyclean:
	rm -f $(EXEC_OBJS) $(OBJS) $(EX) $(TEST_OBJS) unit lib/*o src/*o
