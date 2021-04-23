.PHONY=all tests clean obj update
CXX=g++
CC=gcc

MAKE?=make

GMATCH=$(findstring g++,$(CXX))
GIT_VERSION := $(shell git describe --abbrev=4 --always)

WARNINGS=-Wall -Wextra -Wno-char-subscripts \
		 -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
		 -Wformat -Wcast-align -Wno-unused-function -Wno-unused-parameter \
		 -pedantic -DUSE_PDQSORT -Wunused-variable -Wno-attributes -Wno-cast-align \
        -Wno-ignored-attributes -Wno-missing-braces
EXTRA?=
INCPLUS?=
EXTRA_LD?=
DBG?=-DNDEBUG
OS:=$(shell uname)
FLAGS=

OPT_MINUS_OPENMP= -O3 -funroll-loops\
	  -pipe -fno-strict-aliasing -march=native -mpclmul $(FLAGS) $(EXTRA)
OPT=$(OPT_MINUS_OPENMP) -fopenmp
XXFLAGS=-fno-rtti
CXXFLAGS=$(OPT) $(XXFLAGS) -std=c++14 $(WARNINGS) -DBONSAI_VERSION=\"$(GIT_VERSION)\"
CXXFLAGS_MINUS_OPENMP=$(OPT_MINUS_OPENMP) $(XXFLAGS) -std=c++1z $(WARNINGS) -Wno-cast-align -Wno-gnu-zero-variadic-macro-arguments -DBONSAI_VERSION=\"$(GIT_VERSION)\"
CCFLAGS=$(OPT) $(CFLAGS) -std=c11 $(WARNINGS) -DBONSAI_VERSION=\"$(GIT_VERSION)\"
LIB=-lz
LD=-L. $(EXTRA_LD)

ifneq (,$(findstring g++,$(CXX)))
	ifeq ($(shell uname),Darwin)
		ifeq (,$(findstring clang,$(CXX)))
			POPCNT_CXX:=clang
		else
			POPCNT_CXX:=$(CXX)
		endif
	else
		POPCNT_CXX:=$(CXX)
	endif
endif

OBJS=$(patsubst %.c,%.o,$(wildcard src/*.c) klib/kthread.o) $(patsubst %.cpp,%.o,$(wildcard src/*.cpp)) klib/kstring.o clhash.o
DOBJS=$(patsubst %.c,%.do,$(wildcard src/*.c) klib/kthread.o) $(patsubst %.cpp,%.do,$(wildcard src/*.cpp)) klib/kstring.o clhash.o
ZOBJS=$(patsubst %.c,%.zo,$(wildcard src/*.c) klib/kthread.o) $(patsubst %.cpp,%.zo,$(wildcard src/*.cpp)) klib/kstring.o clhash.o

TEST_OBJS=$(patsubst %.cpp,%.o,$(wildcard test/*.cpp))
ZTEST_OBJS=$(patsubst %.cpp,%.zo,$(wildcard test/*.cpp))

EXEC_OBJS=$(patsubst %.cpp,%.o,$(wildcard bin/*.cpp))
ZW_OBJS=$(patsubst %.c,%.o,$(wildcard zstd/zlibWrapper/*.c)) libzstd.a

EX=$(patsubst bin/%.cpp,bin/%,$(wildcard bin/*.cpp)) fsetsketcher fopsetsketcher opsetsketcher
DEX=$(patsubst %_d,%,$(EX))
STATEX=$(patsubst %,%_s,$(EX))
STATEXZ=$(patsubst %,%_sz,$(EX))


ZSTD_INCLUDE_DIRS=zstd/zlibWrapper zstd/lib/common zstd/lib
ZSTD_INCLUDE=$(patsubst %,-I%,$(ZSTD_INCLUDE_DIRS))
ZFLAGS=-DZWRAP_USE_ZSTD=1
ZCOMPILE_FLAGS= $(ZFLAGS)
ALL_ZOBJS=$(ZOBJS) $(ZW_OBJS)
INCLUDE=-Iclhash/include -I. -I.. -Ihll/include -Ihll/libpopcnt -I.. -Iinclude -Icircularqueue $(ZSTD_INCLUDE) $(INCPLUS) -Ihll/vec -Ihll -Ipdqsort -Iinclude/bonsai -Iinclude \
    -Ihll/vec/blaze



all: $(OBJS) $(EX) $(ZW_OBJS) unit lib/libbonsai.a $(DEX)
ex: $(EX)


update:
	for i in circularqueue clhash hll klib kspp libpopcnt linear pdqsort tinythreadpp; \
        do cd $$i && git checkout master && git pull && cd ..; done; \
	cd hll/vec && git checkout master && git pull && cd ../..

obj: $(OBJS) $(DOBJS) $(ZOBJS) $(ZW_OBJS)

libzstd.a:
	+cd zstd && $(MAKE) lib && cp lib/libzstd.a .. && cd ..

clhash.o: clhash/src/clhash.c
	ls $@ 2>/dev/null || mv clhash/clhash.o . 2>/dev/null || (cd clhash && $(MAKE) && cd .. && ln -s clhash/clhash.o .)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -DNDEBUG -c $< -o $@ $(LIB)

%.do: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -g -c $< -o $@ $(LIB)

%.zo: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@ $(LIB) $(ZCOMPILE_FLAGS)

test/%.o: test/%.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

test/%.zo: test/%.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(ZCOMPILE_FLAGS)

%.o: %.cpp $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -DNDEBUG -c $< -o $@ $(LIB)

%.do: %.cpp $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -g -c $< -o $@ $(LIB)

%.zo: %.cpp $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD)  -c $< -o $@ -DNDEBUG $(LIB) $(ZCOMPILE_FLAGS)

src/popcnt.o: src/popcnt.cpp
	$(POPCNT_CXX) $(CXXFLAGS_MINUS_OPENMP) $(DBG) $(INCLUDE) $(LD) -DNDEBUG -c $< -o $@ $(LIB)

src/popcnt.do: src/popcnt.cpp
	$(POPCNT_CXX) $(CXXFLAGS_MINUS_OPENMP) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

src/popcnt.zo: src/popcnt.cpp
	$(POPCNT_CXX) $(CXXFLAGS_MINUS_OPENMP) $(DBG) $(INCLUDE) $(LD)  -c $< -o $@ $(LIB)

adsetsketcher: bin/setsketcher.cpp clhash.o klib/kthread.o $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h) zlib/libz.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) clhash.o klib/kthread.o -DNDEBUG $< -o $@ $(LIB) zlib/libz.a -fsanitize=address -DBLAZE_SHAE
pgsetsketcher: bin/setsketcher.cpp clhash.o klib/kthread.o $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h) zlib/libz.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) clhash.o klib/kthread.o -DNDEBUG $< -o $@ $(LIB) zlib/libz.a -pg
setsketcher: bin/setsketcher.cpp clhash.o klib/kthread.o $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h) zlib/libz.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) clhash.o klib/kthread.o -DNDEBUG $< -o $@ $(LIB) zlib/libz.a
fsetsketcher: bin/setsketcher.cpp clhash.o klib/kthread.o $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h) zlib/libz.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) clhash.o klib/kthread.o -DNDEBUG $< -o $@ $(LIB) zlib/libz.a -DCSETFT=float
opsetsketcher: bin/setsketcher.cpp clhash.o klib/kthread.o $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h) zlib/libz.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) clhash.o klib/kthread.o -DNDEBUG $< -o $@ $(LIB) zlib/libz.a -DUSE_OPH
fopsetsketcher: bin/setsketcher.cpp clhash.o klib/kthread.o $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h) zlib/libz.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) clhash.o klib/kthread.o -DNDEBUG $< -o $@ $(LIB) zlib/libz.a -DCSETFT=float -DUSE_OPH
setsketchindexer: bin/setsketchindexer.cpp clhash.o klib/kthread.o $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h) zlib/libz.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) clhash.o klib/kthread.o -DNDEBUG $< -o $@ $(LIB) zlib/libz.a
fsetsketchindexer: bin/setsketchindexer.cpp clhash.o klib/kthread.o $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h) zlib/libz.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) clhash.o klib/kthread.o -DNDEBUG $< -o $@ $(LIB) zlib/libz.a -DCSETFT=float
ad%: bin/%.cpp clhash.o klib/kthread.o $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h) zlib/libz.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) clhash.o klib/kthread.o -DNDEBUG $< -o $@ $(LIB) zlib/libz.a -fsanitize=address
%: bin/%.cpp clhash.o klib/kthread.o $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h) zlib/libz.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) clhash.o klib/kthread.o -DNDEBUG $< -o $@ $(LIB) zlib/libz.a -Ihll/include
bin/errexp: bin/errexp.cpp clhash.o klib/kthread.o $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h) zlib/libz.a
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) clhash.o klib/kthread.o -DNDEBUG $< -o $@ $(LIB) zlib/libz.a -Ihll/include

bin/kmercnt: bin/kmercnt.cpp clhash.o klib/kthread.o
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) clhash.o klib/kthread.o -DNDEBUG $< -o $@ $(LIB) zlib/libz.a

bin/%_z: bin/%.cpp $(ALL_ZOBJS) $(wildcard include/bonsai/*.h) $(wildcard hll/include/sketch/*.h)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)

%_d: bin/%.cpp $(DOBJS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(DOBJS) $< -g -o $@ $(LIB)

bonsai/fahist: src/fahist.cpp $(OBJS) klib/kthread.o
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) klib/kthread.o $< -o $@ -lz

lib/libbonsai.a: $(OBJS)
	ar cr $@ $(OBJS)

zlib/libz.a:
	+cd zlib && (./configure || echo "no configure needed") && $(MAKE)

tests: clean unit

unit: $(DOBJS) $(TEST_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TEST_OBJS) $(LD) $(DOBJS) -o $@ $(LIB)

zunit: $(ZOBJS) $(ZTEST_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ZTEST_OBJS) $(LD) $(ALL_ZOBJS) -o $@ $(LIB) $(ZCOMPILE_FLAGS)

clean:
	rm -f $(ZOBJS) $(ZTEST_OBJS) $(ZW_OBJS) $(OBJS) $(DEX) $(ZEX) $(EX) $(TEST_OBJS) $(DOBJS) unit lib/*o src/*o

rolling: rolling_multk rolling_multk_sketch

mostlyclean:
	rm -f $(ZOBJS) $(ZTEST_OBJS) $(ZW_OBJS) $(OBJS) $(DEX) $(ZEX) $(EX) $(TEST_OBJS) $(DOBJS) unit lib/*o src/*o
