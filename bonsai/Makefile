.PHONY=all tests clean obj update
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
ifndef EXTRA_LD
	EXTRA_LD:=
endif
DBG:=
OS:=$(shell uname)
FLAGS=

OPT_MINUS_OPENMP= -O3 -funroll-loops\
	  -pipe -fno-strict-aliasing -march=native -mpclmul $(FLAGS) $(EXTRA)
OPT=$(OPT_MINUS_OPENMP) -fopenmp
XXFLAGS=-fno-rtti
CXXFLAGS=$(OPT) $(XXFLAGS) -std=c++17 $(WARNINGS)
CXXFLAGS_MINUS_OPENMP=$(OPT_MINUS_OPENMP) $(XXFLAGS) -std=c++1z $(WARNINGS) -Wno-cast-align -Wno-gnu-zero-variadic-macro-arguments
CCFLAGS=$(OPT) $(CFLAGS) -std=c11 $(WARNINGS)
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

OBJS=$(patsubst %.c,%.o,$(wildcard src/*.c) ../klib/kthread.o) $(patsubst %.cpp,%.o,$(wildcard src/*.cpp)) ../klib/kstring.o clhash.o
DOBJS=$(patsubst %.c,%.do,$(wildcard src/*.c) ../klib/kthread.o) $(patsubst %.cpp,%.do,$(wildcard src/*.cpp)) ../klib/kstring.o clhash.o
ZOBJS=$(patsubst %.c,%.zo,$(wildcard src/*.c) ../klib/kthread.o) $(patsubst %.cpp,%.zo,$(wildcard src/*.cpp)) ../klib/kstring.o clhash.o

TEST_OBJS=$(patsubst %.cpp,%.o,$(wildcard test/*.cpp))
ZTEST_OBJS=$(patsubst %.cpp,%.zo,$(wildcard test/*.cpp))

EXEC_OBJS=$(patsubst %.cpp,%.o,$(wildcard bin/*.cpp))
ZW_OBJS=$(patsubst %.c,%.o,$(wildcard ../zstd/zlibWrapper/*.c)) libzstd.a

EX=$(patsubst bin/%.cpp,%,$(wildcard bin/*.cpp))
DEX=$(patsubst %_d,%,$(EX))
ZEX=$(patsubst %_z,%,$(EX))
STATEX=$(patsubst %,%_s,$(EX))
STATEXZ=$(patsubst %,%_sz,$(EX))


ZSTD_INCLUDE_DIRS=../zstd/zlibWrapper ../zstd/lib/common ../zstd/lib
ZSTD_INCLUDE=$(patsubst %,-I%,$(ZSTD_INCLUDE_DIRS))
ZFLAGS=-DZWRAP_USE_ZSTD=1
ZCOMPILE_FLAGS= $(ZFLAGS) -lzstd
ALL_ZOBJS=$(ZOBJS) $(ZW_OBJS)
INCLUDE=-Iclhash/include -I. -I.. -I../libpopcnt -I.. -Iinclude -I../circularqueue $(ZSTD_INCLUDE) $(INCPLUS) -I../filtered_hll/hll/vec -I../filtered_hll -I../pdqsort



all: update $(OBJS) $(EX) $(ZW_OBJS) unit lib/libbonsai.a $(DEX) $(ZEX)

#libhll.a:
#	cd hll && make CXX=$(CXX) && cp libhll.a ..

update:
	for i in circularqueue clhash cppitertools filtered_hll klib kraken kspp libpopcnt linear pdqsort tinythreadpp zstd; \
        do cd ../$$i && git checkout master && git pull && cd ../bonsai; done; \
    cd ../filtered_hll/hll &&  git checkout master && git pull && cd .. && cd cruel_bloom && git checkout master && git pull && cd ../..

obj: $(OBJS) $(DOBJS) $(ZOBJS) $(ZW_OBJS)

libzstd.a:
	cd ../zstd && make && mv lib/libzstd.a ../bonsai && cd ../bonsai

clhash.o: ../clhash/src/clhash.c
	ls $@ 2>/dev/null || mv ../clhash/clhash.o . 2>/dev/null || (cd ../clhash $(CLHASH_CHECKOUT) && make && cd ../bonsai && ln -s ../clhash/clhash.o .)

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

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -DNDEBUG -c $< -o $@ $(LIB)

%.do: %.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) -g -c $< -o $@ $(LIB)

%.zo: %.cpp
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD)  -c $< -o $@ -DNDEBUG $(LIB) $(ZCOMPILE_FLAGS)

src/popcnt.o: src/popcnt.cpp
	$(POPCNT_CXX) $(CXXFLAGS_MINUS_OPENMP) $(DBG) $(INCLUDE) $(LD) -DNDEBUG -c $< -o $@ $(LIB)

src/popcnt.do: src/popcnt.cpp
	$(POPCNT_CXX) $(CXXFLAGS_MINUS_OPENMP) $(DBG) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

src/popcnt.zo: src/popcnt.cpp
	$(POPCNT_CXX) $(CXXFLAGS_MINUS_OPENMP) $(DBG) $(INCLUDE) $(LD)  -c $< -o $@ $(LIB)

%: bin/%.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(OBJS) -DNDEBUG $< -o $@ $(LIB)

%_s: bin/%.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(OBJS) -static-libstdc++ -static-libgcc -DNDEBUG $< -o $@ $(LIB)


%_z: bin/%.cpp $(ALL_ZOBJS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)

%_sz: bin/%.cpp $(ALL_ZOBJS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(ALL_ZOBJS) -static-libstdc++ -static-libgcc -DNDEBUG $< -o $@ $(ZCOMPILE_FLAGS) $(LIB)

%_d: bin/%.cpp $(DOBJS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) $(DOBJS) $< -g -o $@ $(LIB)

bonsai/fahist: src/fahist.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(DBG) $(INCLUDE) $(LD) ../klib/kthread.o $< -o $@ -lz

lib/libbonsai.a: $(OBJS)
	ar cr $@ $(OBJS)

tests: clean unit

unit: $(DOBJS) $(TEST_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TEST_OBJS) $(LD) $(DOBJS) -o $@ $(LIB)

zunit: $(ZOBJS) $(ZTEST_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(ZTEST_OBJS) $(LD) $(ALL_ZOBJS) -o $@ $(LIB) $(ZCOMPILE_FLAGS)

clean:
	rm -f $(ZOBJS) $(ZTEST_OBJS) $(ZW_OBJS) $(OBJS) $(DEX) $(ZEX) $(EX) $(TEST_OBJS) $(DOBJS) unit lib/*o src/*o \
	&& cd ../hll && make clean && cd ../bonsai

mostlyclean:
	rm -f $(ZOBJS) $(ZTEST_OBJS) $(ZW_OBJS) $(OBJS) $(DEX) $(ZEX) $(EX) $(TEST_OBJS) $(DOBJS) unit lib/*o src/*o
