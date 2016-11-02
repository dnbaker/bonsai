.PHONY=all tests clean obj unit_tests
CXX=g++
WARNINGS=-Wall -pedantic -Wextra -Wno-char-subscripts \
         -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
         -Wformat -Wcast-align
OPT=-O3 -funroll-loops
CXXFLAGS=$(OPT) -std=c++17 $(WARNINGS)
LIB=-lz -pthread
LD=-L. -L/usr/lib/gcc/x86_64-redhat-linux/6.2.1/


OBJS=lib/encoder.o lib/spacer.o lib/cms.o \
     lib/feature_min.o lib/util.o lib/db.o

TEST_OBJS=test/test_encoding.o \
          test/test_feature_min.o \
		  test/test_hll.o \
		  test/test_main.o \
		  test/test_qmap.o \
		  test/test_util.o

EXEC_OBJS=src/taxmap.o

EX=taxmap mapmake

HEADERS=lib/encoder.h lib/kmerutil.h lib/spacer.h lib/misc.h \
        lib/cms.h lib/kseq_declare.h lib/feature_min.h lib/hll.h lib/hash.h lib/db.h


INCLUDE=-I. -Ihtslib -Ilib -Ijellyfish/include/

all: $(OBJS) $(EX)

obj: $(OBJS) $(EXEC_OBJS) libcityhash.a
libcityhash.a:
	+cd cityhash/ && ./configure && make && make install prefix=$$PWD && cp lib/libcityhash.a ../ && cd ../


%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LD) -c $< -o $@ $(LIB)

%: src/%.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LD) $(OBJS) $< -o $@ $(LIB)

tests: clean unit_tests


unit_tests: $(OBJS) $(TEST_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TEST_OBJS) $(LD) $(OBJS) -o $@ $(LIB) ;\
	 ./unit_tests

clean:
	rm -f $(EXEC_OBJS) $(OBJS) $(EX) $(TEST_OBJS)
