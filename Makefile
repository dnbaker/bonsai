CXX=g++
CXXFLAGS=-O3 -flto -funroll-loops -std=c++17 -DNDEBUG -Wall -pedantic -Wno-char-subscripts
LD=-lz -pthread


OBJS=lib/encoder.o lib/kmerutil.o lib/spacer.o lib/cms.o \
     lib/feature_min.o lib/ncbi.o lib/util.o

EXEC_OBJS=src/taxmap.o

EX=taxmap

HEADERS=lib/encoder.h lib/kmerutil.h lib/spacer.h lib/misc.h lib/ncbi.h \
        lib/cms.h lib/kseq_declare.h lib/feature_min.h lib/hll.h lib/hash.h



INCLUDE=-I. -Ihtslib -Ilib -Ithird_party

all: $(OBJS) $(EX)

obj: $(OBJS) $(EXEC_OBJS)

lib/%.o: lib/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LD)

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LD)

%: src/%.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OBJS) $< -o $@ $(LD)

clean:
	rm -f $(EXEC_OBJS) $(OBJS) $(EX)
