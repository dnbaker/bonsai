CXX=g++
CXXFLAGS=-O3 -flto -funroll-loops -std=c++17 -DNDEBUG
LD=-lz -pthread
INCLUDE=-I. -Ihtslib -Ilib -Ithird_party

OBJS=lib/encoder.o lib/kmerutil.o lib/spacer.o lib/cms.o third_party/quickfile.o \
     lib/feature_min.o

obj: $(OBJS)

%.o: %.cpp lib/encoder.h lib/kmerutil.h lib/spacer.h
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LD)

