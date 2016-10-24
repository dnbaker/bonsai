CXX=g++
CXXFLAGS=-O3 -flto -std=c++17 -DNDEBUG
LD=-lz -pthread
INCLUDE=-I. -Ihtslib -Ilib



%.o: %.cpp lib/encoder.h lib/kmerutil.h lib/spacer.h
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LD)

