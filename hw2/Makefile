CXX = g++
CXXFLAGS = -std=c++11 -O3 -g

# Default target
all: rasterizer

rasterizer: *.cpp
	$(CXX) $^ $(CXXFLAGS) -o $@

clean:
	rm -f rasterizer
