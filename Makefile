# Compiler settings
CXX = g++
CXXFLAGS = -std=c++17 -Wall

# Targets
all: KDTree

KDTree: KDTree.o
	ar rvs libKDTree.a KDTree.o

KDTree.o: KDTree.cpp KDTree.hpp
	$(CXX) $(CXXFLAGS) -c KDTree.cpp

install: KDTree
	mkdir -p lib
	mv libKDTree.a lib/

clean:
	rm -f KDTree.o libKDTree.a
	rm -rf lib
