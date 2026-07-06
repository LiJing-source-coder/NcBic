# Chandler_1753353634222 # Set the version and distribution
VER=1.0
DIST=NcBic$(VER)
PROGS=NcBic
SRCS=struct.cpp read_array.cpp get_options.cpp fib.cpp construct_ranking_matrix.cpp enum_gene_pvalue.cpp make_graph.cpp cluster.cpp write_block.cpp checkpoint.cpp main.cpp
OBJS=$(SRCS:.cpp=.o)
CC=g++

# Set environment variable R_HOME (good practice, though not always required if paths are absolute)
#export R_HOME=/usr/lib/R

# --- Compiler and Linker Flags ---

# Base C++ flags for most files.
# -std=c++11 is the modern equivalent of -std=c++0x
# -g adds debugging info, -Wall enables all warnings, -I. lets you #include "header.h"
CPPFLAGS= -std=c++11 -DVER=$(VER) -O3 -march=native  

# Flags for files that need OpenMP
OMPFLAGS= -fopenmp

# Flags for linking the final executable.
# We include all paths and libraries needed for the final program here.
LDFLAGS= -lm -fopenmp -I/usr/include

# Libraries to link against (-l flags)
LDLIBS= -lR

# Special include flags for cluster.cpp compilation
CLUSTER_CPPFLAGS= -I/usr/include

# for Boost
Boost_CPPFLAGS= -I/usr/include

# --- Rules ---

# Default target: build the program
all: $(PROGS)

# Rule to link the final executable from all object files
${PROGS}: $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

# --- Special Compilation Rules ---

# Special rule for enum_gene_pvalue.cpp, which needs Boost
enum_gene_pvalue.o: enum_gene_pvalue.cpp
	$(CC) -c $< -o $@ $(CPPFLAGS) $(Boost_CPPFLAGS)

# Special rule for cluster.cpp, which needs OpenMP,Boost and R-specific include paths
cluster.o: cluster.cpp
	$(CC) -c $< -o $@ $(CPPFLAGS) $(OMPFLAGS) $(CLUSTER_CPPFLAGS)

# Special rule for make_graph.cpp, which needs OpenMP
make_graph.o: make_graph.cpp
	$(CC) -c $< -o $@ $(CPPFLAGS) $(OMPFLAGS)

# Special rule for struct.cpp, which needs Boost
struct.o: struct.cpp
	$(CC) -c $< -o $@ $(CPPFLAGS) $(Boost_CPPFLAGS)


# --- Generic Compilation Rule ---

# Generic rule for all other .cpp files.
# Make will use the specific rules above when they match, and this one for everything else.
# The `-c` flag means "compile, but do not link".
%.o: %.cpp
	$(CC) -c $< -o $@ $(CPPFLAGS)

# --- Cleanup and Distribution ---

clean:
	rm -f $(PROGS)
	rm -f *.o
	rm -f *.blocks

dist:
	$(MAKE) clean
	cd .. && tar czvf $(DIST).tar.gz $(DIST)/


