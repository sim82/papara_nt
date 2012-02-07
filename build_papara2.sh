#!/bin/sh

# minimal build script for papara 2.0.
# please use the supplied cmake files for anything more than just building papara 2.0.

g++ -o papara -O3 papara.cpp -I ivy_mike/src/ -I ublasJama-1.0.2.3 sequence_model.cpp parsimony.cpp ivy_mike/src/time.cpp ivy_mike/src/tree_parser.cpp ivy_mike/src/getopt.cpp ivy_mike/src/demangle.cpp ivy_mike/src/multiple_alignment.cpp ublasJama-1.0.2.3/EigenvalueDecomposition.cpp -lpthread

