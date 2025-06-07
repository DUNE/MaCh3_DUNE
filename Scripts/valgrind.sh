#!/bin/bash

#valgrind --tool=memcheck --leak-check=full ${MACH3}/bin/nue

valgrind --tool=callgrind --dump-instr=yes --simulate-cache=yes --collect-jumps=yes ${MACH3}/bin/nue
