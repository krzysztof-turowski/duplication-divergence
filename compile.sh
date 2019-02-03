#!/bin/bash

# Creates C++-executable ./FILENAME from provided code files and libraries
# Usage: ./compile.sh FILENAME.cpp FILENAMES

if [ "$#" -gt 0 ]; then
  EXE=${1%.*}
  g++ $@ -O3 -std=c++11 -pthread -Wall -Wextra -Wstrict-aliasing -Wpedantic -fmax-errors=5 -Werror -Wunreachable-code -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option -o ./$EXE
  if [ $? -ne 0 ]; then
    echo "ERROR: couldn't compile ./$EXE"
    exit
  fi
  echo "Successfully generated ./$EXE"
else
  echo "ERROR: no input files given"
fi
