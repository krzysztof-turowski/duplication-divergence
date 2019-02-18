#!/bin/bash

# Script for libraries configuration

if [ "$#" -gt 0 ]; then
  for lib in "$@"
  do
    if [[ "$lib" == "nauty" ]]; then
      cd lib/nauty
      bash -e ./configure
      if [ $? -ne 0 ]; then
        echo "ERROR: couldn't configure nauty"
        exit
      fi
      make
      if [ $? -ne 0 ]; then
        echo "ERROR: couldn't build nauty"
      fi
      cd ../..
    fi
  done
else
  echo "ERROR: no libraries provided"
fi
