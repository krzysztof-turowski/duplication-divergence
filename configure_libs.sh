#!/bin/bash

# Script for libraries configuration

if [ "$#" -gt 0 ]; then
  CURRENT_DIR=`pwd`
  for lib in "$@"
  do
    if [[ "$lib" == "glpk" ]]; then
      GLPK_LIB=$(find -L /usr/lib $HOME -path */glpk/lib/libglpk.a -print0 2>/dev/null | head -n 1)
      if [ $? -ne 0 ]; then
        # TODO: find GLPK sources nad install
        echo "ERROR: GLPK not found"
        exit
      fi
      GLPK_DIR=$(dirname $(dirname $GLPK_LIB))
      echo "export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$GLPK_DIR/include" >> $HOME/.bashrc
      ld -lglpk 2>/dev/null
      if [ $? -eq 0 ]; then
        echo "GLPK already present"
      else
        echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GLPK_DIR/lib" >> $HOME/.bashrc
        echo "export LIBRARY_PATH=$LIBRARY_PATH:$GLPK_DIR/lib" >> $HOME/.bashrc
      fi
    fi
    if [[ "$lib" == "gurobi" ]]; then
      if [[ -z "$GUROBI_HOME" ]]; then
        GUROBI_LIB=$(find -L /usr/lib $HOME $CURRENT_DIR -path */gurobi/lib/libgurobi_c++.a -print0 2>/dev/null | head -n 1)
        if [ $? -ne 0 ]; then
          echo "ERROR: Gurobi not found"
          exit
        fi
        GUROBI_DIR=$(dirname $(dirname $GUROBI_LIB))
        echo "export GUROBI_HOME=$GUROBI_DIR" >> $HOME/.bashrc
      else
        GUROBI_DIR=$GUROBI_HOME
      fi
      if [[ $PATH != *"$GUROBI_DIR"* ]]; then
        echo "export PATH=$PATH:$GUROBI_DIR/bin" >> $HOME/.bashrc
      fi
      if [[ -z "$GRB_LICENSE_FILE" ]]; then
        GUROBI_LICENSE=$(find -L $HOME $CURRENT_DIR -name gurobi.lic -print0 2>/dev/null | head -n 1)
        if [ $? -ne 0 ]; then
          echo "ERROR: Gurobi license not found"
          exit
        fi
        echo "export GRB_LICENSE_FILE=$GUROBI_LICENSE" >> $HOME/.bashrc
      fi
      cd $GUROBI_DIR/src/build
      make clean
      make c++
      if [ $? -ne 0 ]; then
        echo "ERROR: couldn't build Gurobi"
        exit
      fi
      ln -sf $GUROBI_DIR/src/build/libgurobi_c++.a $GUROBI_DIR/lib/libgurobi_c++.a
      echo "export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$GUROBI_DIR/include" >> $HOME/.bashrc
      ld -lgurobi_c++ -lgurobi81 2>/dev/null
      if [ $? -eq 0 ]; then
        echo "Gurobi already present"
      else
        echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GUROBI_DIR/lib" >> $HOME/.bashrc
        echo "export LIBRARY_PATH=$LIBRARY_PATH:$GUROBI_DIR/lib" >> $HOME/.bashrc
      fi
    fi
    if [[ "$lib" == "nauty" ]]; then
      cd $CURRENT_DIR/lib/nauty
      chmod +x ./configure
      ./configure
      if [ $? -ne 0 ]; then
        echo "ERROR: couldn't configure nauty"
        exit
      fi
      make
      if [ $? -ne 0 ]; then
        echo "ERROR: couldn't build nauty"
        exit
      fi
    fi
    if [[ "$lib" == "networkit" ]]; then
      mkdir $CURRENT_DIR/lib/networkit/build_lib
      cd $CURRENT_DIR/lib/networkit/build_lib
      cmake ..
      if [ $? -ne 0 ]; then
        echo "ERROR: couldn't create Makefile for NetworKit"
        exit
      fi
      make
      if [ $? -ne 0 ]; then
        echo "ERROR: couldn't build NetworKit"
        exit
      fi
      echo "export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$CURRENT_DIR/lib/networkit/networkit/cpp" >> $HOME/.bashrc
      ld -networkit 2>/dev/null
      if [ $? -eq 0 ]; then
        echo "NetworKit already present"
      else
        echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CURRENT_DIR/lib/networkit/build_lib" >> $HOME/.bashrc
        echo "export LIBRARY_PATH=$LIBRARY_PATH:$CURRENT_DIR/lib/networkit/build_lib" >> $HOME/.bashrc
      fi
    fi
    if [[ "$lib" == "snap" ]]; then
      cd $CURRENT_DIR/lib/snap/snap-core
      make lib
      if [ $? -ne 0 ]; then
        echo "ERROR: couldn't build SNAP"
        exit
      fi
      echo "export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$CURRENT_DIR/lib/snap/snap-core:$CURRENT_DIR/lib/snap/glib-core" >> $HOME/.bashrc
      ld -lsnap 2>/dev/null
      if [ $? -eq 0 ]; then
        echo "SNAP already present"
      else
        echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CURRENT_DIR/lib/snap/snap-core:$CURRENT_DIR/lib/snap/glib-core" >> $HOME/.bashrc
        echo "export LIBRARY_PATH=$LIBRARY_PATH:$CURRENT_DIR/lib/snap/snap-core:$CURRENT_DIR/lib/snap/glib-core" >> $HOME/.bashrc
      fi
    fi
  done
  cd $CURRENT_DIR
  printf "\nAll libraries set. Reload the shell or run 'source ~/.bashrc'.\n"
else
  echo "ERROR: no libraries provided"
fi
