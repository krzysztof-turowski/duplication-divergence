#!/bin/bash

# Script for libraries configuration

if [ "$#" -gt 0 ]; then
  CURRENT_DIR=`pwd`
  INCLUDE_DIRS=""
  LIBRARY_DIRS=""
  for lib in "$@"
  do
    if [[ "$lib" == "glpk" ]]; then
      ld -L `echo $LD_LIBRARY_PATH | tr ':' '\n' | grep glpk | head -1` -lglpk 2>/dev/null
      if [ $? -eq 0 ]; then
        echo "GLPK already present"
      else
        echo "GLPK not present"
        if [[ "$lib" == "glpk" ]]; then
          GLPK_LIB=$(find -L /usr/lib $HOME -path */glpk/lib/libglpk.a -print0 2>/dev/null | head -1)
          if [ $? -ne 0 ]; then
            # TODO: find GLPK sources and install
            echo "ERROR: GLPK not found"
            exit
          fi
          GLPK_DIR=$(dirname $(dirname $GLPK_LIB))
          INCLUDE_DIRS="$INCLUDE_DIRS:$GLPK_DIR/include"
          LIBRARY_DIRS="$LIBRARY_DIRS:$GLPK_DIR/lib"
        fi
      fi
    fi
    if [[ "$lib" == "gurobi" ]]; then
      ld -L `echo $LD_LIBRARY_PATH | tr ':' '\n' | grep gurobi | head -1` -lgurobi_c++ -lgurobi81 2>/dev/null
      if [ $? -eq 0 ]; then
        echo "Gurobi already present"
      else
        echo "Gurobi not present"
        if [[ -z "$GUROBI_HOME" ]]; then
          GUROBI_LIB=$(find -L /usr/lib $HOME $CURRENT_DIR -path */gurobi/lib/libgurobi_c++.a -print0 2>/dev/null | head -1)
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
          echo "export PATH=\$PATH:$GUROBI_DIR/bin" >> $HOME/.bashrc
        fi
        if [[ -z "$GRB_LICENSE_FILE" ]]; then
          GUROBI_LICENSE=$(find -L $HOME $CURRENT_DIR -name gurobi.lic -print0 2>/dev/null | head -1)
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
        INCLUDE_DIRS="$INCLUDE_DIRS:$GUROBI_DIR/include"
        LIBRARY_DIRS="$LIBRARY_DIRS:$GUROBI_DIR/lib"
      fi
    fi
    if [[ "$lib" == "nauty" ]]; then
      ld -L `echo $LD_LIBRARY_PATH | tr ':' '\n' | grep nauty | head -1` -l:nauty.a 2>/dev/null
      if [ $? -eq 0 ]; then
        echo "nauty already present"
      else
        echo "nauty not present"
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
        INCLUDE_DIRS="$INCLUDE_DIRS:$CURRENT_DIR/lib/nauty"
        LIBRARY_DIRS="$LIBRARY_DIRS:$CURRENT_DIR/lib/nauty"
      fi
    fi
    if [[ "$lib" == "networkit" ]]; then
      ld -L `echo $LD_LIBRARY_PATH | tr ':' '\n' | grep networkit | head -1` -lnetworkit 2>/dev/null
      if [ $? -eq 0 ]; then
        echo "NetworKit already present"
      else
        echo "NetworKit not present"
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
        INCLUDE_DIRS="$INCLUDE_DIRS:$CURRENT_DIR/lib/networkit/networkit/cpp"
        LIBRARY_DIRS="$LIBRARY_DIRS:$CURRENT_DIR/lib/networkit/build_lib"
      fi
    fi
    if [[ "$lib" == "snap" ]]; then
      ld -L `echo $LD_LIBRARY_PATH | tr ':' '\n' | grep snap | head -1` -lsnap 2>/dev/null
      if [ $? -eq 0 ]; then
        echo "SNAP already present"
      else
        echo "SNAP not present"
        cd $CURRENT_DIR/lib/snap/snap-core
        make lib
        if [ $? -ne 0 ]; then
          echo "ERROR: couldn't build SNAP"
          exit
        fi
        INCLUDE_DIRS="$INCLUDE_DIRS:$CURRENT_DIR/lib/snap/snap-core:$CURRENT_DIR/lib/snap/glib-core"
        LIBRARY_DIRS="$LIBRARY_DIRS:$CURRENT_DIR/lib/snap/snap-core:$CURRENT_DIR/lib/snap/glib-core"
      fi
    fi
  done
  cd $CURRENT_DIR
  echo "export CPLUS_INCLUDE_PATH=\$CPLUS_INCLUDE_PATH$INCLUDE_DIRS" >> $HOME/.bashrc
  echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH$LIBRARY_DIRS" >> $HOME/.bashrc
  echo "export LIBRARY_PATH=\$LIBRARY_PATH$LIBRARY_DIRS" >> $HOME/.bashrc
  printf "\nAll libraries set. Run a new shell or execute 'source ~/.bashrc'.\n"
else
  echo "ERROR: no libraries provided"
fi
