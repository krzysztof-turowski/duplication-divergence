#!/bin/bash

# Script for libraries configuration. Run: bash ./configure_lish.sh

set -e
if [ "$#" -gt 0 ]; then
  CURRENT_DIR=`pwd`
  INCLUDE_DIRS=""
  LIBRARY_DIRS=""
  for lib in "$@"
  do
    if [[ "$lib" == "glpk" ]]; then
      DIR=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep glpk | head -1)
      ld -L ${DIR:-.} -lglpk 2>/dev/null || CHECK=`echo $?`
      if [ ${CHECK:-0} -eq 0 ]; then
        echo "GLPK already present"
      else
        echo "GLPK not present"
        if [[ "$lib" == "glpk" ]]; then
          GLPK_LIB=$(find -L /usr/lib /usr/local/lib $HOME -name libglpk.a -print0 2>/dev/null | head -1)
          GLPK_DIR="$(dirname "$(dirname "$GLPK_LIB")")"
          INCLUDE_DIRS="$INCLUDE_DIRS:$GLPK_DIR/include"
          LIBRARY_DIRS="$LIBRARY_DIRS:$GLPK_DIR/lib"
        fi
      fi
    fi
    if [[ "$lib" == "gurobi" ]]; then
      DIR=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep gurobi | head -1)
      ld -L ${DIR:-.} -lgurobi_c++ -lgurobi81 2>/dev/null || CHECK=`echo $?`
      if [ ${CHECK:-0} -eq 0 ]; then
        echo "Gurobi already present"
      else
        echo "Gurobi not present"
        echo $(uname) >> MAC_PATH
        if [[ $UNAME == "Darwin" ]]; then
          MAC_PATH = "/Library"
        fi
        if [[ -z "$GUROBI_HOME" ]]; then
          GUROBI_LIB=$(find -L /usr/lib $HOME $CURRENT_DIR $MAC_PATH -path */lib/libgurobi_c++.a -print0 2>/dev/null | head -1)
          GUROBI_DIR="$(dirname "$(dirname "$GUROBI_LIB")")"
          echo "export GUROBI_HOME=$GUROBI_DIR" >> $HOME/.bashrc
        else
          GUROBI_DIR=$GUROBI_HOME
        fi
        if [[ $PATH != *"$GUROBI_DIR"* ]]; then
          echo "export PATH=\$PATH:$GUROBI_DIR/bin" >> $HOME/.bashrc
        fi
        if [[ -z "$GRB_LICENSE_FILE" ]]; then
          GUROBI_LICENSE=$(find -L $HOME $CURRENT_DIR -name gurobi.lic -print0 2>/dev/null | head -1)
          echo "export GRB_LICENSE_FILE=$GUROBI_LICENSE" >> $HOME/.bashrc
        fi
        cd $GUROBI_DIR/src/build
        make clean
        make c++
        ln -sf $GUROBI_DIR/src/build/libgurobi_c++.a $GUROBI_DIR/lib/libgurobi_c++.a
        INCLUDE_DIRS="$INCLUDE_DIRS:$GUROBI_DIR/include"
        LIBRARY_DIRS="$LIBRARY_DIRS:$GUROBI_DIR/lib"
      fi
    fi
    if [[ "$lib" == "koala" ]]; then
      INCLUDE_DIRS="$INCLUDE_DIRS:$CURRENT_DIR/lib/koala"
    fi
    if [[ "$lib" == "nauty" ]]; then
      DIR=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep nauty | head -1)
      ld -L ${DIR:-.} -l:nauty.a 2>/dev/null || CHECK=`echo $?`
      if [ ${CHECK:-0} -eq 0 ]; then
        echo "nauty already present"
      else
        echo "nauty not present"
        cd $CURRENT_DIR/lib/nauty
        chmod +x ./configure
        ./configure
        make
        INCLUDE_DIRS="$INCLUDE_DIRS:$CURRENT_DIR/lib/nauty"
        LIBRARY_DIRS="$LIBRARY_DIRS:$CURRENT_DIR/lib/nauty"
      fi
    fi
    if [[ "$lib" == "networkit" ]]; then
      DIR=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep networkit | head -1)
      ld -L ${DIR:-.} -lnetworkit 2>/dev/null || CHECK=`echo $?`
      if [ ${CHECK:-0} -eq 0 ]; then
        echo "NetworKit already present"
      else
        echo "NetworKit not present"
        mkdir $CURRENT_DIR/lib/networkit/build_lib
        cd $CURRENT_DIR/lib/networkit/build_lib
        cmake ..
        make
        INCLUDE_DIRS="$INCLUDE_DIRS:$CURRENT_DIR/lib/networkit/networkit/cpp"
        LIBRARY_DIRS="$LIBRARY_DIRS:$CURRENT_DIR/lib/networkit/build_lib"
      fi
    fi
    if [[ "$lib" == "snap" ]]; then
      DIR=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep snap | head -1)
      ld -L ${DIR:-.} -lsnap 2>/dev/null || CHECK=`echo $?`
      if [ ${CHECK:-0} -eq 0 ]; then
        echo "SNAP already present"
      else
        echo "SNAP not present"
        cd $CURRENT_DIR/lib/snap/snap-core
        make lib
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
