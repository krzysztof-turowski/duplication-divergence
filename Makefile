CC = g++
COMPILER_FLAGS =-fmax-errors=5 -Wlogical-op -Wstrict-null-sentinel -Wnoexcept -fopenmp
FLAGS = -std=c++17 -lstdc++ -Wall -Wextra -Wstrict-aliasing -Wpedantic -Werror -Wunreachable-code -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-overflow=2 -Wswitch-default -Wundef -fdiagnostics-show-option -O3 -pthread
GRAPH_LIB = koala
NAUTY_LIB = lib/nauty/nauty.a -Wno-unused-variable -DTHREADS=1
LP_SOLVER = glpk

GRAPH_FLAGS = -D$(GRAPH_LIB) -DNDEBUG
ifeq ($(GRAPH_LIB),snap)
	GRAPH_FLAGS += -lsnap -Wno-error -fpermissive
else ifeq ($(GRAPH_LIB),networkit)
	GRAPH_FLAGS += -lnetworkit
endif

LP_FLAGS = -D$(LP_SOLVER) -DNDEBUG
ifeq ($(LP_SOLVER),glpk)
	LP_FLAGS += -lglpk -lgmp -lgmpxx
else ifeq ($(LP_SOLVER),gurobi)
	LP_FLAGS += -lgurobi_c++ -lgurobi81 -lgmp -lgmpxx -Wno-error
else
endif

CPP_HDRS := $(wildcard *.h)
CPP_SRCS := $(wildcard *.cpp)
PY_SRCS := $(wildcard *.py)
EXEC := $(patsubst %.cpp,%,$(CPP_SRCS))

all: g++ check

g++: $(EXEC)

clang++: CC = clang++
clang++: COMPILER_FLAGS = -ferror-limit=5 -Wno-unused-const-variable
clang++: $(EXEC)

debug: COMPILER_FLAGS += -ggdb -fsanitize=thread,undefined
debug: $(EXEC)

%: %.cpp
	@$(CC) $(FLAGS) $(COMPILER_FLAGS) $< $(GRAPH_FLAGS) -o $@

dd_automorphisms: dd_automorphisms.cpp
	@$(CC) $(FLAGS) $(COMPILER_FLAGS) $< $(NAUTY_LIB) -o $@

dd_temporal_bound: dd_temporal_bound.cpp
	@$(CC) $(FLAGS) $(COMPILER_FLAGS) $< $(GRAPH_FLAGS) $(LP_FLAGS) -o $@

check:
	cpplint --linelength=100 --extensions=cpp,h --filter=-legal/copyright,-build/c++11,-build/namespaces,-runtime/references,-runtime/string $(CPP_SRCS) $(CPP_HDRS)
	pylint --disable=bad-whitespace,invalid-name,missing-docstring,too-many-locals,star-args,no-member,fixme,superfluous-parens --max-line-length=100 --extension-pkg-whitelist=numpy $(PY_SRCS)

clean:
	@rm -vf $(EXEC)

.PHONY: all check clean
