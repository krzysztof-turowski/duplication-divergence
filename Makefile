CC = g++
COMPILER_FLAGS =-fmax-errors=5 -Wlogical-op -Wstrict-null-sentinel -Wnoexcept
FLAGS = -std=c++17 -lstdc++ -Wall -Wextra -Wstrict-aliasing -Wpedantic -Werror -Wunreachable-code -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-overflow=2 -Wswitch-default -Wundef -fdiagnostics-show-option -O3 -pthread
NAUTY_LIB = lib/nauty/nauty.a -Wno-unused-variable
LP_SOLVER = gurobi

ifeq ($(LP_SOLVER),glpk)
	LP_FLAGS = -Dglpk -lglpk -lgmp -lgmpxx
else
	LP_FLAGS = -Dgurobi -Ilib/gurobi/include -Llib/gurobi/lib -lgurobi_c++ -lgurobi81 -lgmp -lgmpxx
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

%: %.cpp
	@$(CC) $(FLAGS) $(COMPILER_FLAGS) $< -o $@

dd_automorphisms: dd_automorphisms.cpp
	@$(CC) $(FLAGS) $(COMPILER_FLAGS) $< $(NAUTY_LIB) -o $@

dd_temporal_order: dd_temporal_order.cpp
	@$(CC) $(FLAGS) $(COMPILER_FLAGS) $< $(LP_FLAGS) -o $@

check:
	cpplint --linelength=100 --extensions=cpp,h --filter=-legal/copyright,-build/c++11,-build/namespaces,-runtime/references,-runtime/string $(CPP_SRCS) $(CPP_HDRS)
	pylint --disable=bad-whitespace,invalid-name,missing-docstring,too-many-locals,star-args,no-member,fixme --max-line-length=100 --extension-pkg-whitelist=numpy $(PY_SRCS)

clean:
	@rm -vf $(EXEC)

.PHONY: all check clean
