CC = g++
FLAGS = -std=c++17 -lstdc++ -Wall -Wextra -Wstrict-aliasing -Wpedantic -Werror -Wunreachable-code -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -fdiagnostics-show-option -O3 -pthread

NAUTY_LIB = lib/nauty/nauty.a
LP_FLAGS = -lgmp -lgmpxx -lglpk

CPP_SRCS := $(wildcard *.cpp)
PY_SRCS := $(wildcard *.py)
EXEC := $(patsubst %.cpp,%,$(CPP_SRCS))

all: g++ check

g++: COMPILER_FLAGS =-fmax-errors=5
g++: $(EXEC)

clang++: CC = clang++
clang++: COMPILER_FLAGS = -ferror-limit=5 -Wno-unknown-warning-option -Wno-c++11-extensions -Wno-unused-const-variable
clang++: $(EXEC)

%: %.cpp
	@$(CC) $(FLAGS) $(COMPILER_FLAGS) $< -o $@

dd_automorphisms: dd_automorphisms.cpp
	@$(CC) $(FLAGS) $(COMPILER_FLAGS) $< $(NAUTY_LIB) -o $@

dd_temporal_order: dd_temporal_order.cpp
	@$(CC) $(FLAGS) $(COMPILER_FLAGS) $< $(NAUTY_LIB) $(LP_FLAGS) -o $@

check:
	cppcheck --enable=all --force --suppress=*:lib/* $(CPP_SRCS)
	pylint --disable=bad-whitespace,invalid-name,missing-docstring,too-many-locals,star-args,no-member,fixme --max-line-length=100 --extension-pkg-whitelist=numpy $(PY_SRCS)

clean:
	@rm -vf $(EXEC)

.PHONY: all check clean
