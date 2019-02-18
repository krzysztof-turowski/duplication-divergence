CC = g++
CC_FLAGS = -O3 -std=c++17 -Wall -Wextra -Wstrict-aliasing -Wpedantic -fmax-errors=5 -Werror -Wunreachable-code -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Wno-unused -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option

PARALLEL_FLAG = -pthread
FLAGS = $(PARALLEL_FLAG)

NAUTY_LIB = lib/nauty/nauty.a

SRCS := $(wildcard *.cpp)
EXEC := $(patsubst %.cpp,%,$(SRCS))

all: $(EXEC)

%: %.cpp
	@$(CC) $(CC_FLAGS) $(FLAGS) $< -o $@

dd_automorphisms: dd_automorphisms.cpp
	@$(CC) $(CC_FLAGS) $(FLAGS) $< $(NAUTY_LIB) -o $@
