UNAME := $(shell uname)
CC = g++
COMPILER_FLAGS =-fmax-errors=5 -Wlogical-op -Wstrict-null-sentinel -Wnoexcept
FLAGS = -std=c++17 -lstdc++ -Wall -Wextra -Wstrict-aliasing -Wpedantic -Werror -Wunreachable-code -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-overflow=2 -Wswitch-default -Wundef -fdiagnostics-show-option -O3 -pthread

ifeq ($(UNAME),Darwin)
  COMPILER_FLAGS += -Xclang -fopenmp -lomp
else
  COMPILER_FLAGS += -fopenmp
endif

INC = -isystem lib/nauty -Llib/nauty
GRAPH_LIB = snap
NAUTY_LIB = -l:nauty.a -Wno-unused-variable -DTHREADS=1
LP_SOLVER = glpk

INPUT_FLAGS = -lsnap -Wno-error -fpermissive

GRAPH_FLAGS = -D$(GRAPH_LIB) -DNDEBUG
ifeq ($(GRAPH_LIB),snap)
	GRAPH_FLAGS += -lsnap -Wno-error -fpermissive
	INC += -isystem lib/snap/snap-core -isystem lib/snap/glib-core -Llib/snap/snap-core -Llib/snap/glib-core 
else ifeq ($(GRAPH_LIB),networkit)
	GRAPH_FLAGS += -lnetworkit
	INC += -isystem lib/networkit/networkit/cpp -Llib/networkit/build_lib
endif

MATH_FLAGS = -lgmp -lgmpxx

LP_FLAGS = -D$(LP_SOLVER) -DNDEBUG
ifeq ($(LP_SOLVER),glpk)
	LP_FLAGS += -lglpk
else ifeq ($(LP_SOLVER),gurobi)
	LP_FLAGS += -lgurobi_c++ -lgurobi81 -Wno-error
else
endif

SRC_DIR = src
BUILD_DIR = .

CPP_HDRS := $(wildcard $(SRC_DIR)/*.h $(SRC_DIR)/*/*.h)
CPP_SRCS := $(wildcard $(SRC_DIR)/*.cpp $(SRC_DIR)/*/*.cpp)
OBJ_FILES := $(patsubst %.cpp,%.o,$(wildcard $(SRC_DIR)/*/*.cpp))
PY_SRCS := $(wildcard $(SRC_DIR)/*.py)
EXEC := $(patsubst $(SRC_DIR)/%.cpp,%,$(wildcard $(SRC_DIR)/*.cpp))

all: g++ check

g++: $(EXEC)

clang++: CC = clang++
clang++: COMPILER_FLAGS = -ferror-limit=5 -Wno-unused-const-variable -fopenmp=libomp
clang++: $(EXEC)

debug: COMPILER_FLAGS += -ggdb -fsanitize=thread,undefined
debug: $(EXEC)

%: $(SRC_DIR)/%.cpp $(CPP_HDRS) $(OBJ_FILES)
	@$(CC) $(INC) $(FLAGS) $(COMPILER_FLAGS) $< $(OBJ_FILES) $(INPUT_FLAGS) $(NAUTY_LIB) $(GRAPH_FLAGS) -o $(BUILD_DIR)/$@

%.o: %.cpp
	@$(CC) -c $(INC) $(FLAGS) $(COMPILER_FLAGS) $< $(INPUT_FLAGS) $(GRAPH_FLAGS) -o $(BUILD_DIR)/$@

dd_automorphisms: $(SRC_DIR)/dd_automorphisms.cpp $(CPP_HDRS) $(OBJ_FILES)
	@$(CC) $(INC) $(FLAGS) $(COMPILER_FLAGS) $< $(OBJ_FILES) $(INPUT_FLAGS) $(NAUTY_LIB) -o $(BUILD_DIR)/$@

dd_temporal_%: $(SRC_DIR)/dd_temporal_%.cpp $(CPP_HDRS) $(OBJ_FILES)
	@$(CC) $(INC) $(FLAGS) $(COMPILER_FLAGS) $< $(OBJ_FILES) $(INPUT_FLAGS) $(GRAPH_FLAGS) $(LP_FLAGS) $(MATH_FLAGS) -o $(BUILD_DIR)/$@

check:
	cpplint --linelength=100 --extensions=cpp,h --filter=-legal/copyright,-build/c++11,-build/include,-build/namespaces,-runtime/references,-runtime/string $(CPP_SRCS) $(CPP_HDRS)
	pylint --disable=bad-whitespace,invalid-name,missing-docstring,too-many-locals,star-args,no-member,fixme,superfluous-parens --max-line-length=100 --extension-pkg-whitelist=numpy $(PY_SRCS)

clean:
	@rm -vf $(addprefix $(BUILD_DIR)/,$(EXEC))

.PHONY: all check clean
