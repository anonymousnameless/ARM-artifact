CXX = g++
CPLEX_DIR ?= /opt/ibm/ILOG/CPLEX_Studio2211
PREFIX ?= /usr/local
CXXFLAGS ?= -O2 -std=c++17 -fopenmp
CPPFLAGS += -DIL_STD -Isrc -Isrc/data -Isrc/algorithm -Isrc/geometry -Isrc/spatial -Isrc/utils -I$(PREFIX)/include -I$(CPLEX_DIR)/cplex/include -I$(CPLEX_DIR)/concert/include
LDFLAGS += -fopenmp -L$(PREFIX)/lib -L$(CPLEX_DIR)/cplex/lib/x86-64_linux/static_pic -L$(CPLEX_DIR)/concert/lib/x86-64_linux/static_pic
# Supply the full path to a locally built lp_solve archive if needed.
LP_SOLVE_LIB ?= -llpsolve55
LDLIBS = $(LP_SOLVE_LIB) -lqhull_r -lortools -lprotobuf -labsl_base -labsl_strings -labsl_synchronization -labsl_time -labsl_int128 -labsl_log_internal_check_op -labsl_log_internal_message -lgmpxx -lgmp -lilocplex -lcplex -lconcert -lpthread -ldl
SOURCES = $(wildcard src/*.cpp src/algorithm/*.cpp src/data/*.cpp src/geometry/*.cpp src/spatial/*.cpp src/utils/*.cpp)
OBJECTS = $(patsubst src/%.cpp,build/%.o,$(SOURCES))

.PHONY: all clean
all: bin/arf
bin/arf: $(OBJECTS)
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)
build/%.o: src/%.cpp
	mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@
-include $(OBJECTS:.o=.d)
clean:
	rm -rf build bin
