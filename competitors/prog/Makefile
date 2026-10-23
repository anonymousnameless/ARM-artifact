CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17
PREFIX ?= /usr/local
CPLEX_DIR ?= /opt/ibm/ILOG/CPLEX_Studio2211
LP_SOLVE_LIB ?= -llpsolve55
CPPFLAGS += -Isrc/headers -I$(PREFIX)/include
LDFLAGS += -L$(PREFIX)/lib
LDLIBS += $(LP_SOLVE_LIB) -lqhull_r -ldl

SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(patsubst src/%.cpp,build/%.o,$(SOURCES))

.PHONY: all clean
all: bin/prog
bin/prog: $(OBJECTS)
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)
build/%.o: src/%.cpp
	mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@
-include $(OBJECTS:.o=.d)
clean:
	rm -rf build bin
