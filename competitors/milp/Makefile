CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17
PREFIX ?= /usr/local
CPLEX_DIR ?= /opt/ibm/ILOG/CPLEX_Studio2211
CPPFLAGS += -Isrc/headers -I$(PREFIX)/include
LDFLAGS += -L$(PREFIX)/lib
LDLIBS += -lqhull_r -ldl
CPPFLAGS += -DIL_STD -I$(CPLEX_DIR)/cplex/include -I$(CPLEX_DIR)/concert/include
LDFLAGS += -L$(CPLEX_DIR)/cplex/lib/x86-64_linux/static_pic -L$(CPLEX_DIR)/concert/lib/x86-64_linux/static_pic
LDLIBS += -lilocplex -lcplex -lconcert -lpthread
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(patsubst src/%.cpp,build/%.o,$(SOURCES))

.PHONY: all clean
all: bin/milp
bin/milp: $(OBJECTS)
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)
build/%.o: src/%.cpp
	mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@
-include $(OBJECTS:.o=.d)
clean:
	rm -rf build bin
