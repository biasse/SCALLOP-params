
.PHONY: all
all: main

SOURCES = main.cpp Fp2k.cpp ec.cpp ecp.cpp isog.cpp utils.cpp id2iso.cpp quaternions.cpp
HEADERS = Fp2k.hpp ec.hpp ecp.hpp isog.hpp utils.hpp id2iso.hpp quaternions.hpp

CXXFLAGS = -std=c++17 -pedantic -Wall -Wextra -O2
LDFLAGS = -lm -lntl -lgmp -lfplll -lmpfr

objs-main/%.o: CXXFLAGS += -march=native -DNDEBUG
objs-debug/%.o: CXXFLAGS += -g

objs-main/%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

objs-debug/%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

main: $(patsubst %.cpp,objs-main/%.o,$(SOURCES))
	$(CXX) $(LDFLAGS) $^ -o $@

debug: $(patsubst %.cpp,objs-debug/%.o,$(SOURCES))
	$(CXX) $(LDFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -f main objs-main/*
	rm -f debug objs-debug/*

