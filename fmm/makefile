CXX = g++-8
EXT = cpp

DEBUGFLAGS = -g -ggdb -O3 -ffast-math #-fno-inline
OFLAGS = -O3 -ffast-math -fopenmp
LIBS = -lstdc++fs -lgsl -lgslcblas 
CXXFLAGS = -Wall -Wno-unknown-pragmas -fno-builtin -std=c++17 -pedantic 

TARGET = examples/minimal
#TARGET = test/fmm_tree
DEPENDS = makefile $(TARGET).$(EXT) abstract_orthtree.hpp balanced_orthtree.hpp \
		   vector.hpp abstract_fmm_tree.hpp balanced_fmm_tree.hpp series_expansion.hpp \
		  multipole_expansion.hpp local_expansion.hpp \
		  adaptive_fmm_tree.hpp fmm_general.hpp fields.hpp

all: $(TARGET)

$(TARGET): $(DEPENDS)  #link
	$(CXX) $(CXXFLAGS) $(OFLAGS) $(TARGET).$(EXT) -o $(TARGET) $(LIBS) 

debug: 
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) $(TARGET).$(EXT) -o $(TARGET) $(LIBS)

run: all
	./$(TARGET) 

valgrind: debug
	valgrind --leak-check=full --track-origins=yes ./$(TARGET) 

callgrind: debug
	valgrind --tool=callgrind --dump-instr=yes --collect-jumps=yes \
		--collect-atstart=no --callgrind-out-file=./callgrind.out ./$(TARGET) 
	kcachegrind callgrind.out & 

perf: all
	perf stat -r 10 ./$(TARGET)

clean:
	rm -f *.o $(TARGET) logs/*.dat *.out *.dat




