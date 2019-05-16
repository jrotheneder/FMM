CXX = g++
EXT = cpp
DEBUGFLAGS = -g -ggdb -O0
CXXFLAGS = -Wall -std=c++17 -pedantic $(DEBUGFLAGS)
TARGET = expansions_test
DEPENDS = abstract_orthtree.hpp balanced_orthtree.hpp point_orthtree.hpp vector.hpp fmm_tree.hpp series_expansion.hpp multipole_expansion.hpp local_expansion.hpp debugging.hpp

all: $(TARGET)

$(TARGET): $(TARGET).o $(DEPENDS)  #link
	$(CXX) $(CXXFLAGS) $(TARGET).o -o $(TARGET)

$(TARGET).o: $(TARGET).cpp $(DEPENDS) # compile
	$(CXX) -c $(CXXFLAGS) $(TARGET).$(EXT)

run: all
	./$(TARGET) 

valgrind: all
	valgrind --leak-check=full --track-origins=yes ./$(TARGET)

callgrind: 
	valgrind --tool=callgrind --dump-instr=yes --simulate-cache=yes --collect-jumps=yes \
		--callgrind-out-file=./callgrind.out ./$(TARGET) 

clean:
	rm -f *.o $(TARGET) *.dat *.out



