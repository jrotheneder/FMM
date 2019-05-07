CXX = g++
EXT = cpp
DEBUGFLAGS = -g -ggdb -O0
CXXFLAGS = -Wall -std=c++11 -pedantic $(DEBUGFLAGS)
TARGET = test2
DEPENDS = abstract_orthtree.hpp balanced_orthtree.hpp point_orthtree.hpp vector.hpp fmm_tree.hpp multipole_expansion.hpp local_expansion.hpp

all: $(TARGET)

$(TARGET): $(TARGET).o $(DEPENDS)  #link
	$(CXX) $(CXXFLAGS) $(TARGET).o -o $(TARGET)

$(TARGET).o: $(TARGET).cpp $(DEPENDS) # compile
	$(CXX) -c $(CXXFLAGS) $(TARGET).$(EXT)

run: all
	./$(TARGET) 

valgrind: all
	valgrind --leak-check=full --track-origins=yes ./$(TARGET)

clean:
	rm -f *.o $(TARGET) *.dat



