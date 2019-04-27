CXX = g++
EXT = cpp
DEBUGFLAGS = -g -ggdb -O0
CXXFLAGS = -Wall -std=c++11 -pedantic $(DEBUGFLAGS)
TARGET = test

all: $(TARGET)

$(TARGET): $(TARGET).o   #link
	$(CXX) $(CXXFLAGS) $(TARGET).o -o $(TARGET)

$(TARGET).o: $(TARGET).cpp quadtree.hpp # compile
	$(CXX) -c $(CXXFLAGS) $(TARGET).$(EXT)

run: all
	./$(TARGET) 

valgrind: all
	valgrind --leak-check=full --track-origins=yes ./$(TARGET)

clean:
	rm -f *.o $(TARGET) *.dat



