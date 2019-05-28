CXX = g++
EXT = cpp

DEBUGFLAGS = -g -ggdb -O3 #-fno-inline
OFLAGS = -O3 #-fopenmp
#OFLAGS = -O3 -ffast-math #-fopenmp
LIBS = -lgsl -lcblas
CXXFLAGS = -Wall  -std=c++17 -pedantic  

TARGET = test2
DEPENDS = makefile $(TARGET).$(EXT) abstract_orthtree.hpp balanced_orthtree.hpp \
		  point_orthtree.hpp vector.hpp fmm_tree.hpp series_expansion.hpp \
		  multipole_expansion.hpp local_expansion.hpp debugging.hpp

all: $(TARGET)

# $(TARGET): $(TARGET).o $(DEPENDS)  #link
# 	$(CXX) $(CXXFLAGS) $(TARGET).o -o $(TARGET)
# 
# $(TARGET).o: $(TARGET).cpp $(DEPENDS) # compile
# 	$(CXX) -c $(CXXFLAGS) $(TARGET).$(EXT)

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

clean:
	rm -f *.o $(TARGET) *.dat *.out



