CXX := g++
CXXFLAGS := -O -Wall
INC := $(shell root-config --cflags)
LIB := $(shell root-config --libs)

default: compila

compila:
	$(CXX) -c Elements.C -o Elements.o
	$(CXX) $(INC) -c Particle.C -o Particle.o
	$(CXX) $(INC) -c Tracker.C -o Tracker.o
	$(CXX) $(INC) -c Calorimeter.C -o Calorimeter.o
	$(CXX) $(INC) -c programma.cpp -o programma.o
	$(CXX) Elements.o Particle.o Tracker.o Calorimeter.o programma.o $(LIB) -o program

clean:
	rm -f *.o
	rm -f program
	rm -f *.pdf
	rm -f *~
