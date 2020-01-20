
CXX = mpic++
CXXFLAGS = -Wall -g -lstdc++

run: ParallelMain.o Parameters.o DataStorage.o ParallelSolutions.o
	$(CXX) -o run ParallelMain.o Parameters.o DataStorage.o ParallelSolutions.o

Main.o: Main.cpp Parameters.h DataStorage.h ParallelSolutions.h
	$(CXX) -c ParallelMain.cpp

Parameters.o: Parameters.h

DataStorage.o: DataStorage.h

ParallelSolutions.o: ParallelSolutions.h Parameters.h DataStorage.h

clean:
	$(RM) count *.o *~|
