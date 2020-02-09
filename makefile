
CXX = mpicxx
CXXFLAGS = -Wall -g -lstdc++

run: ParallelMain.o Parameters.o DataStorage.o ParallelSolutions.o ParallelSolvers.o
	$(CXX) -o run ParallelMain.o Parameters.o DataStorage.o ParallelSolutions.o ParallelSolvers.o

Main.o: Main.cpp Parameters.h DataStorage.h ParallelSolutions.h ParallelSolvers.h
	$(CXX) -c ParallelMain.cpp

Parameters.o: Parameters.h

DataStorage.o: DataStorage.h

ParallelSolutions.o: ParallelSolutions.h Parameters.h DataStorage.h

ParallelSolvers.o: ParallelSolvers.h

clean:
	$(RM) count *.o *~|
