
CXX = g++
CXXFLAGS = -Wall -g -lstdc++

a.out: Main.o DataStorage.o SetVariables.o Solutions.o Utils.o LinearSolvers.o
	$(CXX) $(CXXFLAGS) -o a.out Main.o DataStorage.o SetVariables.o Solutions.o Utils.o LinearSolvers.o

Main.o: Main.cpp DataStorage.h SetVariables.h Solutions.h Utils.h LinearSolvers.h
	$(CXX) $(CXXFLAGS) -c Main.cpp

DataStorage.o: DataStorage.h

SetVariables.o: SetVariables.h

Solutions.o: Solutions.h DataStorage.h SetVariables.h LinearSolvers.h

Utils.o: Utils.h

LinearSolvers.o: LinearSolvers.h Utils.h

clean:
	$(RM) count *.o *~|
