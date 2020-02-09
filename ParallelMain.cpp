#include <chrono>

#include "DataStorage.h"
#include "Parameters.h"
#include "ParallelSolutions.h"
#include "Utils.h"
#include "LinearSolvers.h"
#include "fstream"
#include <string.h>

int main(int argc, char* argv[]) {
	double Time;
	ofstream times;
	string space, time, name, time_name;

	space = argv[1];
	time = argv[2];
	name = argv[3];
	name = name + "_" + space + "_" + time; 
	time_name = "time_" + name + ".dat";
	times.open(time_name);
	//Definition of the problem's parameters
	Parameters params(31, 93, 5, 149, 38);
	
	//Initialisation of the container object to hold the data of the computed results
	DataStorage data(atof(argv[1]), atof(argv[2]));


	// Call each function to solve the problem using the studied schemes
	if (strcmp(argv[3], "analytical") == 0){
		cout << "Computing Analytical solution" << endl;
		Time = ParallelAnalytical(data, params);
		times << "analytical solving " << argv[1] << ";" << argv[2] << ": " << Time << "\n";
		cout << "Saving results." << endl;
		data.saveDataToFile(name);
	}
	else if (strcmp(argv[3], "ftcs") == 0) {
		cout << "Computing Forward Time/Central Space solution" << endl;
		Time = FTCS(data, params);
		times << "ftcs solving " << argv[1] << ";" << argv[2] << ": " << Time << "\n";
		cout << "Saving results." << endl;
		data.saveDataToFile(name);
	}
	else if (strcmp(argv[3], "laasonen")  == 0) {
		cout << "Computing Laasonen Implicit solution" << endl;
		Time = laasonenImplicit(data, params);
		times << "laasonen solving " << argv[1] << ";" << argv[2] << ": " << Time << "\n";
		cout << "Saving results." << endl;
		data.saveDataToFile(name);
	}
	else if (strcmp(argv[3], "crank") == 0){
		cout << "Computing Crank-Nicholson Implicit solution" << endl;
		Time = crankNicholson(data, params);
		times << "crank-nicholson solving " << argv[1] << ";" << argv[2] << ": " << Time << "\n";
		cout << "Saving results." << endl;
		data.saveDataToFile(name);
	}
	/**/
}
