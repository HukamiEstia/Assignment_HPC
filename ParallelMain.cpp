#include <chrono>

#include "DataStorage.h"
#include "Parameters.h"
#include "ParallelSolutions.h"
#include "Utils.h"
#include "LinearSolvers.h"

int main() {
	//Definition of the problem's parameters
	Parameters params(31, 93, 0.5, 149, 38);
	
	//Initialisation of the container object to hold the data of the computed results
	DataStorage analyt(0.5, 0.001);
	DataStorage ftcs(0.5, 0.001);

	
	//Call each function to solve the problem using the studied schemes
	//cout << "Computing Analytical solution" << endl;
	//ParallelAnalytical(analyt, params);

	cout << "Computing Forward Time/Central Space solution" << endl;
	FTCS(ftcs, params);
	
	cout << "Saving results." << endl;
	//analyt.saveDataToFile("Analytical");
	ftcs.saveDataToFile("FTCS");

	/*
	cout << "Computing Laasonen Implicit solution" << endl;
	laasonenImplicit(Laas, params, Thomas);
	cout << "Computing Crank-Nicholson Implicit solution" << endl;
	crankNicholson(CN, params, Thomas);
	//Saving the results
	
	Laas.saveDataToFile("Laasonen");
	CN.saveDataToFile("Crank-Nicholson");
	/**/
}
