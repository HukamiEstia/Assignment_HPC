#include <chrono>

#include "DataStorage.h"
#include "Parameters.h"
#include "Solutions.h"
#include "Utils.h"
#include "LinearSolvers.h"

int main() {
	//Definition of the problem's parameters
	Parameters params(31, 93, 0.5, 149, 38);

	//Initialisation of the container object to hold the data of the computed results
	DataStorage analyt(0.05, 0.01), ftcs(0.05, 0.01), Laas(0.05, 0.01), CN(0.05, 0.01);

	//Call each function to solve the problem using the studied schemes
	cout << "Computing Analytical solution" << endl;
	analytical(analyt, params);
	cout << "Computing Forward Time/Central Space solution" << endl;
	FTCS(ftcs, params);
	cout << "Computing Laasonen Implicit solution" << endl;
	laasonenImplicit(Laas, params, Thomas);
	cout << "Computing Crank-Nicholson Implicit solution" << endl;
	crankNicholson(CN, params, Thomas);
	/*
	//Saving the results
	cout << "Saving results." << endl;
	analyt.saveDataToFile("Analytical");
	ftcs.saveDataToFile("FTCS");
	Laas.saveDataToFile("Laasonen");
	CN.saveDataToFile("Crank-Nicholson");
	/**/
}
