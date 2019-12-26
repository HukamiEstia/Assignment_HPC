#include <chrono>

#include "DataStorage.h"
#include "SetVariables.h"
#include "Solutions.h"
#include "Utils.h"
#include "LinearSolvers.h"

void printDataForLatex(vector<vector<double>> matrix, double dT) {

	cout << "0 & 0.1 & 0.2 & 0.3 & 0.4 & 0.5" << endl;
	for (int i = 0; i < matrix[0].size(); i++) {
		cout << matrix[0][i];
		for (int j = 1; j < 6; j++) {
			cout << " & " << matrix[j * (0.1 / dT)][i] ;
		}
		cout << " \\\\" << endl;
	}
}

int main() {
	//Definition of the problem's parameters
	SetVariables init(31, 93, 0.5, 149, 149, 38);

	//Initialisation of the container object to hold the data of the computed results
	DataStorage analyt(0.05, 0.01), DFF(0.05, 0.01), Rich(0.05, 0.01), Laas(0.05, 0.01), CN(0.05, 0.01);

	//Call each function to solve the problem using the studied schemes
	cout << "Computing Analytical solution" << endl;
	analytical(analyt, init);
	cout << "Computing Richardson solution" << endl;
	richardson(Rich, init);
	cout << "Computing duFort Frankel solution" << endl;
	duFortFrankel(DFF, init);
	cout << "Computing Laasonen Implicit solution" << endl;
	laasonenImplicit(Laas, init, Thomas);
	crankNicholson(CN, init, Thomas);

	//Saving the results
	cout << "Saving results." << endl;
	analyt.saveDataToFile("Analytical");
	Rich.saveDataToFile("Richardson");
	DFF.saveDataToFile("DuFort-Frankel");
	Laas.saveDataToFile("Laasonen");
	CN.saveDataToFile("Crank-Nicholson");

	//updating the discretization steps
	Laas.setDX(0.05);
	Laas.setDT(0.01);

	//Then we compute The solution using the Laasonen method for different time resolutions
	cout << "Computing Laasonen at dX=0.05, dT=0.01" << endl;
	for (int i = 0; i < 10; i++) {
		laasonenImplicit(Laas, init, Thomas);
	}
	laasonenImplicit(Laas, init, Thomas);
	Laas.saveSpecifiedTimeDataToFile("Laasonen0.01");

	cout << "Computing Laasonen at dX=0.05, dT=0.025" << endl;
	Laas.setDT(0.025);
	for (int i = 0; i < 10; i++) {
		laasonenImplicit(Laas, init, Thomas);
	}
	laasonenImplicit(Laas, init, Thomas);
	Laas.saveSpecifiedTimeDataToFile("Laasonen0.025");

	cout << "Computing Laasonen at dX=0.05, dT=0.05" << endl;
	Laas.setDT(0.05);
	for (int i = 0; i < 10; i++) {
		laasonenImplicit(Laas, init, Thomas);
	}
	laasonenImplicit(Laas, init, Thomas);
	Laas.saveSpecifiedTimeDataToFile("Laasonen0.05");

	cout << "Computing Laasonen at dX=0.05, dT=0.1" << endl;
	Laas.setDT(0.1);
	for (int i = 0; i < 10; i++) {
		laasonenImplicit(Laas, init, Thomas);
	}
	laasonenImplicit(Laas, init, Thomas);
	Laas.saveSpecifiedTimeDataToFile("Laasonen0.1");
	/**/
}
