#include "Solutions.h"

void debugSimple(DataStorage &storage, const SetVariables variables) {
	/*
	This method is forward time forward space
	used to debug during the first stages of development
	/**/
	double left, right, central;
	double dT = sqrt(pow(storage.getDX(), 2) / (variables.getD()*2));
	unsigned int numberOfT = variables.getT() / storage.getDT() + 1;
	unsigned int numberOfX = variables.getX() / storage.getDX() + 1;
	vector<vector<double>> finalVector(numberOfT);
	vector<double> previousTime(numberOfX);
	finalVector[0] = variables.getT0(numberOfX);

	//define r = D.(dt/dx^2)
	double r = (variables.getD() * dT) / (pow(storage.getDX(), 2));

	for (unsigned int time = 1; time < numberOfT; time++) {
		previousTime = finalVector[time - 1];
		vector<double> construct(numberOfX);
		for (unsigned int space = 0; space < numberOfX; space++) {
			if (space == 0) {
				construct[0] = variables.getX0();
			}else if (space == numberOfX - 1) {
				construct[space] = variables.getXMax();
			}else {
				////T(i-1,n)
				left = previousTime[space - 1];
				//T(i,n)
				central = previousTime[space];
				//T(i+1,n)
				right = previousTime[space + 1];
				//compute T(i,n+1)
				construct[space] = central+r*(left-2*central+right);
			}
		}
		finalVector[time] = construct;
	}
	storage.setData(finalVector);
}

void duFortFrankel(DataStorage &storage, const SetVariables variables){
	/*
	Solve using the DuFort-Frankel scheme 
	/**/
	double right, left, down;
	unsigned int numberOfT = (variables.getT() / storage.getDT())+1;	
	unsigned int numberOfX = (variables.getX() / storage.getDX())+1;
	vector<vector<double>> finalVector(numberOfT);
	vector<double> previousTime(numberOfX);
	finalVector[0] = variables.getT0(numberOfX);
	double r = 2*(variables.getD() * storage.getDT()) / (powf(storage.getDX(), 2));

	// Method for getting time == 1;
	//simply duplicating the initial conditions
	
	finalVector[1] = variables.getT0(numberOfX);
	/**/
	/*
	//Computing the first timestep with The laasonen method
	SetVariables varT0(variables, storage.getDT() * 2);
	DataStorage storT0(storage.getDX(), storage.getDT());
	laasonenImplicit(storT0, varT0, Thomas);
	finalVector[1] = storT0.getData()[1];
	/**/
	
	for (unsigned int time = 2; time < numberOfT; time++) {
		previousTime = finalVector[time - 1];
		vector<double> construct(numberOfX);
		for (unsigned int space = 0; space < numberOfX; space++) {
			//set boundary solution
			if (space == 0) {
				construct[space] = variables.getX0();
			}else if (space == numberOfX - 1) {
				construct[space] = variables.getXMax();
			}else{
				/*cout << space << "\n";
				cout << "size " << previousTime.size() << "\n";
				cout << "size " << finalVector.size() << "\n";*/
				
				//T(i-1,n)
				left = previousTime[space - 1];
				//T(i+1,n)
				right = previousTime[space + 1];
				//T(i,n-1)
				down = finalVector[time - 2][space];
				//compute T(i,n+1)
				construct[space] = (down + r * (right - down + left))/(1 + r);
			}
		}
		finalVector[time] = construct;

	}
	storage.setData(finalVector);
}

void richardson(DataStorage &storage, const SetVariables variables){
	/*
	Solve using the Richardson scheme
	/**/
	double right, left, central, down;
	unsigned int numberOfT = (variables.getT() / (double)storage.getDT())+1;
	unsigned int numberOfX = ((double)variables.getX() / (double)storage.getDX())+1;
	vector<vector<double>> finalVector(numberOfT);
	vector<double> previousTime(numberOfX);
	finalVector[0] = variables.getT0(numberOfX);
	double r = (variables.getD() * storage.getDT()) / (pow(storage.getDX(), 2));

	// Method of getting time == 1

	//simply duplicating the initial conditions
	finalVector[1] = variables.getT0(numberOfX);
	/**/
	
	//Computing the first timestep with The Crank-Niholson method
	SetVariables varT0(variables, storage.getDT()*2);
	DataStorage storT0(storage.getDX(), storage.getDT());
	crankNicholson(storT0, varT0, Thomas);
	finalVector[1] = storT0.getData()[1];
	/**/
	

	for (unsigned int time = 2; time < numberOfT; time++) {
		previousTime = finalVector[time - 1];
		vector<double> construct(numberOfX);
		for (unsigned int space = 0; space < numberOfX; space++) {
			//Set boundary conditions
			if (space == 0) {
				construct[space] = variables.getX0();
			}
			else if (space == numberOfX - 1) {
				construct[space] = variables.getXMax();
			}
			else {
				//T(i-1,n)
				left = previousTime[space - 1];
				//T(i+1,n)
				right = previousTime[space + 1];
				//T(i,n)
				central = previousTime[space];
				//T(i,n-1)
				down = finalVector[time - 2][space];
				//compute T(i,n+1)
				construct[space] = (down + 2 * r * (right - 2*central + left));
			}
		}
		finalVector[time] = construct;
	}
	storage.setData(finalVector);
}

void laasonenImplicit(DataStorage& storage, const SetVariables variables, vector<double>(*f)(vector<vector<double>>, vector<double>)) {
	/*
	Solve using the laasonen scheme
	this is an implicit scheme, at each timestep, all the points in space are computed by solving a linear system
	We formulate our problem as the matricial equation:
	A*T(n+1) = T(n)
	which we solve for T(n+1)
	/**/
	unsigned int numberOfT = (variables.getT() / (double)storage.getDT()) + 1;
	unsigned int numberOfX = ((double)variables.getX() / (double)storage.getDX()) + 1;
	vector<vector<double>> finalVector(numberOfT);
	vector<double> previousTime(numberOfX);
	finalVector[0] = variables.getT0(numberOfX);
	//compute the coeficients of the linear system
	double r = 1 * (variables.getD() * storage.getDT()) / powf(storage.getDX(), 2.0);
	double a = 2 * r + 1;
	double b = -r;

	for (int time = 1; time < numberOfT; time++) {

		// Construct matrix A
		vector<vector<double>> A(numberOfX);
		for (int timeS = 1; timeS < numberOfX - 1; timeS++) {
			vector<double> construct(numberOfX);
			for (int space = 0; space < numberOfX; space++) {
				//define tri-diagonal matrix
				if (space == timeS) {
					construct[space] = a;
				}
				else if (space == (timeS + 1) || space == (timeS - 1)) {
					construct[space] = b;
				}
				else {
					construct[space] = 0;
				}
			}
			A[timeS] = construct;
		}
		//extend A to include the boundary conditions into the resolution of the system
		vector<double> temp1(numberOfX), temp2(numberOfX);
		for (int i = 0; i < numberOfX; i++) {
			if (i == 0) {
				temp1[i] = 1;
				temp2[i] = 0;
			}
			else if (i == numberOfX - 1) {
				temp1[i] = 0;
				temp2[i] = 1;
			}
			else {
				temp1[i] = 0;
				temp2[i] = 0;
			}
		}
		A[0] = temp1;
		A[A.size() - 1] = temp2;

		for (int time = 1; time < numberOfT; time++) {

			// Construct vector b
			for (int i = 0; i < numberOfX; i++) {
				previousTime[i] = finalVector[time - 1][i];
			}

			// Solver
			vector<double> constr = (f)(A, previousTime);
			//compute all the unknown of T(i,n+1)
			vector<double> construct(numberOfX);

			for (int i = 0; i < numberOfX; i++) {
				//set the boundary conditions
				if (i == 0) {
					construct[i] = variables.getX0();
				}
				else if (i == numberOfX - 1) {
					construct[i] = variables.getXMax();
				}
				//copy the results of the resolution
				else {
					construct[i] = constr[i];
				}
			}
			finalVector[time] = construct;
		}
		storage.setData(finalVector);
	}
}

void crankNicholson(DataStorage &storage, const SetVariables variables, vector<double>(*f)(vector<vector<double>>, vector<double>)){
	/*
	Solve using the Crank-Nicholson scheme
	this is an implicit scheme, at each timestep, all the points in space are computed by solving a linear system
	We formulate our problem as the matricial equation:
	A*T(n+1) = B*T(n)
	which we solve for T(n+1)
	/**/
	unsigned int numberOfT = variables.getT() / (double)storage.getDT() + 1;
	unsigned int numberOfX = (double)variables.getX() / (double)storage.getDX() + 1;
	vector<vector<double>> finalVector(numberOfT);
	vector<double> previousTime = variables.getT0(numberOfX);
	finalVector[0] = previousTime;
	//compute r
	double r = variables.getD() * storage.getDT() /(2 * pow(storage.getDX(), 2));
	//compute the coeficients of the linear system
	double a = 1 + 2 * r;
	double b = 1 - 2 * r;
	double c = - r;

	// Construct matrix A and B
	vector<vector<double>> A(numberOfX);
	vector<vector<double>> B(numberOfX);

	for (int timeStep = 0; timeStep < numberOfX; timeStep++) {
		vector<double> constructA(numberOfX);
		vector<double> constructB(numberOfX);
		if (timeStep == 0) {
			constructA[timeStep] = 1;
			constructB[timeStep] = 1;
			for (int space = 1; space < numberOfX; space++) {
				constructA[space] = 0;
				constructB[space] = 0;
			}
		}
		else if (timeStep == numberOfX - 1) {
			constructA[timeStep] = 1;
			constructB[timeStep] = 1;
			for (int space = 0; space < numberOfX - 1; space++) {
				constructA[space] = 0;
				constructB[space] = 0;
			}
		}
		else {
			for (int space = 0; space < numberOfX; space++) {
				if (space == timeStep) {
					constructA[space] = a;
					constructB[space] = b;
				}
				else if (space == (timeStep + 1) || space == (timeStep - 1)) {
					constructA[space] = c;
					constructB[space] = -c;
				}
				else {
					constructA[space] = 0;
					constructB[space] = 0;
				}
			}
		}
		A[timeStep] = constructA;
		B[timeStep] = constructB;
	}

	for (int time = 1; time < numberOfT; time++) {

		// Construct vector b
		previousTime = finalVector[time - 1];
		// b = B * T(n)
		vector <double> righthandSide = dotProduct(previousTime, B);

		// Solver
		vector<double> construct = (f)(A, righthandSide);

		finalVector[time] = construct;
	}
	storage.setData(finalVector);
}

void analytical(DataStorage &storage, const SetVariables variables) {
	/*
	Compute the analytical solution for the specified discretization step and problem parameters
	/**/
	unsigned int numberOfT = (variables.getT() / (double)storage.getDT())+1;
	unsigned int numberOfX = ((double)variables.getX() / (double)storage.getDX())+1;
	vector<vector<double>> finalVector(numberOfT);
	finalVector[0] = variables.getT0(numberOfX);
	vector<double> previousTime;

	for (unsigned int time = 1; time < numberOfT; time++) {
		previousTime = finalVector[time-1];
		vector<double> construct(numberOfX);
		for (unsigned int space = 0; space < numberOfX; space++) {
			double sum{ 0 };
			for (int i = 1; i < 200; i++) {
				sum += exp(-variables.getD() * powf((i * PI / variables.getX()), 2.0) * storage.getDT() * time) * (1 - (pow(-1, i))) * sin(i * PI * storage.getDX() * space / variables.getX()) / (i * PI);
			}
			construct[space] = variables.getX0() + 2 * (38.0 - variables.getX0())* sum;
		}

		finalVector[time] = construct;
	}
	storage.setData(finalVector);
}
