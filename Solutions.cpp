#include "Solutions.h"


void analytical(DataStorage &storage, const Parameters params) {
	/*
	Compute the analytical solution for the specified discretization step and problem parameters
	/**/
	unsigned int numberOfT = (params.getDuration() / (double)storage.getDT()) + 1;
	unsigned int numberOfX = ((double)params.getL() / (double)storage.getDX()) + 1;
	vector<vector<double>> solutionMat(numberOfT);
	solutionMat.push_back(params.getTime_0(numberOfX));

	for (unsigned int n = 1; n < numberOfT; n++) {
		vector<double> currentTime(numberOfX);
		for (unsigned int space = 0; space < numberOfX; space++) {
			double sum{ 0 };
			for (int i = 1; i < 200; i++) {
				sum += exp(-params.getD()
						* powf((i * PI / params.getL()), 2.0)
						* storage.getDT() * n) 
					* (1 - (pow(-1, i))) 
					* sin(i * PI * storage.getDX() * space / params.getL())
					/ (i * PI);
			}
			currentTime.push_back(params.getT_in() + 2 * (38.0 - params.getT_in())* sum);
		}

		solutionMat.push_back(currentTime);
	}
	storage.setData(solutionMat);
}

void FTCS(DataStorage &storage, const Parameters params) {
	/*
	This method is forward time central space
	/**/
	double left, right, central;
	
	//set dT so that the scheme is stable
	double dT = sqrt(pow(storage.getDX(), 2) / (params.getD()*2));

	//compute the number of time and space discretisation steps
	unsigned int numberOfT = params.getDuration() / storage.getDT() + 1;
	unsigned int numberOfX = params.getL() / storage.getDX() + 1;

	//Initialiazing the matrix holding the solution
	vector<vector<double>> solutionMat;
	solutionMat.push_back(params.getTime_0(numberOfX));

	//define r = D.(dt/dx^2)
	double r = (params.getD() * dT) / (pow(storage.getDX(), 2));

	for (unsigned int n = 1; n < numberOfT; n++) {
		
		//Initializing vector holding the solution at t=n+1
		vector<double> currentTime(numberOfX);

		//set boundary conditions
		currentTime[0] = params.getT_sur();
		currentTime[numberOfX - 1] = params.getT_in();
		for (unsigned int i = 1; i < numberOfX - 1; i++) {
			////T(n, i-1)
			left = solutionMat[n - 1][i - 1];
			//T(n, i)
			central = solutionMat[n - 1][i];
			//T(n, i+1)
			right = solutionMat[n - 1][i+1];
			//compute T(n+1, i)
			currentTime[i] = central+r*(left-2*central+right);
			}
		solutionMat.push_back(currentTime);
	}
	storage.setData(solutionMat);
}

void laasonenImplicit(DataStorage& storage, const Parameters params, vector<double>(*f)(vector<vector<double>>, vector<double>)) {
	/*
	Solve using the laasonen scheme
	this is an implicit scheme, at each timestep, all the points in space are computed by solving a linear system
	We formulate our problem as the matricial equation:
	A*T(n+1) = T(n)
	which we solve for T(n+1)
	/**/
	unsigned int numberOfT = (params.getDuration() / (double)storage.getDT()) + 1;
	unsigned int numberOfX = ((double)params.getL() / (double)storage.getDX()) + 1;
	vector<vector<double>> finalVector(numberOfT);
	vector<double> previousTime(numberOfX);
	finalVector[0] = params.getTime_0(numberOfX);
	//compute the coeficients of the linear system
	double r = 1 * (params.getD() * storage.getDT()) / powf(storage.getDX(), 2.0);
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
					construct[i] = params.getT_sur();
				}
				else if (i == numberOfX - 1) {
					construct[i] = params.getT_sur();
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

void crankNicholson(DataStorage &storage, const Parameters params, vector<double>(*f)(vector<vector<double>>, vector<double>)){
	/*
	Solve using the Crank-Nicholson scheme
	this is an implicit scheme, at each timestep, all the points in space are computed by solving a linear system
	We formulate our problem as the matricial equation:
	A*T(n+1) = B*T(n)
	which we solve for T(n+1)
	/**/
	unsigned int numberOfT = params.getDuration() / (double)storage.getDT() + 1;
	unsigned int numberOfX = (double)params.getL() / (double)storage.getDX() + 1;
	vector<vector<double>> finalVector(numberOfT);
	vector<double> previousTime = params.getTime_0(numberOfX);
	finalVector[0] = previousTime;
	//compute r
	double r = params.getD() * storage.getDT() /(2 * pow(storage.getDX(), 2));
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