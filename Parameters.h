#pragma once
#include <vector>
using namespace std;
class Parameters
{
/*
Container to hold the parameters of the problem
- initial conditions
- boundary conditions
- diffusivity
/**/
private:
	double L, D, t_max, T_sur, T_in;
public:
	Parameters(double L, double D, double t_max, double T_sur, double T_in);
	Parameters(Parameters other, double t);
	double getL() const;
	double getDuration() const;
	double getD() const;
	double getT_sur() const;
	double getT_in() const;
	vector<double> getTime_0(int size) const;
};

