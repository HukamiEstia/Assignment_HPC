#pragma once
#include <vector>
using namespace std;
class SetVariables
{
/*
Container to hold the parameters of the problem
- initial conditions
- boundary conditions
- diffusivity
/**/
private:
	double x, D, t, x0, xmax, t0;
public:
	SetVariables(double x, double D, double t, double x0, double xmax, double t0);
	SetVariables(SetVariables other, double t);
	double getX() const;
	double getT() const;
	double getD() const;
	double getX0() const;
	double getXMax() const;
	vector<double> getT0(int size) const;
};

