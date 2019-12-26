#include "SetVariables.h"

SetVariables::SetVariables(double x, double D, double t, double x0, double xmax, double t0){
	SetVariables::x = x;
	SetVariables::D = D;
	SetVariables::t = t;
	SetVariables::x0 = x0;
	SetVariables::xmax = xmax;
	SetVariables::t0 = t0;
}

SetVariables::SetVariables(SetVariables other, double t) {
	SetVariables::x = other.x;
	SetVariables::D = other.D;
	SetVariables::t = t;
	SetVariables::x0 = other.x0;
	SetVariables::xmax = other.xmax;
	SetVariables::t0 = other.t0;
}

double SetVariables::getX() const {
	return x;
}

double SetVariables::getD() const {
	return D;
}

double SetVariables::getT() const {
	return t;
}

double SetVariables::getX0() const {
	return x0;
}

double SetVariables::getXMax() const {
	return xmax;
}

vector<double> SetVariables::getT0(int size) const {
	vector<double> data;
	for (int i = 0; i < size; i++) {
		data.push_back(t0);
	}
	data[0] = x0;
	data[data.size() - 1] = xmax;
	return data;
}
