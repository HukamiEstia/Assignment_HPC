#include "Parameters.h"

Parameters::Parameters(double L, double D, double t_max, double T_sur, double T_in){
	Parameters::L = L;
	Parameters::D = D;
	Parameters::t_max = t_max;
	Parameters::T_sur = T_sur;
	Parameters::T_in = T_in;
}

Parameters::Parameters(Parameters other, double t) {
	Parameters::L = other.L;
	Parameters::D = other.D;
	Parameters::t_max = t_max;
	Parameters::T_sur = other.T_sur;
	Parameters::T_in = other.T_in;
}

double Parameters::getL() const {
	return L;
}

double Parameters::getD() const {
	return D;
}

double Parameters::getDuration() const {
	return t_max;
}

double Parameters::getT_sur() const {
	return T_sur;
}

double Parameters::getT_in() const {
	return T_in;
}

vector<double> Parameters::getTime_0(int size) const {
	vector<double> data;
	data.push_back(T_sur);
	for (int i = 1; i < size; i++) {
		data.push_back(T_in);
	}
	data.push_back(T_in);
	return data;
}
