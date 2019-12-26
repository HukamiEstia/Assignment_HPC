#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <math.h> 
using namespace std;

class DataStorage
{
	/*
	Container to hold the computed results and write them to a file
	/**/
private:
	vector<vector<double>> data;
	double dX, dT;
public:
	DataStorage();
	DataStorage(double dX, double dT);
	void setDX(double dX);
	void setDT(double dT);
	void setData(vector<vector<double>> data);

	vector<vector<double>> getData() const;
	double getDX() const;
	double getDT() const;

	void printData();
	void saveDataToFile(string fileName);
	void saveSpecifiedTimeDataToFile(string fileName);
	DataStorage operator-(DataStorage const& obj);
};

