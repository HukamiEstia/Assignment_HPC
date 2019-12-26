#include "DataStorage.h"

DataStorage::DataStorage(){
	dX = 0;
	dT = 0;
}

DataStorage::DataStorage(double dX, double dT){
	DataStorage::dX = dX;
	DataStorage::dT = dT;
}

void DataStorage::setDX(double dX){
	DataStorage::dX = dX;
}

void DataStorage::setDT(double dT){
	DataStorage::dT = dT;
}

void DataStorage::setData(vector<vector<double>> data){
	DataStorage::data = data;
}

vector<vector<double>> DataStorage::getData() const {
	return data;
}

double DataStorage::getDX() const{
	return dX;
}

double DataStorage::getDT() const{
	return dT;
}

void DataStorage::printData(){
	for (int i = 0; i < data.size(); i++) {
		vector<double> innerData = data[i];
		for (int j = 0; j < innerData.size(); j++) {
			cout << setprecision(4) << innerData[j] << " ";
		}
		cout << endl;
	}
}

void DataStorage::saveDataToFile(string fileName){
	ofstream file(fileName + ".dat");
	if (file.is_open()) {
		for (int i = 0; i < data.size(); i++) {
			vector<double> innerData = data[i];
			for (int j = 0; j < innerData.size(); j++) {
				if (isnan(innerData[j])) {
					cout << "nan: " << i << "," << j << "\n";
				}
				file << innerData[j] << " ";
			}
			file << endl;
		}
		file.close();
	}else {
		cout << "Unable to open file!" << endl;
	}
}

void DataStorage::saveSpecifiedTimeDataToFile(string fileName) {
	ofstream file(fileName + ".dat");
	if (file.is_open()) {
		for (int i = 0; i <= 5; i++) {
			vector<double> innerData = data[i*(0.1/dT)];
			for (int j = 0; j < innerData.size(); j++) {
				file << innerData[j] << " ";
			}
			file << endl;
		}
		file.close();
	}else {
		cout << "Unable to open file!" << endl;
	}

}

DataStorage DataStorage::operator-(DataStorage const& obj) {
	DataStorage created(dX - obj.getDX(), dT - obj.getDT());

	if (obj.getData().size() != data.size() || data.empty()) {
		cout << "ERROR: Outer vectors not of same size!";
		return created;
	}

	vector<vector<double>> finalVector;
	for (int i = 0; i < data.size(); i++) {
		vector<double> construct;
		for (int j = 0; j < data[i].size(); j++) {
			construct.push_back(data[i][j] - obj.getData()[i][j]);
		}
		finalVector.push_back(construct);
	}
	created.setData(finalVector);
	return created;
}
