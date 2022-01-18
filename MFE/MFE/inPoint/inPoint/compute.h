#pragma once
#include <vector>
#include <functional>
#include <string>
#include <fstream>

typedef std::function<double(double, double) > ScalFunc2D;

class compute
{
public:
	compute();
	std::pair<double, double> read_point();
	void write(double val);
	double pressureInPoint(double r, double z);

private:
	std::vector<double> q;
	std::vector <std::pair<double, double>> MeshRZ;
	std::vector<std::vector<size_t>> L;
	std::vector<std::vector<size_t>> FE;
	size_t Nel;
	size_t Nuz;

	size_t findElem(double r, double z);
};

