#include "compute.h"

using namespace std;

int main()
{
	compute calc;

	std::pair<double, double> point = calc.read_point();

	double pressure = calc.pressureInPoint(point.first, point.second);

	calc.write(pressure);

	return 0;
}