#include "compute.h"

compute::compute()
{
	std::ifstream inGrid("node.txt");
	if (inGrid.is_open())
	{
		inGrid >> Nuz;
		MeshRZ.resize(Nuz);

		double x, y;
		for (size_t i = 0; i < Nuz; i++)
		{
			inGrid >> x >> y;
			MeshRZ[i] = std::make_pair(x, y);
		}
		inGrid.close();
	}
	else throw 1;

	std::ifstream inElem("elem.txt");
	if (inElem.is_open())
	{
		inElem >> Nel;
		FE.resize(Nel);
		for (size_t i = 0; i < Nel; i++)
			FE[i].resize(4);

		size_t n1, n2, n3, n4;
		for (size_t i = 0; i < Nel; i++)
		{
			inElem >> n1 >> n2 >> n3 >> n4;
			FE[i] = { n1, n2, n3, n4 };
		}

		inElem.close();
	}
	else throw 1;

	std::ifstream qw("q.txt");
	if (qw.is_open())
	{
		qw >> Nuz;
		q.resize(Nuz);

		for (size_t i = 0; i < Nuz; i++)
			qw >> q[i];

		qw.close();
	}
	else throw 1;
}

std::pair<double, double> compute::read_point()
{
	std::ifstream p("point.txt");
	if (p.is_open())
	{
		double r, z;
		p >> r >> z;
		p.close();

		return std::make_pair(r, z);
	}
	else throw 2;
}

void compute::write(double val)
{
	std::ofstream out("point.txt");
	if (out.is_open())
	{
		out << val;
		out.close();
	}
	else throw 3;
}

// Найти элемент, в который попадает заданная точка
size_t compute::findElem(double r, double z)
{
	double rk, zk, rk1, zk1;
	for (size_t i = 0; i < Nel; i++)
	{
		rk = MeshRZ[FE[i][0]].first;
		zk = MeshRZ[FE[i][0]].second;

		rk1 = MeshRZ[FE[i][3]].first;
		zk1 = MeshRZ[FE[i][3]].second;

		if (r >= rk && r <= rk1 && z >= zk && z <= zk1)
			return i;
	}
	throw 0;
}

// Выдать давление в любой заданной точке
double compute::pressureInPoint(double r, double z)
{
	size_t ielem = findElem(r, z);

	double rk = MeshRZ[FE[ielem][0]].first;
	double zk = MeshRZ[FE[ielem][0]].second;

	double rk1 = MeshRZ[FE[ielem][3]].first;
	double zk1 = MeshRZ[FE[ielem][3]].second;

	std::vector<ScalFunc2D> f(4);
	f[0] = [rk, zk, rk1, zk1](double r, double z)
	{
		return (rk1 - r) / abs(rk1 - rk) * (zk1 - z) / abs(zk1 - zk);
	};
	f[1] = [rk, zk, rk1, zk1](double r, double z)
	{
		return (r - rk) / abs(rk1 - rk) * (zk1 - z) / abs(zk1 - zk);
	};
	f[2] = [rk, zk, rk1, zk1](double r, double z)
	{
		return (rk1 - r) / abs(rk1 - rk) * (z - zk) / abs(zk1 - zk);
	};
	f[3] = [rk, zk, rk1, zk1](double r, double z)
	{
		return (r - rk) / abs(rk1 - rk) * (z - zk) / abs(zk1 - zk);
	};

	double pressure = 0;
	for (uint8_t i = 0; i < 4; i++)
		pressure += q[FE[ielem][i]] * f[i](r, z);

	return pressure;
}