#pragma once
#ifndef MFE_H
#define MFE_H
#include <vector>
#include <functional>
#include <fstream>
#include <algorithm>
#include <iomanip>

typedef std::vector <std::pair<double, double>> grid;
typedef std::vector<size_t> pvector;
typedef std::vector<double> mvector;
typedef std::vector<std::vector<double>> matrix;
typedef std::vector<std::vector<size_t>> finiteElem;
typedef std::function<double(double, double) > ScalFunc2D;

enum exceptions { BAD_READ, OUT_OF_AREA, BAD_WRITE };

class mfe
{
	struct material
	{
		double K;
		double Phi;
		double S2;
	};

	struct _bc2
	{
		size_t n_i;
		size_t loc_node1_i;
		size_t loc_node2_i;
		double Tetta_i;
	};

public:
	mfe();

	//========================================================================
	void buildPortraitOfMatrix();
	void buildLocalG(size_t ielem);
	void addElementToGlobal(size_t i, size_t j, double elem);
	void assemblyGlobalMatrix();

	void toDense(const std::string _dense);

	double psi(size_t elem, size_t func, double r, double z);
	double d_psi(size_t elem, size_t func, double r, double z, size_t var);

	void bc_1();
	void bc_2();

	//========================================================================
	void mult(mvector& x, mvector& y);
	void LOS( );
	void MSG( );
	void MSG_ch();
	void calc_cholesky();
	void writeToFile(mvector q);
	double EuclideanNorm(mvector& x);

	//========================================================================
	double integration(ScalFunc2D f, double rk, double rk1, double zk, double zk1);


public:
	matrix G;					// Матрица жесткости
	matrix alfa;
	mvector q;					// Вектор весов
	mvector F;					// Вектор правой части

	size_t Nuz;					// Размер сетки
	grid MeshRZ;				// Сетка

	size_t Nel;					// Количество КЭ
	finiteElem FE;				// Конечные элементы

	size_t Nph;					// Количество фаз
	mvector Mu;					// Вязкости 


	std::vector<material> mats; // Материалы
	double lambda;

	size_t Nbc1;
	std::vector <std::pair<size_t, double>> bc1;

	size_t Nbc2;
	std::vector<_bc2> bc2;

	// Глобальная матрица
	mvector di;
	mvector cdi;
	mvector gg;
	mvector cgg;
	pvector ig;
	pvector jg;

	matrix mat;

	bool isOrdered(const pvector& v);

	double calcSum(size_t ielem);

	mvector r;
	mvector z;
	mvector p;
	size_t maxIter;
	double eps;


	mvector points;
	mvector weights;
};

inline mvector operator+(const mvector& a, const mvector& b)
{
	mvector res = a;
	for (size_t i = 0; i < res.size(); i++)
		res[i] += b[i];
	return res;
}

inline mvector operator-(const mvector& a, const mvector& b)
{
	mvector res = a;
	for (size_t i = 0; i < res.size(); i++)
		res[i] -= b[i];
	return res;
}

inline double operator*(const mvector& a, const mvector& b)
{
	double scalar = 0.0;
	for (size_t i = 0; i < a.size(); i++)
		scalar += a[i] * b[i];
	return scalar;
}

inline mvector operator*(double c, const mvector& a)
{
	std::vector<double> res = a;
	for (size_t i = 0; i < res.size(); i++)
		res[i] *= c;
	return res;
}

#endif