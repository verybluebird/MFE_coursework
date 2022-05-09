#pragma once
#ifndef MFE_H
#define MFE_H
#include <vector>
#include <functional>
#include <fstream>
#include <algorithm>
#include <iomanip>

typedef std::vector <std::pair<double, double>> grid;
typedef std::vector<int> pvector;
typedef std::vector<double> mvector;
typedef std::vector<std::vector<double>> matrix;
typedef std::vector<std::vector<int>> finiteElem;
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

	struct _bc
	{
		int n_i;
		int loc_node1_i;
		int loc_node2_i;
		double Tetta_i;
	};

public:
	mfe();
	void iterationProcess();
	//========================================================================
	void buildPortraitOfMatrix();
	void buildLocalG(int ielem, double r1, double r2, double r3, double z1, double z2, double z3, double detD);
	
	void buildLocalF(int ielem, double r1, double r2, double r3, double z1, double z2, double z3, double detD, double dt, double t, double q0[3]);
	void buildLocalF(int ielem, double r1, double r2, double r3, double z1, double z2, double z3, double detD);
	void buildLocalF(int ielem, double dt, double t, double q0[3]);
	void addElementToGlobal(int i, int j, double elem);
	void assemblyGlobalMatrix();
	void assemblyGlobalMatrix(double t, double dt, std::vector<double>& q_);
	void toDense(const std::string _dense);

	double rightPart(int field, double r, double z);
	double rightPart(int field, double r, double z, double t);
	double Lambda(int field);
	double u_beta(double r, double z);
	double u_t(double r, double z, double t);
	double sigma(int field);

	void make_bc_1(double t);
	void bc_1();
	void bc_2();
	void bc_3();

	//========================================================================
	void mult(mvector& x, mvector& y);
	void MSG( );
	void writeToFile(mvector& q);
	void writeToFile(mvector& q, double t);
	double EuclideanNorm(mvector& x);

public:
	matrix G;	// Матрица жесткости
	matrix alfa;
	matrix c;
	matrix M;

	mvector q;	// текущее решение
	mvector q0; //предыдущее решение
	mvector p; // следующий вектор итерации
	mvector p0;//текущий вектор итерации
	mvector p_1;
	mvector Au;

	mvector b_loc;
	mvector p_loc;
	mvector F;	// Вектор правой части
	
	std::vector<int> bc1nodes;
	std::vector<double> b_2;
	std::vector<double> b_3;
	std::vector<double> ub;
	matrix A3;

	int Nuz;					// Размер сетки
	grid MeshRZ;				// Сетка
	std::vector<double> r_coord; 
	std::vector<double> z_coord;
	int Rsize;
	int Zsize;

	int Nel;					// Количество КЭ
	finiteElem FE;				// Конечные элементы

	int Nph;					// Количество фаз
	mvector Mu;					// Вязкости 


	std::vector<material> mats; // Материалы
	std::vector<std::vector<int>> list; 

	int Nbc1;
	std::vector <std::pair<int, double>> bc1;

	int Nbc2;
	int Nbc3;
	std::vector<_bc> bc2;
	std::vector<_bc> bc3;

	// Глобальная матрица
	mvector di;
	mvector gg;
	mvector gu;
	pvector ig;
	pvector jg;

	matrix mat;

	bool isOrdered(const pvector& v);

	double calcSum(int ielem);


	int maxIter;
	double eps;

	mvector time;
	int n;

	//для мсг
	mvector um;
	mvector r;
	mvector z;

};

inline mvector operator+( mvector& a, const mvector& b)
{
	for (int i = 0; i < a.size(); i++)
		a[i] += b[i];
	return a;
}

inline mvector operator-( mvector& a, const mvector& b)
{
	for (int i = 0; i < a.size(); i++)
		a[i] -= b[i];
	return a;
}

inline double operator*(const mvector& a, const mvector& b)
{
	double scalar = 0.0;
	for (int i = 0; i < a.size(); i++)
		scalar += a[i] * b[i];
	return scalar;
}

inline mvector operator*(double c,  mvector& a)
{
	for (int i = 0; i < a.size(); i++)
		a[i] *= c;
	return a;
}

#endif