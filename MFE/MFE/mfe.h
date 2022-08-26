#pragma once
#ifndef MFE_H
#define MFE_H
#include <vector>
#include <functional>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <iostream>

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
	void buildLocalG(int ielem);
	void buildLocalG(int ielem, double r1, double r2, double r3, double z1, double z2, double z3, double detD);
	void buildLocalG(int ielem, double r[3], double z[3], double detD);
	void buildLocalF(int ielem, double r[3], double z[3], double detD,
		double dt, double t, mvector& q0, mvector& q1, mvector& q2);
	void buildLocalF(int ielem, double r[3], double z[3], double detD, double t, mvector& q0, mvector& q1, mvector& q2, double diffn1, double diffn2, double diffn3);
	void buildLocalF(int ielem, double dt, double t, double q0[3]);
	void buildLocalF(int ielem);
	void addElementToGlobal(int i, int j, double elem);

	double calcSum(int ielem);
	void assemblyGlobalMatrix();
	void assemblyGlobalMatrix(double t, double dt, std::vector<double>& q_);
	void assemblyGlobalMatrix(double t, double t1, double t2, double t3,
		std::vector<double>& q_0, std::vector<double>& q_1, std::vector<double>& q_2);
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
	void LOS();
	void MSG( );
	void writeToFile(mvector& q);
	void writeToFile(mvector& q, double t);
	double EuclideanNorm(mvector& x);

private:
	matrix G;	// stiffness matrix
	matrix alfa;
	matrix c;
	matrix M; //mass matrix

	mvector q;	// current solution
	mvector q0; //previous solution at t-1
	mvector q1; //solution at t-2
	mvector q2; //solution at t-3
	//mvector p; // next iteration vector
	//mvector p0; //current iteration vector
	//mvector p_1;
	//mvector Au;

	mvector b_loc; // local right side vector
	mvector p_loc;
	mvector F;	//Right side vector
	
	std::vector<int> bc1nodes;
	std::vector<double> b_2;
	std::vector<double> b_3;
	std::vector<double> ub;
	matrix A3;

	int Nuz;					// Grid Size
	grid MeshRZ;				// grid
	std::vector<double> r_coord; 
	std::vector<double> z_coord;
	int Rsize;
	int Zsize;

	int Nel;					// Number of FE
	finiteElem FE;				// Finite Elements

	int Nph;					// Number of phases
	mvector Mu;					// Viscosities


	std::vector<material> mats; // Materials
	std::vector<std::vector<int>> list; 

	int Nbc1;
	std::vector <std::pair<int, double>> bc1;

	int Nbc2;
	int Nbc3;
	std::vector<_bc> bc2;
	std::vector<_bc> bc3;

	// Global matrix
	mvector di;
	mvector gg;
	mvector gu;
	pvector ig;
	pvector jg;

	matrix mat;

	bool isOrdered(const pvector& v);

	int maxIter;
	double eps;

	mvector time;
	int n;

	//for CGM
	mvector um;
	mvector r;
	mvector z;

};

inline mvector operator+( mvector& a, const mvector& b)
{
	mvector res = a;
	for (int i = 0; i < a.size(); i++)
		res[i] += b[i];
	return res;
}

inline mvector operator-( mvector& a, const mvector& b)
{
	mvector res = a;
	for (int i = 0; i < a.size(); i++)
		res[i] -= b[i];
	return res;
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