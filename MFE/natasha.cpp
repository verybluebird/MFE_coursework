#include <iostream>
#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <iomanip>
#include <fstream>
#include "Loss.h"
#include "losv2.h"

using namespace std;

typedef vector<vector<double>> Matrixx;


class Net
{
public:
	Net()
	{
	}
	Net(fstream& nodes, fstream& elements, fstream& fields, fstream& condi1, fstream& condi2, fstream& condi3)
	{
		double r, z;
		int a, b, c, field, type;
		while (nodes >> r >> z)
		{
			Node.push_back({ r, z });
		}
		while (elements >> a >> b >> c)
		{
			Elements.push_back({ a,b,c });
		}
		while (fields >> field)
		{
			this->fields.push_back(field);
		}
		while (condi1 >> a >> type)
		{
			this->FirstCond.push_back({ a,type });
		}
		while (condi2 >> a >> b >> type)
		{
			this->SecondCond.push_back({ a, b,type });
		}
		while (condi3 >> a >> b >> type >> field)
		{
			this->ThirdCond.push_back({ a,b,type, field });
		}
	}
	
	void SaveNet(fstream& nodes, fstream& elements, fstream& fields)
	{
		int length = Node.size();
		for (size_t i = 0; i < length; i++)
		{
			nodes << Node[i][0] << " " << Node[i][1] << "\n";
		}
		length = this->Elements.size();
		for (size_t i = 0; i < length; i++)
		{
			elements << this->Elements[i][0]+1 << " " << this->Elements[i][1]+1 << " " << this->Elements[i][2]+1 << "\n";
			fields << this->fields[i] << "\n";
		}
	}
	vector<vector<double>> Node;
	vector<vector<int>> Elements;
	vector<int> fields;
	vector<vector<int>> FirstCond;
	vector<vector<int>> SecondCond;
	vector<vector<int>> ThirdCond;
	
	void DevideBy2Fields()//test
	{
		fields = vector<int>(Elements.size());
		int middle = fields.size() / 2;
		for (int i = middle; i < fields.size(); i++)
		{
			fields[i] = 1;
		}
	}
	
	void AddFirstCond(int nr, int nz)
	{
		for (int j = 0; j < nz; j++)
		{
			for (int i = 0; i < nr; i++)
			{
				if (/*j == 0 || j == nz - 1 || i == 0 || i == nr - 1*//*j == 0 && i==0 || j == nz - 1 && i==0 || i == nr-1 && j == 0 || i == nr - 1 && j == nz - 1*/i == 0 || i == nr-1)
				{
					int k = nr * j + i;
					FirstCond.push_back({ k,0 });
				}
			}
		}
	}
	void AddSecondCond(int nr, int nz)
	{
		for (int j = 0; j < nz; j++)
		{
			for (int i = 0; i < nr; i++)
			{
				if (j == nz - 1 && i != nr - 1)
				{
					int k1 = nr * j + i;
					int k2 = nr * j + i + 1;
					SecondCond.push_back({ k1,k2,0 });
				}
			}
		}
	}
	void AddThirdCond(int nr, int nz)
	{
		for (int j = 0; j < nz; j++)
		{
			for (int i = 0; i < nr; i++)
			{
				if (j == 0 && i != nr - 1)
				{
					int k1 = i;
					int k2 = i + 1;
					ThirdCond.push_back({ k1,k2,0 });
				}
			}
		}
	}
	void BuildNet(double xmin, double xmax, double ymin, double ymax, int nx, int ny)
	{
		double hx = (xmax - xmin) / nx;
		double hy = (ymax - ymin) / ny;
		Node = vector<vector<double>>((nx + 1) * (ny + 1));
		Node[0] = vector<double>{ xmin, ymin };
		for (int i = 0; i < ny; i++)
		{
			double y = ymin + i * hy;
			for (int j = 0; j < nx; j++)
			{
				double x = xmin + j * hx;
				Node[i * (nx + 1) + j + 1] = { x + hx, y };
				Node[(i + 1) * (nx + 1) + j] = { x,y + hy };
				Elements.push_back({ j + i * (nx + 1),j + 1 + i * (nx + 1), j + (nx + 1) * (i + 1) });
			}
		}
		Node[Node.size() - 1] = { xmax,ymax };
		for (int i = ny; i > 0; i--)
		{
			for (int j = nx; j > 0; j--)
			{
				Elements.push_back({ j + i * (nx + 1) - nx - 1,j - 1 + i * (nx + 1), j + i * (nx + 1) });
			}
		}
		int length = Elements.size();
		vector<vector<int>> Elementstmp(length);
		for (int j = 0, i = 0; i < length; j++, i += 2)
		{
			Elementstmp[i] = Elements[j];
		}
		for (int i = 1, j = length - 1; i < length; i += 2, j--)
		{
			Elementstmp[i] = Elements[j];
		}
		Elements = Elementstmp;
		fields.resize(Elements.size());
		
	}
};

class SLU
{
public:
	SLU()
	{
		TheNet.BuildNet(0, 2, 0, 2, 2, 2);
		b = vector<double>(TheNet.Node.size());
		q = vector<double>(TheNet.Node.size());
	}
	SLU(Net net)
	{
		TheNet = net;
		b = vector<double>(TheNet.Node.size());
		q = vector<double>(TheNet.Node.size());
	}
	~SLU() {}
	MatrixProf AProf;
	Net TheNet;
	vector<double> b;
	vector<double> q;
	void BuildProf()
	{

	}
	void BuildProfile()
	{
		vector<vector<int>> profile(TheNet.Node.size());

		for (int i = 0; i < TheNet.Elements.size(); i++)
		{
			for (int j = 1; j < 3; j++)
			{
				for (int k = 0; k < j; k++)
				{
					int current = TheNet.Elements[i][j];
					int node = TheNet.Elements[i][k];
					if (!count(profile[current].begin(), profile[current].end(), node))
					{
						if (profile[current].size() != 0 && pro-file[current][profile[current].size() - 1] > node)
						{
							for (int l = 0; l < profile[current].size(); l++)
							{
								if (node < profile[current][l])
								{
									pro-file[current].insert(profile[current].begin() + l, node);
									break;
								}
							}
						}
						else
						{
							profile[current].push_back(node);
						}
					}
				}
			}
		}

		AProf.IA.push_back(1);
		int count = 0;
		for (int i = 1; i < TheNet.Node.size(); i++)
		{
			AProf.IA.push_back(AProf.IA[i - 1] + count);
			count = 0;
			for (int j = 0; j < profile[i].size(); j++)
			{
				AProf.JA.push_back(profile[i][j]);
				count++;
			}
		}
		AProf.IA.push_back(AProf.IA[AProf.IA.size() - 1] + count);
		AProf.AL = vector<double>(AProf.IA[AProf.IA.size() - 1] - 1);
		AProf.AU = vector<double>(AProf.IA[AProf.IA.size() - 1] - 1);
		AProf.DI = vector<double>(TheNet.Node.size());
		AProf.size = AProf.DI.size();
	}

	vector<double> MVecMult(Matrixx& A, vector<double>& b)
	{
		vector<double> result(A.size());
		int length = A.size();
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < length; j++)
			{
				result[i] += A[i][j] * b[j];
			}
		}
		return result;
	}

	vector<double> VecSum(vector<double>& a, vector<double>& b)
	{
		int length = a.size();
		vector<double> result;
		for (size_t i = 0; i < length; i++)
		{
			result.push_back(a[i] + b[i]);
		}
		return result;
	}

	void LocalToGLobal(Matrixx& A, vector<double>& b, vector<int>& el)
	{
		int length = A.size();
		for (int i = 0; i < length; i++)
		{
			AProf.DI[el[i]] = AProf.DI[el[i]] + A[i][i];
		}

		for (int i = 0; i < length; i++)
		{
			int ibeg = AProf.IA[el[i]] - 1;
			for (int j = 0; j < i; j++)
			{
				int iend = AProf.IA[el[i] + 1] - 1;
				while (AProf.JA[ibeg] != el[j])
				{
					int ind = (ibeg + iend) / 2;
					if (AProf.JA[ind] <= el[j])
					{
						ibeg = ind;
					}
					else
					{
						iend = ind;
					}
				}
				AProf.AL[ibeg] += A[i][j];
				AProf.AU[ibeg] += A[j][i];
				ibeg++;
			}
		}

		for (int i = 0; i < length; i++)
		{
			this->b[el[i]] += b[i];
		}

	}

	vector<vector<double>> BuildG(vector<vector<double>>& D_1, double detD, vector<int>& el, int field)
	{
		vector<vector<double>> G(3);
		double r1 = TheNet.Node[el[0]][0];
		double r2 = TheNet.Node[el[1]][0];
		double r3 = TheNet.Node[el[2]][0];

		double multix = abs(detD) * (r1 + r2 + r3) / 6.;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{

				double L = Lambda(field);
				G[i].push_back(L * multix * (D_1[i][1] * D_1[j][1] + D_1[i][2] * D_1[j][2])); 
			}
		}
		return G;
	}

	double Factorial(int n)
	{
		return (n == 1 || n == 0) ? 1 : Factorial(n - 1) * n;
	}

	double GetIntegral(vector<double> & v, double detD)
	{
		double result = Factorial(v[0]) * Factorial(v[1]) * Factorial(v[2]) / (Factori-al(v[0] + v[1] + v[2] + 2)) * abs(detD);
		return result;
	}

	vector<vector<double>> BuildM(vector<vector<double>>& D_1, double detD, vector<int>& el, int field)
	{
		vector<double> r(6);
		vector<double> z(6);
		for (size_t i = 0; i < 3; i++)
		{
			r[i] = TheNet.Node[el[i]][0];
			z[i] = TheNet.Node[el[i]][1];
		}
		r[3] = (r[0] + r[1]) / 2;
		z[3] = (z[0] + z[1]) / 2;
		r[4] = (r[1] + r[2]) / 2;
		z[4] = (z[1] + z[2]) / 2;
		r[5] = (r[0] + r[2]) / 2;
		z[5] = (z[0] + z[2]) / 2;
		

		vector<vector<double>> M(3);
		for (size_t i = 0; i < 3; i++)
		{
			M[i] = vector<double>(3);
		}

		vector<double> v(3);
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				for (size_t l = 0; l < 6; l++)
				{
					for (size_t k = 0; k < 3; k++)
					{
						v.clear();
						v.resize(3);

						if (l<=2)
						{
							v[i]++;
							v[j]++;
							v[k]++;
							v[l]++;
							double tmp1 = -GetIntegral(v, detD);
							v[l]++;
							double tmp2 = 2 * GetIntegral(v, detD);
							M[i][j] += (tmp1 + tmp2)*Gamma(field,r[l],z[l])*r[k];
						}
						else 
						{
							v[i]++;
							v[j]++;
							v[k]++;
							if (l==3)
							{
								v[0]++;
								v[1]++;
							}
							if (l==4)
							{
								v[1]++;
								v[2]++;
							}
							if (l==5)
							{
								v[0]++;
								v[2]++;
							}
							double tmp1 = 4 * GetIntegral(v, detD);
							M[i][j] += tmp1 * Gamma(field, r[l], z[l]) * r[k];
						}
					}
				}
			}
		}
		return M;
	}

	Matrixx BuildC(double DetD, vector<int>& el)
	{
		double r1 = TheNet.Node[el[0]][0];
		double r2 = TheNet.Node[el[1]][0];
		double r3 = TheNet.Node[el[2]][0];

		Matrixx C = Matrixx{
			{ 6 * r1 + 2 * r2 + 2 * r3, 2 * r1 + 2 * r2 + r3, 2 * r1 + r2 + 2 * r3 },
			{ 2 * r1 + 2 * r2 + r3, 2 * r1 + 6 * r2 + 2 * r3, r1 + 2 * r2 + 2 * r3 },
			{ 2 * r1 + r2 + 2 * r3, r1 + 2 * r2 + 2 * r3,2 * r1 + 2 * r2 + 6 * r3}
		};
		double mult = abs(DetD) / 120;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				C[i][j] *= mult;
			}
		}
		return C;
	}

	Matrixx BuildLocal(vector<int>& el, int field)
	{
		double r1 = TheNet.Node[el[0]][0];
		double r2 = TheNet.Node[el[1]][0];
		double r3 = TheNet.Node[el[2]][0];
		double z1 = TheNet.Node[el[0]][1];
		double z2 = TheNet.Node[el[1]][1];
		double z3 = TheNet.Node[el[2]][1];
		Matrixx D 
		{
			{1,1,1},
			{r1,r2,r3},
			{z1,z2,z3}
		};
		Matrixx D_1
		{
			{r2* z3 - r3 * z2, z2 - z3, r3 - r2},
			{r3* z1 - r1 * z3, z3 - z1, r1 - r3},
			{r1* z2 - r2 * z1, z1 - z2, r2 - r1}
		};
		double detD = (r2 - r1) * (z3 - z1) - (r3 - r1) * (z2 - z1);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				D_1[i][j] /= detD;
			}
		}
		Matrixx G = BuildG(D_1, detD, el, field);
		Matrixx M = BuildC(detD, el);/*BuildM(D_1, detD, el, field);*/
		Matrixx C = BuildC(detD, el);
		
		vector<double> f = { F(r1,z1,field),F(r2,z2,field),F(r3,z3,field) };
		vector<double> b = MVecMult(C, f);
		int length = b.size();
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < length; j++)
			{
				M[i][j] = M[i][j]; //* Gamma(field,0,0);
				G[i][j] = G[i][j] + M[i][j];
			}
		}
		LocalToGLobal(G, b, el);
		return G;
	}

	void AddFirst()
	{
		double max = 0;
		int length = AProf.AL.size();
		for (int i = 0; i < length; i++)
		{
			if (max < abs(AProf.AL[i]))
			{
				max = abs(AProf.AL[i]);
			}
		}
		max *= 1e+30;
		length = TheNet.FirstCond.size();
		for (int i = 0; i < length; i++)
		{
			int n = TheNet.FirstCond[i][0];
			AProf.DI[n] = max;
			double r = TheNet.Node[n][0];
			double z = TheNet.Node[n][1];
			b[n] = max * Ug(TheNet.Node[n], TheNet.FirstCond[i][1]);
			q[n] = Ug(TheNet.Node[n], TheNet.FirstCond[i][1]);
		}
	}
	void AddThird()
	{
		int length = TheNet.ThirdCond.size();
		for (int i = 0; i < length; i++)
		{
			vector<int> Edge = TheNet.ThirdCond[i];
			double r1 = TheNet.Node[Edge[0]][0];
			double z1 = TheNet.Node[Edge[0]][1];
			double r2 = TheNet.Node[Edge[1]][0];
			double z2 = TheNet.Node[Edge[1]][1];
			double hm = sqrt((r2 - r1) * (r2 - r1) + (z2 - z1) * (z2 - z1));
			Matrixx MS3 = GetMS(r1, r2);
			double B = Betta((int)Edge[2]);
			double mult = B * hm / 24.;
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					MS3[i][j] = MS3[i][j] * mult;
				}
			}
			double UB1 = UB(TheNet.Node[Edge[0]], Edge[2]);
			double UB2 = UB(TheNet.Node[Edge[1]], Edge[2]);
			vector<double> Ub = { UB1,UB2 };
			vector<double> b = MVecMult(MS3, Ub);

			LocalToGLobal(MS3, b, Edge);
		}
	}
	void AddSecond()
	{
		int length = TheNet.SecondCond.size();
		for (int i = 0; i < length; i++)
		{
			vector<int> Edge = TheNet.SecondCond[i];
			double r1 = TheNet.Node[Edge[0]][0];
			double z1 = TheNet.Node[Edge[0]][1];
			double r2 = TheNet.Node[Edge[1]][0];
			double z2 = TheNet.Node[Edge[1]][1];
			double hm = sqrt((r2 - r1) * (r2 - r1) + (z2 - z1) * (z2 - z1));
			double mult = hm / 24;
			Matrixx M2 = GetMS(r1, r2);
			vector<double> Tet;
			Tet.push_back(Tetta(TheNet.Node[Edge[0]], Edge[2]));
			Tet.push_back(Tetta(TheNet.Node[Edge[1]], Edge[2]));
			vector<double> b = MVecMult(M2, Tet);
			this->b[Edge[0]] += b[0] * mult;
			this->b[Edge[1]] += b[1] * mult;
		}
	}
	Matrixx GetMS(double r1, double r2)
	{
		Matrixx M2;
		M2.push_back({ 6 * r1 + 2 * r2,2 * (r1 + r2) });
		M2.push_back({ 2 * (r1 + r2),2 * r1 + 6 * r2 });
		return M2;
	}
	void BuildGlobal()
	{
		for (size_t i = 0; i < TheNet.Elements.size(); i++)
		{
			BuildLocal(TheNet.Elements[i], TheNet.fields[i]);
		}
		
		//Добавить 3-е, вторые и первые
	}
	
	/*void Calculate(vector<double>& sol)
	{
		MSG solver;
		sol = solver.Solve(AProf, b);
	}*/

	

	void CalculateLossv2(vector<double>& sol)
	{
		Matrix A;
		A.AL = AProf.AL.data();
		A.AU = AProf.AU.data();
		A.DI = AProf.DI.data();
		A.IA = AProf.IA.data();
		A.JA = AProf.JA.data();
		A.N = AProf.size;

		Matrix LU;
		for (size_t i = 0; i < AProf.IA.size(); i++)
		{
			A.IA[i] -= 1;
		}
		AuxVectors aux;
		aux.Ax = new double[A.N];
		aux.p = new double[A.N];
		aux.r = new double[A.N];
		aux.z = new double[A.N];
		aux.temp = new double[A.N];
		aux.LU = new double[A.N];

		double* x = new double[A.N];
		double res;
		LUFactorization(A, LU);
		LOS_LU(A, q.data(), b.data(),LU, aux, 50000, 1e-15, res);
	}

	
	void FindSolution()
	{
		BuildProfile();
		BuildGlobal();
		AddThird();
		AddSecond();
		AddFirst();
		CalculateLossv2(q);
	}
	double U(vector<double>& node)
	{
		double r = node[0];
		double z = node[1];
		return z;
	}

	double Ug(vector<double> &node, int k)
	{
		double r = node[0];
		double z = node[1];
		return z;
	}


	double UB(vector<double>& node, int k)
	{
		double r = node[0];
		double z = node[1];
		return 1-Lambda(k)/Betta(k);
	}
	double F(double r, double z, int field)
	{
		//return -Lambda(field) / r + Gamma(field, r, z) * r;// U = r
		//return Gamma(field, r, z);
		return z * Gamma(field,r,z);
		//return -4*Lambda(field) + Gamma(field,r,z)*r*r;
		//return z*z*Gamma(field,r,z)-2*Lambda(field);//U=zz
		//return z*z*z*Gamma(field,r,z)-6*Lambda(field)*z;//U=zzz
		// 
		//return Gamma(field, r, z) * z * z+Gamma(field, r, z) * r * r -  6*Lambda(field); //U=z^2
		//return Lambda(field)*sin(r+z)-Lambda(field)*cos(r+z)+Lambda(field)*r*sin(r+z);
		//return Lambda(field) * sin(r + z) - Lambda(field) * cos(r + z) + Lambda(field) * r * sin(r + z) + Gamma(field, r, z) * sin(r + z);

	}

	double Lambda(int field)
	{
		return 1;
	}
	//входные параметры скопировать
	double Gamma(int field, double r, double z)
	{
		return 1;
	}

	double Tetta(vector<double>& node, int field)
	{
		double r = node[0];
		double z = node[1];
		return Lambda(field);
	}

	double Betta(int field)
	{
		return 1;
	}

private:

};



int main()
{
	fstream nodes;
	fstream elements;
	fstream fields;
	fstream condi1;
	fstream condi2;
	fstream condi3;
	fstream result;
	nodes.open("nodes.txt");
	elements.open("elements.txt");
	fields.open("fields.txt");
	condi1.open("condi1.txt");
	condi2.open("condi2.txt");
	condi3.open("condi3.txt");
	result.open("result.txt");
	int nr = 2 , nz = 3;
	Net net = Net();
	int n;
	//Net net = Net(nodes, elements, fields, condi1, condi2, condi3);
	net.BuildNet(1, 2, 1, 2, nr, nz);
	net.SaveNet(nodes, elements, fields);
	//net.AddThirdCond(nr + 1, nz + 1);
	//net.AddSecondCond(nr + 1, nz + 1);
	net.AddFirstCond(nr + 1, nz + 1);
	
	SLU eq(net);
	n=eq.Factorial(5);

	eq.FindSolution();
	cout << scientific << setprecision(15);
	double sum_dif = 0;
	double norm = 0;
	for (size_t i = 0; i < net.Node.size(); i++)
	{
		double tmp = eq.U(net.Node[i]);
		sum_dif += pow((tmp - eq.q[i]), 2);
		norm += tmp * tmp;
	}
	sum_dif = sqrt(sum_dif/norm);
	cout << sum_dif << endl;
	for (int j = 0; j < eq.q.size(); j++)
		cout << eq.q[j] <<' '<< eq.U(net.Node[j])<< ' ' << abs((eq.U(net.Node[j])-eq.q[j]))/abs(eq.U(net.Node[j])) << endl;
	
	
}






