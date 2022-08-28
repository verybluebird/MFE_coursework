#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;

class Grid
{
public:
	void input(); // read data
	void nodes();
	void out_coords();
	void elems();
	void material();
	void gr_bc1();
	void gr_bc2();
	void nested_grid(vector<double>& coord); //doubles the number of nodes M times
	void add_if_not_exist_and_sort(double L); // add a new node to the grid if it doesn't exist
	void print_profile(); // print matrix profile
	void time();
private:
	/*calculation area parameters*/
	double Rw; //initial coordinate in r (left)  
	double Rb; //final coordinate in r (right) 
	int NL; //number of layers in z  
	vector<double> H; // thickness of layers in z
	vector<double> K; //structural permeability of layers
	vector<double> Phi; // porosity of layers
	vector<double> S2; //oil saturation of layers
	
	int nr; //number of partitions over r
	vector<int> nz; //number of partitions over z in each layer
	double kr; //discharging coefficient for r 
	int M; // grid nesting in space 

	int Nzp; //number of perforation zones  
	vector<double> Pu; //upper coordinate of the perforation zone
	vector<double> Pd; //lower coordinate of the perforation zone
	vector<double> Tetta; // power of the perforation zone

	int Nph; // number of phases  
	vector<double> Mu; //viscosity of phases

	double Plast; //reservoir pressure

	vector<vector<int>> global_numbers; //global node numbers
	double z0 = 0; //where the first layer from above start in z
	vector<double> r_coord; // r-coordinates 
	vector<double> z_coord; // z-coordinates 
	
	vector<double> time_coord; // t-coordinates 
	int M_t; // grid nesting in time
	int n_t; //number of time layers
	int t1; //start time 
	int t2; //end time 
	double k_t; //coefficient of discharging over time 

	int Nuz; //number of nodes
	vector<vector<double>> coords;//{r_i z_i}
	int Nel; //
	vector<vector<int>> numbers; // global numbers
	vector<vector<double>> materials;
	int Nbc1; //number of first boundary conditions
	vector<double> bc1;
	int Nbc2;//number of second boundary conditions
	vector<vector<double>> bc2;

};