#include "../include/mfe.h"


//===========================================================================================
// Constructor: read the mesh, FE with their node numbers, 1st and 2nd boundary conditions, materials
// and phase properties, parameters for SLAE solution.
mfe::mfe()
{
	{
		FILE* file;
		fopen_s(&file, "data/q.txt", "w");
		fclose(file);
	}
	//-----------------------------------------------------
	std::ifstream inCoord("data/coords.txt");
	if (inCoord.is_open())
	{
		inCoord >> Rsize;
		r_coord.resize(Rsize);
		for (int i = 0; i < Rsize; i++)
		{
			inCoord >> r_coord[i];
		}
		inCoord >> Zsize;
		z_coord.resize(Zsize);
		for (int i = 0; i < Zsize; i++)
		{
			inCoord >> z_coord[i];
		}
		inCoord.close();
	}
	else throw BAD_READ;

	//-----------------------------------------------------
	std::ifstream inGrid("data/node.txt");
	if (inGrid.is_open())
	{
		inGrid >> Nuz;
		MeshRZ.resize(Nuz);

		double x, y;
		for (int i = 0; i < Nuz; i++)
		{
			inGrid >> x >> y;
			MeshRZ[i] = std::make_pair(x, y);
		}
		inGrid.close();
	}
	else throw BAD_READ;

	//-----------------------------------------------------
	// FE is defined by the numbers of its nodes
	std::ifstream inElem("data/elem.txt");
	if (inElem.is_open())
	{
		inElem >> Nel;
		FE.resize(Nel);
		for (int i = 0; i < Nel; i++)
			FE[i].resize(3);

		int n1, n2, n3;
		for (int i = 0; i < Nel; i++)
		{
			inElem >> n1 >> n2 >> n3;
			FE[i] = { n1, n2, n3 };
		}

		inElem.close();
	}
	else throw BAD_READ;

	//-----------------------------------------------------
	std::ifstream inMat("data/mat.txt");
	if (inMat.is_open())
	{
		mats.resize(Nel);

		double K, Phi, S2;
		for (int i = 0; i < Nel; i++)
		{
			inMat >> K >> Phi >> S2;
			mats[i].K = K;
			mats[i].Phi = Phi;
			mats[i].S2 = S2;
		}
		inMat.close();
	}
	else throw BAD_READ;

	//-----------------------------------------------------
	std::ifstream inPhase("data/phaseprop.txt");
	if (inPhase.is_open())
	{
		inPhase >> Nph;
		Mu.resize(Nph);
		
		for (int i = 0; i < Nph; i++)
			inPhase >> Mu[i];

		inPhase.close();
	}
	else throw BAD_READ;


	//-----------------------------------------------------
	std::ifstream bc1File("data/bc1.txt");
	if (bc1File.is_open())
	{
		bc1File >> Nbc1;
		bc1.resize(Nbc1);
		for (int i = 0; i < Nbc1; i++)
		{
			int n_i;
			double value_bc1;

			bc1File >> n_i >> value_bc1;
			bc1[i] = std::make_pair(n_i, value_bc1);
		}
		bc1File.close();
	}
	else throw BAD_READ;

	//-----------------------------------------------------
	std::ifstream bc2File("data/bc2.txt");
	if (bc2File.is_open())
	{
		bc2File >> Nbc2;
		bc2.resize(Nbc2);

		int el_i, loc_node1_i, loc_node2_i;
		double Tetta_i;

		for (int i = 0; i < Nbc2; i++)
		{
			bc2File >> el_i >> loc_node1_i >> loc_node2_i >> Tetta_i;
			bc2[i] = { el_i, loc_node1_i, loc_node2_i, Tetta_i };
		}
		bc2File.close();
	}
	else throw BAD_READ;
	//-----------------------------------------------------
	std::ifstream bc3File("data/bc3.txt");
	if (bc3File.is_open())
	{
		bc3File >> Nbc3;
		bc3.resize(Nbc3);

		int el_i, loc_node1_i, loc_node2_i;
		double Betta_i;

		for (int i = 0; i < Nbc3; i++)
		{
			bc3File >> el_i >> loc_node1_i >> loc_node2_i >> Betta_i;
			bc3[i] = { el_i, loc_node1_i, loc_node2_i, Betta_i };
		}
		bc3File.close();
	}
	else throw BAD_READ;
	//-----------------------------------------------------
	maxIter = eps = 0;
	std::ifstream slau("data/kuslau.txt");
	if (slau.is_open())
	{
		slau >> maxIter >> eps;
		slau.close();
	}
	else throw BAD_READ;

	//-----------------------------------------------------
	std::ifstream inTime("data/time.txt");
	if (inTime.is_open())
	{
		inTime >> n;
		time.resize(n);
		for (int i = 0; i < n; i++)
		{
			inTime >> time[i];
		}
		inTime.close();
	}
	else throw BAD_READ;
	//-----------------------------------------------------
	b_loc.resize(3);
	p_loc.resize(3);
	
		
	F.resize(Nuz);

	c.resize(3);
	M.resize(3);
	alfa.resize(3);
	G.resize(3);
	for (int i = 0; i < 3; i++) {
		c[i].resize(3);
		M[i].resize(3);
		alfa[i].resize(3);
		G[i].resize(3);
	}

		
	list.resize(Nuz);
	ig.resize(Nuz + 1);
	di.resize(Nuz);
	mat.resize(Nuz);
	for (int i = 0; i < mat.size(); i++)
		mat[i].resize(Nuz, 0);
	bc1nodes.resize(Nuz, -1);

	b_2.resize(2);
	b_3.resize(2);
	ub.resize(2);

	A3.resize(2);
	A3[0].resize(2);
	A3[1].resize(2);

	um.resize(Nuz);
	z.resize(Nuz);
	r.resize(Nuz);
	q.resize(Nuz);
	q0.resize(Nuz);
	q1.resize(Nuz);
	q2.resize(Nuz);
}

//set the parametrs 
double mfe::Lambda(int field)									{return 2/2.75;}
double mfe::rightPart(int field, double r, double z, double t)	{return -4 * Lambda(field) + 4*t*t*t*sigma(field); }
double mfe::u_beta(double r, double z)							{return 20 * z - 27;}
double mfe::u_t(double r, double z, double t)					{ return r * r + z + t * t * t * t; }
double mfe::sigma(int field)									{return 1;}



// Build local stiffness matrix
void mfe::buildLocalG(int ielem, double r[3], double z[3], double detD)
{
	double r1 = r[0];
	double r2 = r[1];
	double r3 = r[2];

	double z1 = z[0];
	double z2 = z[1];
	double z3 = z[2];

	alfa[0] = { r2 * z3 - r3 * z2,		z2 - z3,		r3 - r2 };
	alfa[1] = { r3 * z1 - r1 * z3,		z3 - z1,		r1 - r3 };
	alfa[2] = { r1 * z2 - r2 * z1,		z1 - z2,		r2 - r1 };
	
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			alfa[i][j] /= detD;
	
	double L = Lambda(ielem);
	//double multix = L * abs(detD) * (r1 + r2 + r3) / 6.;
	double multix = L * abs(detD) /2;

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			G[i][j] = multix * (alfa[i][1] * alfa[j][1] + alfa[i][2] * alfa[j][2]); 
	
}
// Construct a local vector of the right side
void mfe::buildLocalF(int ielem, double r[3], double z[3], double detD, double t, mvector& q0, mvector& q1, mvector& q2, double diffn1, double diffn2, double diffn3)
{
	double r1 = r[0];
	double r2 = r[1];
	double r3 = r[2];

	double z1 = z[0];
	double z2 = z[1];
	double z3 = z[2];

	p_loc[0] = rightPart(ielem, r1, z1, t);
	p_loc[1] = rightPart(ielem, r2, z2, t);
	p_loc[2] = rightPart(ielem, r3, z3, t);
	double sig = sigma(ielem);
	
	c[0] = { 6 * r1 + 2 * r2 + 2 * r3,		2 * r1 + 2 * r2 + r3,		2 * r1 + r2 + 2 * r3 };
	c[1] = { 2 * r1 + 2 * r2 + r3,			2 * r1 + 6 * r2 + 2 * r3,   r1 + 2 * r2 + 2 * r3 };
	c[2] = { 2 * r1 + r2 + 2 * r3,			r1 + 2 * r2 + 2 * r3,		2 * r1 + 2 * r2 + 6 * r3 };
	
	double mult = abs(detD) / 120;
	
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			c[i][j] *= mult;

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			M[i][j] = sig * c[i][j];
		
	b_loc[0] = c[0] * p_loc - diffn1 * (M[0] * q0) - diffn2 * (M[0] * q1) - diffn3 * (M[0] * q2);
	b_loc[1] = c[1] * p_loc - diffn1 * (M[1] * q0) - diffn2 * (M[1] * q1) - diffn3 * (M[1] * q2);
	b_loc[2] = c[2] * p_loc - diffn1 * (M[2] * q0) - diffn2 * (M[2] * q1) - diffn3 * (M[2] * q2);
	
}

// Add an element to the global matrix
void mfe::addElementToGlobal(int i, int j, double elem)
{
	if (i == j)
	{
		di[i] += elem;
		return;
	}
	else
	{
		for (int ind = ig[i]; ind < ig[i + 1]; ind++)
			if (jg[ind] == j)
			{
				gg[ind] += elem;
				//gu[ind] += elem;
				return;
			}
	}
}

// Calculate the sum in front of the stiffness matrix
double mfe::calcSum(int ielem)
{
	double sum = 0.0;

	if (Nph == 2)
	{
		sum = (1 - mats[ielem].S2) / Mu[0] + mats[ielem].S2 / Mu[1];
		sum *= mats[ielem].K;
	}
	else
	{
		sum = mats[ielem].S2 / Mu[0];
		sum *= mats[ielem].K;
	}
	return sum;
}

// Assembly of the global matrix
void mfe::assemblyGlobalMatrix(double t, double t1, double t2, double t3, 
	std::vector<double>& q_0, std::vector<double>& q_1, std::vector<double>& q_2)
{
	
	mvector q0_(3, 0);
	mvector q1_(3, 0);
	mvector q2_(3, 0);
	for (int ielem = 0; ielem < Nel; ielem++)
	{
		double r1 = MeshRZ[FE[ielem][0]].first;
		double r2 = MeshRZ[FE[ielem][1]].first;
		double r3 = MeshRZ[FE[ielem][2]].first;
		double z1 = MeshRZ[FE[ielem][0]].second;
		double z2 = MeshRZ[FE[ielem][1]].second;
		double z3 = MeshRZ[FE[ielem][2]].second;

		double r[3] = { r1,r2,r3 };
		double z[3] = { z1,z2,z3 };

		double detD = (r2 - r1) * (z3 - z1) - (r3 - r1) * (z2 - z1);

		int n1 = FE[ielem][0];
		int n2 = FE[ielem][1];
		int n3 = FE[ielem][2];

		q0_[0] = q_0[n1];
		q0_[1] = q_0[n2];
		q0_[2] = q_0[n3];


		q1_[0] = q_1[n1];
		q1_[1] = q_1[n2];
		q1_[2] = q_1[n3];

		q2_[0] = q_2[n1];
		q2_[1] = q_2[n2];
		q2_[2] = q_2[n3];

		double diffn0 = (t2 * (t1 + t3 - 2 * t) + t1 * (t3 - 2 * t) + t * (3 * t - 2 * t3)) / ((t - t3) * (t - t2) * (t - t1));
		double diffn1 = (t2 * (t3 - t) + t * (t3 - 2 * t) + t * (3 * t - 2 * t3)) / ((t1 - t3) * (t1 - t2) * (t1 - t));
		double diffn2 = (t1 * (t3 - t) + t * (t3 - 2 * t) + t * (3 * t - 2 * t3)) / ((t2 - t3) * (t2 - t1) * (t2 - t));
		double diffn3 = (t2 * (t1 - t) - t1 * t + t * t) / ((t3 - t2) * (t3 - t1) * (t3 - t));
		

		buildLocalG(ielem, r, z, detD);
		buildLocalF(ielem, r, z, detD, t, q0_, q1_, q2_, diffn1, diffn2, diffn3);
	
		F[n1] += b_loc[0];
		F[n2] += b_loc[1];
		F[n3] += b_loc[2];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j <= i; j++) {
				int n1 = FE[ielem][i];
				int n2 = FE[ielem][j];
				addElementToGlobal(n1, n2, G[i][j] + diffn0 * M[i][j]);
			}
	}
}

// Construct a portrait of the global matrix
void mfe::buildPortraitOfMatrix()
{

	list[0].push_back(0);

	// Go through all FEs
	for (int ielem = 0; ielem < Nel; ielem++)
	{
		// Take the 1st basis function corresponding to the element
		for (int i = 0; i < FE[ielem].size() - 1; i++) 
			// Go through all the other functions, starting with the second
			for (int j = i + 1; j < FE[ielem].size(); j++)
			{
				// We need to add the first function to the list to all
				// functions related to the FE
				// Therefore, we determine the position where we will add (to which list)
				int insertPos = FE[ielem][j];
				// We take the element itself, which we will insert
				int element = FE[ielem][i];

				bool isIn = false;

				// Check if this element is already in the list
				for (int k = 0; k < list[insertPos].size() && !isIn; k++)
					if (element == list[insertPos][k])
						isIn = true;

				// If it is not found in the list, then add it
				if (!isIn)
					list[insertPos].push_back(element);
			}
	}

	// Sort all the resulting lists (in ascending order)
	for (int i = 0; i < Nuz; i++)
		if (!isOrdered(list[i]))
			sort(list[i].begin(), list[i].end());
	//----------------------------------------------------------------
	// Form the array ig


	// 1st and 2nd elements are always equal 1, but we will number from 0
	ig[0] = 0;
	ig[1] = 0;
	for (int i = 1; i < list.size(); i++)
		ig[i + 1] = ig[i] + list[i].size();

	//----------------------------------------------------------------
	// Form jg array
	jg.resize(ig.back());

	for (int i = 1, j = 0; i < Nuz; i++)
	{
		for (int k = 0; k < list[i].size(); k++)
			jg[j++] = list[i][k];
	}
}

// Check list for ascending order
bool mfe::isOrdered(const pvector& v)
{
	if(v.size() == 0)
		return true;
	for (int i = 0; i < v.size() - 1; i++)
		if (v[i + 1] < v[i])
			return false;
	return true;
}

// Convert matrix to dense format
void mfe::toDense(const std::string _dense)
{
	for (int i = 0; i < mat.size(); i++)
	{
		mat[i][i] = di[i];
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			mat[i][jg[j]] = gg[j];
			//mat[jg[j]][i] = gu[j];
		}
	}
	std::ofstream dense(_dense);
	dense.precision(5);
	if (dense.is_open())
	{
		for (int i = 0; i < mat.size(); i++)
		{
			dense << std::left << std::setw(20) << F[i];
			for (int j = 0; j <= i; j++)
				dense << std::left << std::setw(10) << mat[i][j];
			
			dense << std::endl << std::endl;
		}

	}
	
}

void mfe::make_bc_1(double t)
{
	bc1.clear();
	for (int i = 1; i <= z_coord.size(); i++) {

		bc1.push_back(	std::make_pair(		(i - 1) * r_coord.size(),	u_t(r_coord[0], z_coord[i - 1], t)));
		bc1.push_back(	std::make_pair(		i * r_coord.size() - 1,		u_t(r_coord.back(), z_coord[i - 1], t)));

	}

}


// Accounting for the first boundary conditions
void mfe::bc_1()
{
	for (int i = 0; i < bc1.size(); i++)
		bc1nodes[bc1[i].first] = i; // The boundary condition is set in the node

	int k;
	for (int i = 0; i < Nuz; i++)
	{
		if (bc1nodes[i] != -1)
		{
			di[i] = 1.0;
			F[i] = bc1[bc1nodes[i]].second;
			for (int j = ig[i]; j < ig[i + 1]; j++)
			{
				k = jg[j];
				if (bc1nodes[k] == -1)
					F[k] -= gg[j] * F[i];
				gg[j] = 0.0;
			}
		}
		else
		{
			for (int j = ig[i]; j < ig[i + 1]; j++)
			{
				k = jg[j];
				if (bc1nodes[k] != -1)
				{
					F[i] -= gg[j] * F[k];
					gg[j] = 0.0;
				}
			}
		}
	}
}

// Accounting for the second boundary conditions

void mfe::bc_2()
{
	for (int i = 0; i < Nbc2; i++)
	{
		int el_i = bc2[i].n_i;
		int loc_node1_i = bc2[i].loc_node1_i;
		int loc_node2_i = bc2[i].loc_node2_i;
		double Tetta_i = bc2[i].Tetta_i;

		int ind_1 = FE[el_i][loc_node1_i];
		int ind_2 = FE[el_i][loc_node2_i];
		double r1 = MeshRZ[ind_1].first;
		double r2 = MeshRZ[ind_2].first;
		double z1 = MeshRZ[ind_1].second;
		double z2 = MeshRZ[ind_2].second;
		
		double hm = sqrt((r2 - r1) * (r2 - r1) + (z2 - z1) * (z2 - z1));
		//b_2[0] = b_2[1] = hm * Tetta_i / 2;

		b_2[0] = (8 * r1 + 4 * r2) * hm * Tetta_i / 24;
		b_2[1] = (8 * r2 + 4 * r1) * hm * Tetta_i / 24 ;
		F[ind_1] += b_2[0];
		F[ind_2] += b_2[1];
	}
}
// Accounting for third edge

void mfe::bc_3()
{
	for (int i = 0; i < Nbc3; i++)
	{
		int el_i = bc3[i].n_i;
		int loc_node1_i = bc3[i].loc_node1_i;
		int loc_node2_i = bc3[i].loc_node2_i;
		double Betta_i = bc3[i].Tetta_i;

		int ind_1 = FE[el_i][loc_node1_i];
		int ind_2 = FE[el_i][loc_node2_i];
		double r1 = MeshRZ[ind_1].first;
		double r2 = MeshRZ[ind_2].first;
		double z1 = MeshRZ[ind_1].second;
		double z2 = MeshRZ[ind_2].second;
		double hm = sqrt((r2 - r1) * (r2 - r1) + (z2 - z1) * (z2 - z1));
		double ub1 = u_beta(r1, z1);
		double ub2 = u_beta(r2, z2);
		ub[0] = ub1;
		ub[1] = ub2;
		
		
		
		double k = Betta_i * hm / 24.;
		A3[0] = { k * (6 * r1 + 2 * r2), k * 2 * (r1 + r2) };
		A3[1] = { k * 2 * (r1 + r2), k * (2 * r1 + 6 * r2) };

		b_3[0] = A3[0] * ub;
		b_3[1] = A3[1] * ub;

		addElementToGlobal(ind_1, ind_1, A3[0][0]);
		addElementToGlobal(ind_1, ind_2, A3[0][1]);
		addElementToGlobal(ind_2, ind_1, A3[1][0]);
		addElementToGlobal(ind_2, ind_2, A3[1][1]);
		F[ind_1] += b_3[0];
		F[ind_2] += b_3[1];
		
	}
}



//===========================================================================================
// SLAE solution
void mfe::mult(mvector& x, mvector& y) {

	for (int i = 0; i < y.size(); i++)
		y[i] = 0;
	

	for (int i = 0; i < Nuz; i++) {
		y[i] = di[i] * x[i];
		for (int k = ig[i]; k < ig[i + 1]; k++) {
			int j = jg[k];
			y[i] += gg[k] * x[j];
			y[j] += gg[k] * x[i];
		}
	}

}


double mfe::EuclideanNorm(mvector & x) {
	double scalar = 0;
	for (int i = 0; i < Nuz; i++)
		scalar += x[i] * x[i];
	scalar = sqrt(scalar);
	return scalar;
}


// Conjugate gradient method
void mfe::MSG() {
	
	for (int i = 0; i < Nuz; i++)
		um[i] = z[i] = r[i] = 0;
	
	double scal1 = 0;
	double scal2 = 0;
	double scal3 = 0;
	double alfa = 0;
	double beta = 0;
	int k = 0;
	mult(q, um);
	for (int i = 0; i < Nuz; i++)
		r[i] = F[i] - um[i];
	z = r;
	double bnorm = EuclideanNorm(F);
	double residual = EuclideanNorm(r) / bnorm;
	if (residual > eps) {
		scal1 = 0;
		scal2 = 0;
		scal3 = 0;
		alfa = 0;
		beta = 0;
		mult(z, um);
		for (int i = 0; i < Nuz; i++) {
			scal1 += r[i] * r[i];
			scal2 += um[i] * z[i];
		}
		alfa = scal1 / scal2;
		for (int i = 0; i < Nuz; i++) {
			q[i] += alfa * z[i];
			r[i] -= alfa * um[i];
		}
		for (int i = 0; i < Nuz; i++)
			scal3 += r[i] * r[i];
		beta = scal3 / scal1;
		for (int i = 0; i < Nuz; i++)
			z[i] = r[i] + beta * z[i];
		residual = EuclideanNorm(r) / bnorm;
	}

	for (k = 1; k < maxIter && residual>eps; k++) {
		scal1 = scal3;
		scal2 = 0;
		scal3 = 0;
		alfa = 0;
		beta = 0;
		mult(z, um);
		for (int i = 0; i < Nuz; i++) {
			scal2 += um[i] * z[i];
		}
		alfa = scal1 / scal2;
		for (int i = 0; i < Nuz; i++) {
			q[i] += alfa * z[i];
			r[i] -= alfa * um[i];
		}
		for (int i = 0; i < Nuz; i++)
			scal3 += r[i] * r[i];
		beta = scal3 / scal1;
		for (int i = 0; i < Nuz; i++)
			z[i] = r[i] + beta * z[i];
		residual = EuclideanNorm(r) / bnorm;

	}
}

// Write result to file
void mfe::writeToFile(mvector& q, double t)
{
	FILE* File;
	fopen_s(&File, "data/q.txt", "a");


	if (File) {
		//fprintf(File, "%-20s%-20s%-20s%-20s%\n", "t", "q", "u", "|q-u|");
		fprintf(File, "%-20.2f%-20.5f%-20.5f%-20.10f%\n", t, q[0], u_t(r_coord[0], z_coord[0], t), abs(q[0] - u_t(r_coord[0], z_coord[0], t)));
		for (int i = 1; i < Nuz; i++) {
			int rn = i % Rsize;
			int zn = i / Rsize;
			fprintf(File, "%-20s%-20.5f%-20.5f%-20.10f%\n", " ", q[i], u_t(r_coord[rn], z_coord[zn], t), abs(q[i] - u_t(r_coord[rn], z_coord[zn], t)));
		}
		mvector u(Nuz);
		for (int i = 0; i < Nuz; i++) {
			int rn = i % Rsize;
			int zn = i / Rsize;
			u[i]= u_t(r_coord[rn], z_coord[zn], t);
		}
		mvector diff = u - q;
		double e1 = EuclideanNorm(diff);
		double e2 = EuclideanNorm(u);
		//fprintf(File, "%-20e%\n", e1/e2);
		fclose(File);
	}
	else
		std::cout << "File q.txt is not opened!";
}


void mfe::iterationProcess() {

	for (int i = 0; i < Nuz; i++) {
		int rn = i % Rsize;
		int zn = i / Rsize;
		q2[i] = u_t(r_coord[rn],z_coord[zn], time[0]);//U(t0)
		q1[i] = u_t(r_coord[rn], z_coord[zn], time[1]);//U(t1)
		q0[i] = u_t(r_coord[rn], z_coord[zn], time[2]);//U(t2)
	}
	gg.resize(ig.back(), 0);
	
	for (int s = 3; s < n; s++)
	{
		for (int i = 0; i < gg.size(); i++)
			gg[i] = 0;
		for (int i = 0; i < Nuz; i++)
			di[i] = F[i] = 0;
		double t = time[s]; //tj
		double t1 = time[s-1]; //tj-1
		double t2 = time[s-2]; //tj-2
		double t3 = time[s-3]; //tj-3

		// q0 = q(s-1)
		// q = q(s)
		assemblyGlobalMatrix(t, t1, t2, t3, q0, q1, q2);
		toDense("data/matrix.txt");
		make_bc_1(t);
		bc_1();
		toDense("data/matrix1.txt");
		MSG();
		q2 = q1;
		q1 = q0;
		q0 = q;
		writeToFile(q, t);
	}
}