#include "mfe.h"

//===========================================================================================
// Конструктор: читаем сетку, КЭ с номерами их узлов, 1-ые и 2-ые краевые, материалы
// и свойства фаз, параметры для решения СЛАУ. Задаем точки и веса для интегрирования
mfe::mfe()
{
	//-----------------------------------------------------
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
	else throw BAD_READ;

	//-----------------------------------------------------
	// Конечный элемент задается номерами его узлов
	std::ifstream inElem("elem.txt");
	if (inElem.is_open())
	{
		inElem >> Nel;
		FE.resize(Nel);
		for (size_t i = 0; i < Nel; i++)
			FE[i].resize(3);

		size_t n1, n2, n3;
		for (size_t i = 0; i < Nel; i++)
		{
			inElem >> n1 >> n2 >> n3;
			FE[i] = { n1, n2, n3 };
		}

		inElem.close();
	}
	else throw BAD_READ;

	//-----------------------------------------------------
	std::ifstream inMat("mat.txt");
	if (inMat.is_open())
	{
		mats.resize(Nel);

		double K, Phi, S2;
		for (size_t i = 0; i < Nel; i++)
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
	std::ifstream inPhase("phaseprop.txt");
	if (inPhase.is_open())
	{
		inPhase >> Nph;
		Mu.resize(Nph);
		
		for (size_t i = 0; i < Nph; i++)
			inPhase >> Mu[i];

		inPhase.close();
	}
	else throw BAD_READ;


	//-----------------------------------------------------
	std::ifstream bc1File("bc1.txt");
	if (bc1File.is_open())
	{
		bc1File >> Nbc1;
		bc1.resize(Nbc1);
		for (size_t i = 0; i < Nbc1; i++)
		{
			size_t n_i;
			double value_bc1;

			bc1File >> n_i >> value_bc1;
			bc1[i] = std::make_pair(n_i, value_bc1);
		}
		bc1File.close();
	}
	else throw BAD_READ;

	//-----------------------------------------------------
	std::ifstream bc2File("bc2.txt");
	if (bc2File.is_open())
	{
		bc2File >> Nbc2;
		bc2.resize(Nbc2);

		size_t el_i, loc_node1_i, loc_node2_i;
		double Tetta_i;

		for (size_t i = 0; i < Nbc2; i++)
		{
			bc2File >> el_i >> loc_node1_i >> loc_node2_i >> Tetta_i;
			bc2[i] = { el_i, loc_node1_i, loc_node2_i, Tetta_i };
		}
		bc2File.close();
	}
	else throw BAD_READ;
	//-----------------------------------------------------
	std::ifstream bc3File("bc3.txt");
	if (bc3File.is_open())
	{
		bc3File >> Nbc3;
		bc3.resize(Nbc3);

		size_t el_i, loc_node1_i, loc_node2_i;
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
	std::ifstream slau("kuslau.txt");
	if (slau.is_open())
	{
		slau >> maxIter >> eps;
		slau.close();
	}
	else throw BAD_READ;

	b_loc.resize(3);
	G.resize(3);
	for (size_t i = 0; i < 3; i++)
		G[i].resize(3);

	alfa.resize(3);
	for (size_t i = 0; i < 3; i++)
		alfa[i].resize(3);
	lambda = 1;
	F.resize(Nuz);
}
double Lambda(int field)
{
	if (field == 0)
		return 10;
	else 
	{
		return 1;
	}	
}
double rightPart(int field)
{
	if (field == 0)
		return 20;
	else
	{
		return 0;
	}

}
double u_beta(double r, double z)
{
	return 20 * z - 27;

}
// Построить локальную матрицу жесткости*
void mfe::buildLocalG(size_t ielem)
{
	double r1 = MeshRZ[FE[ielem][0]].first;
	double r2 = MeshRZ[FE[ielem][1]].first;
	double r3 = MeshRZ[FE[ielem][2]].first;
	double z1 = MeshRZ[FE[ielem][0]].second;
	double z2 = MeshRZ[FE[ielem][1]].second;
	double z3 = MeshRZ[FE[ielem][2]].second;
	double detD = (r2 - r1) * (z3 - z1) - (r3 - r1) * (z2 - z1);
	
	alfa[0] = { r2 * z3 - r3 * z2, z2 - z3, r3 - r2 };
	alfa[1] = { r3 * z1 - r1 * z3, z3 - z1, r1 - r3 };
	alfa[2] = { r1 * z2 - r2 * z1, z1 - z2, r2 - r1 };
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
// Построить локальный вектор правой части*
void mfe::buildLocalF(size_t ielem)
{
	double r1 = MeshRZ[FE[ielem][0]].first;
	double r2 = MeshRZ[FE[ielem][1]].first;
	double r3 = MeshRZ[FE[ielem][2]].first;
	double z1 = MeshRZ[FE[ielem][0]].second;
	double z2 = MeshRZ[FE[ielem][1]].second;
	double z3 = MeshRZ[FE[ielem][2]].second;
	double detD = (r2 - r1) * (z3 - z1) - (r3 - r1) * (z2 - z1);
	mvector p;
	p.resize(3);
	for (int i = 0; i < 3; i++)
	{
		p[i] = rightPart(ielem);
	}

	matrix c;
	c.resize(3);
	for (int i = 0; i < 3; i++)
	{
		c.resize(3);
	}
	c[0] = { 2,1,1 };
	c[1] = { 1,2,1 };
	c[2] = { 1,1,2 };
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			c[i][j] *= abs(detD) / 24;
	for (int i = 0; i < 3; i++)
	{
		b_loc[i] = c[0] * p;
	}
}

// Добавить элемент в глобальную матрицу
void mfe::addElementToGlobal(size_t i, size_t j, double elem)
{
	if (i == j)
	{
		di[i] += elem;
		return;
	}
	
	else
	{
		for (size_t ind = ig[i]; ind < ig[i + 1]; ind++)
			if (jg[ind] == j)
			{
				gg[ind] += elem;
				gu[ind] += elem;
				return;
			}
	}
}

// Посчитать сумму, которая стоит перед матрицей жесткости
double mfe::calcSum(size_t ielem)
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

// Сборка глобальной матрицы
void mfe::assemblyGlobalMatrix()
{
	di.resize(Nuz);
	gg.resize(ig.back());
	gu.resize(ig.back());

	for (size_t ielem = 0; ielem < Nel; ielem++)
	{
		buildLocalG(ielem);
		buildLocalF(ielem);

		//// Домножаем на коэффициент
		//double sum = calcSum(ielem);
		//for (size_t i = 0; i < G.size(); i++)
		//for (size_t j = 0; j < G[i].size(); j++)
		//		G[i][j] *= sum;
		int n1 = FE[ielem][0];
		int n2 = FE[ielem][1];
		int n3 = FE[ielem][2];
		F[n1] += b_loc[0];
		F[n2] += b_loc[1];
		F[n3] += b_loc[2];
		for (size_t i = 0; i < FE[ielem].size(); i++)
			for (size_t j = 0; j <= i; j++) {
				int n1 = FE[ielem][i];
				int n2 = FE[ielem][j];
				addElementToGlobal(n1, n2, G[i][j]);
				
			}
 		toDense("matrix2.txt");
		
	}
	G.clear();
}

// Построить портрет глобальной матрицы*правильно
void mfe::buildPortraitOfMatrix()
{
	std::vector<std::vector<size_t>> list;

	list.resize(Nuz);

	list[0].push_back(0);

	// Идем по всем КЭ
	for (size_t ielem = 0; ielem < Nel; ielem++)
	{
		// Берем 1-ую соответствующую элементу базисную функцию
		for (size_t i = 0; i < FE[ielem].size() - 1; i++) 
			// Идем по всем остальным функциям, начиная со второй 
			for (size_t j = i + 1; j < FE[ielem].size(); j++)
			{
				// Нужно добавить первую функцию(меньшую) в список ко всем 
				// функциям, относящимся к КЭ
				// Поэтому определяем позицию, куда будем добавлять (в какой список)
				size_t insertPos = FE[ielem][j];
				// Берем сам элемент, который будем вставлять
				size_t element = FE[ielem][i];

				bool isIn = false;

				// Проверим, есть ли уже этот элемент в списке
				for (size_t k = 0; k < list[insertPos].size() && !isIn; k++)
					if (element == list[insertPos][k])
						isIn = true;

				// Если он в списке не найден, то добавляем его
				if (!isIn)
					list[insertPos].push_back(element);
			}
	}

	// Сортируем все получившиеся списки (по возрастанию номеров)
	for (size_t i = 0; i < Nuz; i++)
		if (!isOrdered(list[i]))
			sort(list[i].begin(), list[i].end());
	//----------------------------------------------------------------
	// Формируем массив ig
	ig.resize(Nuz + 1);

	// 1-ый и 2-ой элементы всегда равны 1, но мы будем нумеровать с 0
	ig[0] = 0;
	ig[1] = 0;
	for (size_t i = 1; i < list.size(); i++)
		ig[i + 1] = ig[i] + list[i].size();

	//----------------------------------------------------------------
	// Формируем массив jg
	jg.resize(ig.back());

	for (size_t i = 1, j = 0; i < Nuz; i++)
	{
		for (size_t k = 0; k < list[i].size(); k++)
			jg[j++] = list[i][k];
	}
}

// Проверка списка на упорядоченность по возрастанию
bool mfe::isOrdered(const pvector& v)
{
	if(v.size() == 0)
		return true;
	for (size_t i = 0; i < v.size() - 1; i++)
		if (v[i + 1] < v[i])
			return false;
	return true;
}

// Перевод матрицы в плотный формат
void mfe::toDense(const std::string _dense)
{
	mat.resize(Nuz);
	for (size_t i = 0; i < mat.size(); i++)
		mat[i].resize(Nuz, 0);
	for (size_t i = 0; i < mat.size(); i++)
	{
		mat[i][i] = di[i];
		for (size_t j = ig[i]; j < ig[i + 1]; j++)
		{
			mat[i][jg[j]] = gg[j];
			mat[jg[j]][i] = gu[j];
		}
	}
	std::ofstream dense(_dense);
	dense.precision(5);
	if (dense.is_open())
	{
		for (size_t i = 0; i < mat.size(); i++)
		{
			
			for (size_t j = 0; j < mat.size(); j++)
				dense << std::left << std::setw(10) << mat[i][j];
			dense << std::left << std::setw(20) << F[i];
			dense << std::endl << std::endl;
		}

	}
	mat.clear();
}




// Учет первых краевых
void mfe::bc_1()
{
	std::vector<int> bc1nodes(Nuz, -1);
	for (size_t i = 0; i < bc1.size(); i++)
		bc1nodes[bc1[i].first] = i; // В узле задано краевое

	size_t k;
	for (size_t i = 0; i < Nuz; i++)
	{
		if (bc1nodes[i] != -1)
		{
			di[i] = 1;
			F[i] = bc1[bc1nodes[i]].second;
			for (size_t j = ig[i]; j < ig[i + 1]; j++)
			{
				gg[j] = 0.0;
			}
			for (size_t j = 0; j < gu.size(); j++)
			{
				if (jg[j] == i)
					gu[j] = 0;
			}
		}
	}
}

// Учет вторых краевых

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
		std::vector<double> b_2(2);
		double hm = sqrt((r2 - r1) * (r2 - r1) + (z2 - z1) * (z2 - z1));
		b_2[0] = b_2[1] = hm * Tetta_i / 2;
		/*b_2[0] = (8 * r1 + 4 * r2) * hm * Tetta_i / 24;
		b_2[1] = (8 * r2 + 4 * r1) * hm * Tetta_i/24 ;*/
		F[ind_1] += b_2[0];
		F[ind_2] += b_2[1];
	}
}
// Учет третьих краевых

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
		std::vector<double> b_3(2);
		double hm = sqrt((r2 - r1) * (r2 - r1) + (z2 - z1) * (z2 - z1));
		double ub1 = u_beta(r1, z1);
		double ub2 = u_beta(r2, z2);
		b_3[0] = Betta_i * hm * (2 * ub1 + ub2)/6;
		b_3[1] = Betta_i * hm * (ub1 + 2 * ub2) / 6;
		F[ind_1] += b_3[0];
		F[ind_2] += b_3[1];
		
		matrix A3;
		A3.resize(2);
		A3[0].resize(2);
		A3[1].resize(2);
		double k = Betta_i * hm / 6;
		A3[0] = { 2 * k, k };
		A3[1] = { k, 2 * k };

		addElementToGlobal(ind_1, ind_1, A3[0][0]);
		addElementToGlobal(ind_1, ind_2, A3[0][1]);
		addElementToGlobal(ind_2, ind_1, A3[1][0]);
		addElementToGlobal(ind_2, ind_2, A3[1][1]);
		
	}
}



//===========================================================================================
// Решение СЛАУ
void mfe::mult(mvector& x, mvector& y) {
	for (int i = 0; i < y.size(); i++)
	{
		y[i] = 0;
	}

	for (int i = 0; i < Nuz; i++) {
		y[i] = di[i] * x[i];
		for (int k = ig[i]; k < ig[i + 1]; k++) {
			int j = jg[k];
			y[i] += gg[k] * x[j];
			y[j] += gu[k] * x[i];
		}
	}

}

void mfe::LOS()
{
	double alpha, beta;
	std::vector<double> Ar;
	std::vector<double> Aq;
	q.resize(di.size(), 0);
	r.resize(di.size(), 0);
	p.resize(di.size(), 0);
	Ar.resize(di.size(), 0);
	Aq.resize(di.size(), 0);
	mult(q, Aq);
	r = F - Aq;
	z = r;
	mult(z, p);

	double rr = (r * r);
	for (size_t k = 0; k < maxIter && rr >= eps; k++)
	{
		double pp = (p * p);
		alpha = (p * r) / pp;
		q = q + alpha * z;
		rr = (r * r) - alpha * alpha * pp;
		r = r - alpha * p;
		mult(r, Ar);
		beta = -(p * Ar) / pp;
		z = r + beta * z;
		p = Ar + beta * p;
	}

	r.clear();
	p.clear();
	z.clear();
}



double mfe::EuclideanNorm(mvector & x) {
	double scalar = 0;
	for (int i = 0; i < Nuz; i++)
		scalar += x[i] * x[i];
	scalar = sqrt(scalar);
	return scalar;
}


// Метод сопряженных градиентов
void mfe::MSG() {
	mvector um;
	mvector r;
	
	um.resize(Nuz);
	r.resize(Nuz);
	q.resize(Nuz);
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

// Записать результат в файл
void mfe::writeToFile(mvector q)
{
	std::ofstream res("q.txt");
	if (res.is_open())
	{
		res.precision(15);
		res << q.size() << std::endl << std::endl;
		for (size_t i = 0; i < q.size(); i++)
			res << std::setprecision(15) << q[i] << std::endl;
		res.close();
		toDense("matrix.txt");
	}
	else
		throw BAD_WRITE;
}