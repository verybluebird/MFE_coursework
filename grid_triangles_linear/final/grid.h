#pragma once
#include <vector>
using namespace std;
class Grid
{
private:
	/*параметры расчетной области*/
	double Rw; //радиус скважины
	double Rb; //расстояние до удаленной границы (“бак”)
	int NL; //количество слоев
	vector<double> H; // толщина слоев
	vector<double> K; //структурная проницаемость слоев
	vector<double> Phi; // пористость слоев
	vector<double> S2; //нефтенасыщенность слоев
	int nr; //кол-во разбиений по r
	double kr; //коэф. разрядки по r
	vector<int> nz; //кол-во разбиений по z в каждом слое
	int M; // вложение сетки
	int Nzp; //количество зон перфораций
	vector<double> Pu; //верхняя граница з.п.
	vector<double> Pd; //нижняя граница з.п.
	vector<double> Tetta; //мощность
	int Nph; // количество фаз
	vector<double> Mu; //вязкость фаз
	double Plast; //пластовое давление
	vector<vector<int>> global_numbers; //глобальные номера узлов
	double z0 = 0; //откуда начинается первый сверху слой по z
	vector<double> r_coord; // координаты по r
	vector<double> z_coord; // координаты по z

public:
	void input();// ввод  данных
	int Nuz; //количество узлов
	vector<vector<double>> coords;//{r_i y_i}
	int Nel; //
	vector<vector<int>> numbers; // глобальные номера
	vector<vector<double>> materials;
	int Nbc1; //кол-во первых краевых условий
	vector<double> bc1;
	int Nbc2;//кол-во вторых краевых условий
	vector<vector<double>> bc2;
	void nodes();
	void elems();
	void material();
	void gr_bc1();
	void gr_bc2();
	void nested_grid(vector<double>& coord); //увеличивает кол-во узлов вдвое
	void add_if_not_exist_and_sort(double L); // добавить новый узел в сетку, если его не существует
	void print_profile(); // распечатать профиль матрицы


};