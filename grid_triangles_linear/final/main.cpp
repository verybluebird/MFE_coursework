#include <iostream>
#include <fstream>
#include "grid.h"

int main()
{
	Grid grid;
	grid.input();
	grid.nodes();
	grid.elems();
	grid.material();
	grid.gr_bc1();
	grid.gr_bc2();
	grid.time();
	grid.out_coords();
	//grid.print_profile();

	return 0;
}