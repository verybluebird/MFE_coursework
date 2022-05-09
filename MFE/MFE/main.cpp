#include "mfe.h"
#include <iostream>

using namespace std;

int main()
{
	try
	{
		mfe m;

		m.buildPortraitOfMatrix();
		//m.assemblyGlobalMatrix();
		//m.toDense("matrix2.txt");
		////m.bc_2();
		////m.bc_3();
		////m.toDense("matrix2.txt");
		//m.bc_1();
		//m.toDense("matrix2.txt");
		//m.MSG();
		//m.writeToFile(m.q);
		m.iterationProcess();
	}
	catch (int error)
	{
		switch (error)
		{
		case 0:
			cout << "Unable to read file!" << endl;
			break;

		case 1:
			cout << "Point is out of range!" << endl;
			break;

		case 2:
			cout << "Unable to write result to file!" << endl;
			break;
		}		
	}
	
	return 0;
}