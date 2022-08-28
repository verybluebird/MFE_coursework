#include "include/mfe.h"


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
			std::cout << "Unable to read file!" << std::endl;
			break;

		case 1:
			std::cout << "Point is out of range!" << std::endl;
			break;

		case 2:
			std::cout << "Unable to write result to file!" << std::endl;
			break;
		}		
	}
	
	return 0;
}