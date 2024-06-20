#include "inputQuantities_HLL.h"
#include <iostream>
#include <string>
#include <fstream>


inputQuantities::inputQuantities (const std::string fileName) {
	std::ifstream in(fileName);
	in >> xmin;
	in >> xmax;
	in >> tf;
	in >> cfl;
	in >> ncells;
	in >> rho1;
	in >> u1;
	in >> p1;
	in >> rho2;
	in >> u2;
	in >> p2;
	in >> task_size;
	in.close();
}
void inputQuantities::printInputData () {
	std::cout << xmin << " " << xmax << " " << tf << " " << cfl << " " << ncells << '\n'
	 << rho1 << " " << u1 << " " << p1 << '\n' << rho2 << " " << u2 << " " << p2 << '\n'
	 << task_size << '\n';
}
