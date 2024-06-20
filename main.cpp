#include <iostream>
#include <string>

#include "grid.h"
#include "solver.h"


int main() {
	grid mesh("input_HLL.txt");
	mesh.setInitialCondition();

	std::string output_file1 = "data_weno2_n640.csv";
	std::string output_file2 = "data_HLLC_minmod_n640.csv";
	
	solver euler;
	euler.solve_weno(mesh, output_file1);
	//euler.solve(mesh, output_file2, "HLLC_minmod");
	
	
	return 0;
}
