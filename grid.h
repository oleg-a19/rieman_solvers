#ifndef GRID_H
#define GRID_H

#include "inputQuantities_HLL.h"
#include "cell.h"
#include <vector>
#include <string>

class grid {
public:
	grid() = default;
	grid(const std::string fileName);

	void setInitialCondition();
	inputQuantities get_input_quantities();
	std::vector<cell> get_cells();
	double get_dx();

private:
	inputQuantities quantities;
	std::vector<cell> cells;
	double dx;
};

#endif
