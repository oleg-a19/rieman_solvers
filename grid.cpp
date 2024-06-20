#include "grid.h"
#include "utilities.h"

grid::grid (const std::string fileName) : quantities(fileName) {
}

void grid::setInitialCondition () {
    quantities.printInputData();
    cells.resize(quantities.ncells+2);
    dx = (quantities.xmax-quantities.xmin)/quantities.ncells;

    for (int i=0; i<quantities.ncells+2; ++i) {
		cells[i].w.resize(quantities.get_task_size());
        cells[i].gamma = 1.4;
		if (i <= quantities.ncells/2) {
            cells[i].w[0] = quantities.rho1;
            cells[i].w[1] = quantities.u1;
			cells[i].w[2] = quantities.p1;			
		} else {
            cells[i].w[0] = quantities.rho2;
            cells[i].w[1] = quantities.u2;
			cells[i].w[2] = quantities.p2;
		}
		cells[i].u = w2u(cells[i].w, cells[i].gamma);	
		cells[i].xc = quantities.xmin + (i-1)*dx;	
	}
}

inputQuantities grid::get_input_quantities() {
	return quantities;
}

std::vector<cell> grid::get_cells() {
	return cells;
}

double grid::get_dx() {
	return dx;
}