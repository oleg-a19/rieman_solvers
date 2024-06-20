#ifndef CELL_H
#define CELL_H

#include <vector>

class cell {
friend class grid;
public:
	cell();
	double xc;
	std::vector<double> u;
	std::vector<double> w;
	std::vector<double> res;
	std::vector<double> dw;
	double gamma;
	
private:
};

#endif
