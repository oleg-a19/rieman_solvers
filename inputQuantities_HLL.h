#ifndef INPUT_QUANTITIES_H
#define INPUT_QUANTITIES_H
#include <iostream>
#include <string>


class inputQuantities {
public:
	inputQuantities() = default;
	inputQuantities(const std::string fileName);
	void printInputData();
	double xmin, xmax, tf, cfl, ncells;
	double rho1, u1, p1;
	double rho2, u2, p2;

	int get_task_size() {
		return task_size;
	}

private:
	int task_size;
};

#endif
