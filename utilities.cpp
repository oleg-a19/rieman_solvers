#include <vector>
#include "utilities.h"

std::vector<double> w2u (const std::vector<double>& w, const double gamma) {
	std::vector<double> u;
	u.push_back(w[0]);
	u.push_back(w[0]*w[1]);
	u.push_back(w[2]/(gamma-1)+0.5*w[0]*w[1]*w[1]);
	return u;
}

std::vector<double> u2w (const std::vector<double>& u, const double gamma) {
	std::vector<double> w;
	w.push_back(u[0]);
	w.push_back(u[1]/u[0]);
	w.push_back((gamma-1)*(u[2] - 0.5*w[0]*w[1]*w[1]));
	return w;
}