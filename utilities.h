#ifndef UTIL_H
#define UTIL_H

#include <vector>

std::vector<double> w2u (const std::vector<double>& w, const double gamma);

std::vector<double> u2w (const std::vector<double>& u, const double gamma);

#endif