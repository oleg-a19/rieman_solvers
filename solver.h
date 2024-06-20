#ifndef SOLVER_H
#define SOLVER_H

#include <string>
#include <fstream>
#include <vector>
#include "cell.h"
#include "grid.h"

class solver {
public:
    solver ();
    std::vector<double> hll_flux (const std::vector<double>& uL,
                                  const std::vector<double>& uR,
                                  const double gamma);
    
    std::vector<double> hllc_flux (const std::vector<double>& uL,
                                  const std::vector<double>& uR,
                                  const double gamma);

    double minmod(const double a, const double b);

    double timestep (double cfl, double dx, double gamma,
                     int ncells, const std::vector<cell>& cells);

    void solve (grid& mesh, std::string out, std::string scheme);

    void solve_weno (grid& mesh, std::string out);

    std::vector<std::vector<double>> calc_weight(const int node, const std::string side,
                       const std::vector<cell>& cells, const int task_size);

    std::vector<std::vector<double>> calc_face_value (const int node, const std::string side,
                       const std::vector<cell>& cells, const int task_size);

private:
    double t;
    int nsteps;
};

#endif