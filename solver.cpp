#include <cmath>
#include <algorithm>
#include <iomanip>

#include "solver.h"
#include "utilities.h"


solver::solver() {
    t = 0;
    nsteps = 0;
}

std::vector<double> solver::hll_flux(const std::vector<double>& uL,
                 const std::vector<double>& uR, const double gamma) {
    std::vector<double> flux(3);
	std::vector<double> flux_L(3);
	std::vector<double> flux_R(3);
	std::vector<double> wL = u2w(uL, gamma);
	std::vector<double> wR = u2w(uR, gamma);


	double a_L = std::sqrt(gamma*wL[2]/wL[0]);
	double a_R = std::sqrt(gamma*wR[2]/wR[0]);
	
	double u_tilda = (std::sqrt(wL[0])*wL[1] + std::sqrt(wR[0])*wR[1])
					 /
			  		 (std::sqrt(wL[0]) + std::sqrt(wR[0]));
			  		 
	double H_L = a_L*a_L / (gamma-1)
				 +
				 wL[1]*wL[1] / 2;
	
	double H_R = a_R*a_R / (gamma-1)
				 +
				 wR[1]*wR[1] / 2;

	double H_tilda = (std::sqrt(wL[0])*H_L + std::sqrt(wR[0])*H_R)
					 /
					 (std::sqrt(wL[0]) + std::sqrt(wR[0]));			 

	double a_tilda = std::sqrt( (gamma-1)*(H_tilda - 0.5*u_tilda*u_tilda) );
	
	double sL = u_tilda - a_tilda;
	double sR = u_tilda + a_tilda;

	flux_L[0] = wL[0]*wL[1];
	flux_L[1] = wL[0]*wL[1]*wL[1] + wL[2];
	flux_L[2] = wL[0]*wL[1]*H_L;
	
	flux_R[0] = wR[0]*wR[1];
	flux_R[1] = wR[0]*wR[1]*wR[1] + wR[2];
	flux_R[2] = wR[0]*wR[1]*H_R;

	if (sL >= 0) {
		return flux_L;
	} else if (sR <= 0) {
		return flux_R;
	} else {
		for (int i=0; i<3; ++i) {
			flux[i] = ( sR*flux_L[i] - sL*flux_R[i]
					  +
					  sL*sR*(uR[i]-uL[i]) )
					  /
					  (sR - sL);
		}
	}
    return flux;
}


std::vector<double> solver::hllc_flux (const std::vector<double>& uL,
                                  const std::vector<double>& uR,
                                  const double gamma) { // требуются: w, u, rho, T, m, gamma
	std::vector<double> wL = u2w(uL, gamma);
	std::vector<double> wR = u2w(uR, gamma);
	int dim = 3;

	double a_L = std::sqrt(gamma*wL[2]/wL[0]);
	double a_R = std::sqrt(gamma*wR[2]/wR[0]);
	
	double u_tilda = (std::sqrt(wL[0])*wL[1] + std::sqrt(wR[0])*wR[1])
					 /
			  		 (std::sqrt(wL[0]) + std::sqrt(wR[0]));
			  		 
	double H_L = a_L*a_L / (gamma-1)
				 +
				 wL[1]*wL[1] / 2;
	
	double H_R = a_R*a_R / (gamma-1)
				 +
				 wR[1]*wR[1] / 2;

	double H_tilda = (std::sqrt(wL[0])*H_L + std::sqrt(wR[0])*H_R)
					 /
					 (std::sqrt(wL[0]) + std::sqrt(wR[0]));		 

	double a_tilda = std::sqrt( (gamma-1)*(H_tilda - 0.5*u_tilda*u_tilda) );
	
	// compute fastly and slowly wave speed
	double sL = u_tilda - a_tilda;
	double sR = u_tilda + a_tilda;
	
	// comput fields in star region
	double s_star = (wR[2]-wL[2] + wL[0]*wL[1]*(sL-wL[1]) - wR[0]*wR[1]*(sR-wR[1]))
					/
					(wL[0]*(sL-wL[1]) - wR[0]*(sR-wR[1])); // name like Toro
	
	std::vector<double> uL_star(dim);
	double rho_L_star = wL[0] * (sL-wL[1])/(sL-s_star);
	uL_star[0] = uL[0] * (sL-wL[1])/(sL-s_star);
	uL_star[1] = rho_L_star * s_star;
	uL_star[2] = rho_L_star * (uL[2]/wL[0] + (s_star - wL[1])
														 *(s_star + wL[2] / (wL[0]*(sL-wL[1])) )
														 /wL[0]);
	
	std::vector<double> uR_star(dim);
	double rho_R_star = wR[0] * (sR-wR[1])/(sR-s_star);
	uR_star[0] = uR[0] * (sR-wR[1])/(sR-s_star);
	uR_star[1] = rho_R_star * s_star;
	uR_star[2] = rho_R_star * (uR[2]/wR[0] + (s_star - wR[1])
														 *(s_star + wR[2] / (wR[0]*(sR-wR[1])) )
														 /wR[0]);
	
	std::vector<double> flux_L(dim);
	flux_L[0] = wL[0]*wL[1];
	flux_L[1] = wL[0]*wL[1]*wL[1] + wL[2];
	flux_L[2] = wL[0]*wL[1]*H_L;
	
	std::vector<double> flux_R(dim);
	flux_R[0] = wR[0]*wR[1];
	flux_R[1] = wR[0]*wR[1]*wR[1] + wR[2];
	flux_R[2] = wR[0]*wR[1]*H_R;

	std::vector<double> flux_L_star(dim);
	for (unsigned int i=0; i<dim; ++i) flux_L_star[i] = flux_L[i] + sL*(uL_star[i] - uL[i]);
	
	std::vector<double> flux_R_star(dim); 
	for (unsigned int i=0; i<dim; ++i) flux_R_star[i] = flux_R[i] + sR*(uR_star[i] - uR[i]);
	
	std::vector<double> flux(dim);
	if (sL >= 0) for (unsigned int i=0; i<dim; ++i) flux[i] = flux_L[i];
	else if (sL < 0 && 0 <= s_star) for (unsigned int i=0; i<dim; ++i) flux[i] = flux_L_star[i];
	else if (s_star < 0 && 0 < sR) for (unsigned int i=0; i<dim; ++i) flux[i] = flux_R_star[i];
	else for (unsigned int i=0; i<dim; ++i) flux[i] = flux_R[i];
	
	return flux;
}

double solver::minmod(const double a, const double b) {
	double minmod = 1e5;
	if (b > 0) {
		minmod = std::max(0., std::max(std::min(a,b), std::min(a,b)));
	} else {
		minmod = std::min(0., std::min(std::max(a,b), std::max(a,b)));
	}
	//std::cout << a << "               "<< b <<"               " << minmod << '\n';
	return minmod;
}


double solver::timestep(double cfl, double dx, double gamma,
                        int ncells, const std::vector<cell>& cells) {
    double dt;
	double u, c, max_speed;
	
	max_speed = -1.;
	
	for (int i=0; i<ncells+2; ++i) {
		u = cells[i].w[1];
		c = std::sqrt(gamma*cells[i].w[2]/cells[i].w[0]);
		
		max_speed = std::max(max_speed, std::abs(u)+c);
	}
	
	dt = cfl*dx/max_speed;                   
    return dt;
}


void solver::solve(grid& mesh, std::string out, std::string scheme) {
	std::vector<cell> cells = mesh.get_cells();
	double cfl = mesh.get_input_quantities().cfl;
	double gamma = cells[0].gamma;
	int ncells = mesh.get_input_quantities().ncells;
	double tf = mesh.get_input_quantities().tf;
	double dx = mesh.get_dx();

    int task_size = mesh.get_input_quantities().get_task_size();

	std::vector<double> flux(task_size);
    for (auto& elem : cells) elem.res.resize(task_size);

	for (int itime=1; itime<1e6; ++itime) {
		if (t == tf) break;
		double dt = timestep(cfl, dx, gamma, ncells, cells);
		if (t+dt > tf) dt = tf - t;
		t += dt;
		nsteps += 1;

		// compute difference and apply slope limiter
		for (int i=1; i<ncells+1; ++i) {
			std::vector<double> dwl(task_size);
			std::vector<double> dwr(task_size);
			cells[i].dw.resize(task_size);
			for (unsigned int j=0; j<task_size; ++j) {
				// only for uniform grid !!!
				dwl[j] = cells[i].w[j]-cells[i-1].w[j];
				dwr[j] = cells[i+1].w[j]-cells[i].w[j];
				cells[i].dw[j] = minmod(dwl[j], dwr[j]);
			}
		}

		//обнуляем невязку на каждом шаге по времени в каждой ячейке
		for (int j=1; j<ncells+1; ++j) std::fill(cells[j].res.begin(), cells[j].res.end(), 0.);
		
		std::vector<double> wL(task_size);
		std::vector<double> wR(task_size);
		for (int i=1; i<ncells; ++i) {
			for (unsigned int j=0; j<task_size; ++j) {
				wL[j] = cells[i].w[j] + 0.5*cells[i].dw[j];
				wR[j] = cells[i+1].w[j] - 0.5*cells[i+1].dw[j];
			}
			std::vector<double> uL = w2u(wL, gamma);
			std::vector<double> uR = w2u(wR, gamma);
			if (scheme == "HLL") {
				flux = hll_flux(cells[i].u, cells[i+1].u, gamma); //find flux
			} else if (scheme == "HLLC") {
				flux = hllc_flux(cells[i].u, cells[i+1].u, gamma); //find flux
			} else if (scheme == "HLL_minmod") {
				flux = hll_flux(uL, uR, gamma);
			} else if (scheme == "HLLC_minmod") flux = hllc_flux(uL, uR, gamma);

			// residual upadate now
			for (int j=0; j<3; ++j) {
				cells[i].res[j] += flux[j];
				cells[i+1].res[j] -= flux[j];
			}
		}

		// Left most face: left face of cell i=1
		if (scheme == "HLL") {
			flux = hll_flux(cells[2].u, cells[2].u, gamma);
		} else if (scheme == "HLLC") {
			flux = hllc_flux(cells[2].u, cells[2].u, gamma);
		} else if (scheme == "HLLC_minmod") {
			for (unsigned int j=0; j<task_size; ++j) {
				wR[j] = cells[2].w[j] - 0.5*cells[2].dw[j];
				wL[j] = wR[j];
			}
			std::vector<double> uL = w2u(wL, gamma);
			std::vector<double> uR = w2u(wR, gamma);
			flux = hllc_flux(uL, uR, gamma);
		}
		for (int i=0; i<3; ++i) cells[1].res[i] -= flux[i];
		
		// Right most face: right face of cell i=ncells.
		
		if (scheme == "HLL") {
			flux = hll_flux(cells[ncells].u, cells[ncells].u, gamma);
		} else if (scheme == "HLLC") {
			flux = hllc_flux(cells[ncells].u, cells[ncells].u, gamma);
		} else if (scheme == "HLLC_minmod") {
			for (unsigned int j=0; j<task_size; ++j) {
				wL[j] = cells[ncells].w[j] + 0.5*cells[ncells].dw[j];
				wR[j] = wL[j];
			}
			std::vector<double> uL = w2u(wL, gamma);
			std::vector<double> uR = w2u(wR, gamma);
			flux = hllc_flux(uL, uR, gamma);
		}
		for (int i=0; i<3; ++i) cells[ncells].res[i] += flux[i];
		
		// Solution update
		for (int i=1; i<ncells+1; ++i) {
			for (int j=0; j<3; ++j) cells[i].u[j] -= dt/dx * cells[i].res[j];
			cells[i].w = u2w(cells[i].u, gamma);
			
		}
		// Copy the solutions to the ghost cells.
		for (int i=0; i<3; ++i) {
			cells[0].w[i] = cells[1].w[i];
			cells[ncells+1].w[i] = cells[ncells].w[i];
		}
	}

	std::ofstream f(out);
	for (int i=1; i<ncells+1; ++i) f << cells[i].xc << ", " << cells[i].w[0] << ", " 
									 << cells[i].w[1] << ", " << cells[i].w[2] << '\n';
	f.close();
}


std::vector<std::vector<double>> solver::calc_weight(const int node, const std::string side,
                   const std::vector<cell>& cells, const int task_size) {

	int n = cells.size()-2;

	if (node == 1 && side == "L") {
		std::vector<std::vector<double>> w(1);
		w[0].resize(task_size);
		for (int i=0; i<task_size; ++i) w[0][i] = 1;
		return w;
	} else if (node == n-1 && side == "R") {
		std::vector<std::vector<double>> w(1);
		w[0].resize(task_size);
		for (int i=0; i<task_size; ++i) w[0][i] = 1;
		return w;
	} else if (node == 1 && side == "R") {
		std::vector<double> beta0(task_size);
		std::vector<double> beta1(task_size);
		for (int i=0; i<task_size; ++i) {
			beta0[i] = (cells[2].u[i]-cells[1].u[i])*(cells[2].u[i]-cells[1].u[i]);
			beta1[i] = (cells[2].u[i]-cells[3].u[i])*(cells[2].u[i]-cells[3].u[i]);
		}
		std::vector<double> alpha0(task_size);
		std::vector<double> alpha1(task_size);
		for (int i=0; i<task_size; ++i) {
			alpha0[i] = 2./3./((1e-6 + beta0[i])*(1e-6 + beta0[i]));
			alpha1[i] = 1./3./((1e-6 + beta1[i])*(1e-6 + beta1[i]));
		}
		std::vector<double> w0(task_size);
		std::vector<double> w1(task_size);
		for (int i=0; i<task_size; ++i) {
			w0[i] = alpha0[i] / (alpha0[i]+alpha1[i]);
			w1[i] = alpha1[i] / (alpha0[i]+alpha1[i]);
		}
		std::vector<std::vector<double>> w = {w0, w1};
		return w;
	} else if (node == 2 && side == "L") {
		std::vector<double> beta0(task_size);
		std::vector<double> beta1(task_size);
		for (int i=0; i<task_size; ++i) {
			beta0[i] = (cells[3].u[i]-cells[2].u[i])*(cells[3].u[i]-cells[2].u[i]);
			beta1[i] = (cells[2].u[i]-cells[1].u[i])*(cells[2].u[i]-cells[1].u[i]);
		}
		std::vector<double> alpha0(task_size);
		std::vector<double> alpha1(task_size);
		for (int i=0; i<task_size; ++i) {
			alpha0[i] = 2./3./((1e-6 + beta0[i])*(1e-6 + beta0[i]));
			alpha1[i] = 1./3./((1e-6 + beta1[i])*(1e-6 + beta1[i]));
		}
		std::vector<double> w0(task_size);
		std::vector<double> w1(task_size);
		for (int i=0; i<task_size; ++i) {
			w0[i] = alpha0[i] / (alpha0[i]+alpha1[i]);
			w1[i] = alpha1[i] / (alpha0[i]+alpha1[i]);
		}
		std::vector<std::vector<double>> w = {w0, w1};
		return w;
	} else if (node == n-1 && side == "L") {
		std::vector<double> beta0(task_size);
		std::vector<double> beta1(task_size);
		for (int i=0; i<task_size; ++i) {
			beta0[i] = (cells[n].u[i]-cells[n-1].u[i])*(cells[n].u[i]-cells[n-1].u[i]);
			beta1[i] = (cells[n-1].u[i]-cells[n-2].u[i])*(cells[n-1].u[i]-cells[n-2].u[i]);
		}
		std::vector<double> alpha0(task_size);
		std::vector<double> alpha1(task_size);
		for (int i=0; i<task_size; ++i) {
			alpha0[i] = 2./3./((1e-6 + beta0[i])*(1e-6 + beta0[i]));
			alpha1[i] = 1./3./((1e-6 + beta1[i])*(1e-6 + beta1[i]));
		}
		std::vector<double> w0(task_size);
		std::vector<double> w1(task_size);
		for (int i=0; i<task_size; ++i) {
			w0[i] = alpha0[i] / (alpha0[i]+alpha1[i]);
			w1[i] = alpha1[i] / (alpha0[i]+alpha1[i]);
		}
		std::vector<std::vector<double>> w = {w0, w1};
		return w;
	} else if (node == n-2 && side == "R") {
		std::vector<double> beta0(task_size);
		std::vector<double> beta1(task_size);
		for (int i=0; i<task_size; ++i) {
			beta0[i] = (cells[n-2].u[i]-cells[n-1].u[i])*(cells[n-2].u[i]-cells[n-1].u[i]);
			beta1[i] = (cells[n-1].u[i]-cells[n].u[i])*(cells[n-1].u[i]-cells[n].u[i]);
		}
		std::vector<double> alpha0(task_size);
		std::vector<double> alpha1(task_size);
		for (int i=0; i<task_size; ++i) {
			alpha0[i] = 2./3./((1e-6 + beta0[i])*(1e-6 + beta0[i]));
			alpha1[i] = 1./3./((1e-6 + beta1[i])*(1e-6 + beta1[i]));
		}
		std::vector<double> w0(task_size);
		std::vector<double> w1(task_size);
		for (int i=0; i<task_size; ++i) {
			w0[i] = alpha0[i] / (alpha0[i]+alpha1[i]);
			w1[i] = alpha1[i] / (alpha0[i]+alpha1[i]);
		}
		std::vector<std::vector<double>> w = {w0, w1};
		return w;
	} else {
		if (side == "L") {
			std::vector<double> beta0(task_size);
			std::vector<double> beta1(task_size);
			std::vector<double> beta2(task_size);
			for (int i=0; i<task_size; ++i) {
				beta0[i] = 13./12.*(cells[node].u[i]-2*cells[node+1].u[i]+cells[node+2].u[i])
								  *(cells[node].u[i]-2*cells[node+1].u[i]+cells[node+2].u[i])
					+ 1./4.*(3*cells[node].u[i]-4*cells[node+1].u[i]+cells[node+2].u[i])
						   *(3*cells[node].u[i]-4*cells[node+1].u[i]+cells[node+2].u[i]);

				beta1[i] = 13./12.*(cells[node-1].u[i]-2*cells[node].u[i]+cells[node+1].u[i])
								  *(cells[node-1].u[i]-2*cells[node].u[i]+cells[node+1].u[i])
					+ 1./4.*(cells[node-1].u[i]-cells[node+1].u[i])
						   *(cells[node-1].u[i]-cells[node+1].u[i]);

				beta2[i] = 13./12.*(cells[node-2].u[i]-2*cells[node-1].u[i]+cells[node].u[i])
								  *(cells[node-2].u[i]-2*cells[node-1].u[i]+cells[node].u[i])
					+ 1./4.*(cells[node-2].u[i]-4*cells[node-1].u[i]+3*cells[node].u[i])
						   *(cells[node-2].u[i]-4*cells[node-1].u[i]+3*cells[node].u[i]);
			}
			std::vector<double> alpha0(task_size);
			std::vector<double> alpha1(task_size);
			std::vector<double> alpha2(task_size);
			std::vector<double> w0(task_size);
			std::vector<double> w1(task_size);
			std::vector<double> w2(task_size);
			for (int i=0; i<task_size; ++i) {
				alpha0[i] = 3./10./((1e-6 + beta0[i])*(1e-6 + beta0[i]));
				alpha1[i] = 3./5./((1e-6 + beta1[i])*(1e-6 + beta1[i]));
				alpha2[i] = 1./10./((1e-6 + beta2[i])*(1e-6 + beta2[i]));
				w0[i] = alpha0[i] / (alpha0[i]+alpha1[i]+alpha2[i]);
				w1[i] = alpha1[i] / (alpha0[i]+alpha1[i]+alpha2[i]);
				w2[i] = alpha2[i] / (alpha0[i]+alpha1[i]+alpha2[i]);
			}
			std::vector<std::vector<double>> w = {w0, w1, w2};
			return w;
		} else if (side == "R") {
			std::vector<double> beta0(task_size);
			std::vector<double> beta1(task_size);
			std::vector<double> beta2(task_size);
			for (int i=0; i<task_size; ++i) {
				beta0[i] = 13./12.*(cells[node+1].u[i]-2*cells[node].u[i]+cells[node-1].u[i])
								  *(cells[node+1].u[i]-2*cells[node].u[i]+cells[node-1].u[i])
					+ 1./4.*(3*cells[node+1].u[i]-4*cells[node].u[i]+cells[node-1].u[i])
						   *(3*cells[node+1].u[i]-4*cells[node].u[i]+cells[node-1].u[i]);

				beta1[i] = 13./12.*(cells[node+2].u[i]-2*cells[node+1].u[i]+cells[node].u[i])
								  *(cells[node+2].u[i]-2*cells[node+1].u[i]+cells[node].u[i])
					+ 1./4.*(cells[node+2].u[i]-cells[node].u[i])
						   *(cells[node+2].u[i]-cells[node].u[i]);

				beta2[i] = 13./12.*(cells[node+3].u[i]-2*cells[node+2].u[i]+cells[node+1].u[i])
								  *(cells[node+3].u[i]-2*cells[node+2].u[i]+cells[node+1].u[i])
					+ 1./4.*(cells[node+3].u[i]-4*cells[node+2].u[i]+3*cells[node+1].u[i])
						   *(cells[node+3].u[i]-4*cells[node+2].u[i]+3*cells[node+1].u[i]);
			}
			std::vector<double> alpha0(task_size);
			std::vector<double> alpha1(task_size);
			std::vector<double> alpha2(task_size);
			std::vector<double> w0(task_size);
			std::vector<double> w1(task_size);
			std::vector<double> w2(task_size);
			for (int i=0; i<task_size; ++i) {
				alpha0[i] = 3./10./((1e-6 + beta0[i])*(1e-6 + beta0[i]));
				alpha1[i] = 3./5./((1e-6 + beta1[i])*(1e-6 + beta1[i]));
				alpha2[i] = 1./10./((1e-6 + beta2[i])*(1e-6 + beta2[i]));
				w0[i] = alpha0[i] / (alpha0[i]+alpha1[i]+alpha2[i]);
				w1[i] = alpha1[i] / (alpha0[i]+alpha1[i]+alpha2[i]);
				w2[i] = alpha2[i] / (alpha0[i]+alpha1[i]+alpha2[i]);
			}
			std::vector<std::vector<double>> w = {w0, w1, w2};
			return w;
		}
	}
}

std::vector<std::vector<double>> solver::calc_face_value (const int node, const std::string side,
                       const std::vector<cell>& cells, const int task_size) {

	int n = cells.size()-2;
	if (node == 0 && side == "R") {
		std::vector<double> u0(task_size);
		for (int i=0; i<task_size; ++i) {
			u0[i] = 1./3*cells[3].u[i] - 7./6*cells[2].u[i] + 11./6*cells[1].u[i];
		}
		std::vector<std::vector<double>> u = {u0};
		return u;
	} else if (node == n && side == "L") {
		std::vector<double> u0(task_size);
		for (int i=0; i<task_size; ++i) {
			u0[i] = 1./3*cells[n-2].u[i] - 7./6*cells[n-1].u[i] + 11./6*cells[n].u[i];
		}
		std::vector<std::vector<double>> u = {u0};
		return u;
	} else if (node == 1 && side == "L") {
		std::vector<double> u0(task_size);
		for (int i=0; i<task_size; ++i) {
			u0[i] = 1./3*cells[1].u[i] + 5./6*cells[2].u[i] - 1./6*cells[3].u[i];
		}
		std::vector<std::vector<double>> u = {u0};
		return u;
	} else if (node == n-1 && side == "R") {
		std::vector<double> u0(task_size);
		for (int i=0; i<task_size; ++i) {
			u0[i] = 1./3*cells[n].u[i] + 5./6*cells[n-1].u[i] - 1./6*cells[n-2].u[i];
		}
		std::vector<std::vector<double>> u = {u0};
		return u;
	} else if (node == 1 && side == "R") {
		std::vector<double> u0(task_size);
		std::vector<double> u1(task_size);
		for (int i=0; i<task_size; ++i) {
			u0[i] = -1./6*cells[3].u[i] + 5./6*cells[2].u[i] + 1./3*cells[1].u[i];
			u1[i] = 1./3*cells[4].u[i] - 7./6*cells[3].u[i] + 11./6*cells[2].u[i];
		}
		std::vector<std::vector<double>> u = {u0, u1};
		return u;
	} else if (node == 2 && side == "L") {
		std::vector<double> u0(task_size);
		std::vector<double> u1(task_size);
		for (int i=0; i<task_size; ++i) {
			u0[i] = 1./3*cells[2].u[i] + 5./6*cells[3].u[i] - 1./6*cells[4].u[i];
			u1[i] = -1./6*cells[1].u[i] + 5./6*cells[2].u[i] + 1./3*cells[3].u[i];
		}
		std::vector<std::vector<double>> u = {u0, u1};
		return u;
	} else if (node == n-1 && side == "L") {
		std::vector<double> u0(task_size);
		std::vector<double> u1(task_size);
		for (int i=0; i<task_size; ++i) {
			u0[i] = -1./6*cells[n-2].u[i] + 5./6*cells[n-1].u[i] + 1./3*cells[n].u[i];
			u1[i] = 1./3*cells[n-3].u[i] - 7./6*cells[n-2].u[i] + 11./6*cells[n-1].u[i];
		}
		std::vector<std::vector<double>> u = {u0, u1};
		return u;
	} else if (node == n-2 && side == "R") {
		std::vector<double> u0(task_size);
		std::vector<double> u1(task_size);
		for (int i=0; i<task_size; ++i) {
			u0[i] = 1./3*cells[n-1].u[i] + 5./6*cells[n-2].u[i] - 1./6*cells[n-3].u[i];
			u1[i] = -1./6*cells[n].u[i] + 5./6*cells[n-1].u[i] + 1./3*cells[n-2].u[i];
		}
		std::vector<std::vector<double>> u = {u0, u1};
		return u;
	} else {
		if (side == "L") {
			std::vector<double> u0(task_size);
			std::vector<double> u1(task_size);
			std::vector<double> u2(task_size);
			for (int i=0; i<task_size; ++i) {
				u0[i] = 1./3*cells[node].u[i] + 5./6*cells[node+1].u[i] - 1./6*cells[node+2].u[i];
				u1[i] = -1./6*cells[node-1].u[i] + 5./6*cells[node].u[i] + 1./3*cells[node+1].u[i];
				u2[i] = 1./3*cells[node-2].u[i] - 7./6*cells[node-1].u[i] + 11./6*cells[node].u[i];
			}
			std::vector<std::vector<double>> u = {u0, u1, u2};
			return u;
		} else if (side == "R") {
			std::vector<double> u0(task_size);
			std::vector<double> u1(task_size);
			std::vector<double> u2(task_size);
			for (int i=0; i<task_size; ++i) {
				u0[i] = 1./3*cells[node+1].u[i] + 5./6*cells[node].u[i] - 1./6*cells[node-1].u[i];
				u1[i] = -1./6*cells[node+2].u[i] + 5./6*cells[node+1].u[i] + 1./3*cells[node].u[i];
				u2[i] = 1./3*cells[node+3].u[i] - 7./6*cells[node+2].u[i] + 11./6*cells[node+1].u[i];
			}
			std::vector<std::vector<double>> u = {u0, u1, u2};
			return u;
		}
	}
}

void solver::solve_weno(grid& mesh, std::string out) {
	std::vector<cell> cells = mesh.get_cells();
	double cfl = mesh.get_input_quantities().cfl;
	double gamma = cells[0].gamma;
	int ncells = mesh.get_input_quantities().ncells;
	double tf = mesh.get_input_quantities().tf;
	double dx = mesh.get_dx();

    int task_size = mesh.get_input_quantities().get_task_size();

	std::vector<double> flux(task_size);
    for (auto& elem : cells) elem.res.resize(task_size);

	for (int itime=1; itime<1e6; ++itime) {
		if (t == tf) break;
		double dt = timestep(cfl, dx, gamma, ncells, cells);
		if (t+dt > tf) dt = tf - t;
		t += dt;
		nsteps += 1;

		//обнуляем невязку на каждом шаге по времени в каждой ячейке
		for (int j=1; j<ncells+1; ++j) std::fill(cells[j].res.begin(), cells[j].res.end(), 0.);
				
		for (int i=1; i<ncells; ++i) {
			std::vector<double> uL(task_size);
			std::vector<double> uR(task_size);
			
			std::vector<std::vector<double>> weightR = calc_weight(i, "R", cells, task_size);
			std::vector<std::vector<double>> valueR = calc_face_value(i, "R", cells, task_size);
			std::vector<std::vector<double>> weightL = calc_weight(i, "L", cells, task_size);
			std::vector<std::vector<double>> valueL = calc_face_value(i, "L", cells, task_size);

			for (int j=0; j<task_size; ++j) {
				for (int k=0; k<weightR.size(); ++k) {
					uR[j]+=weightR[k][j]*valueR[k][j];
				}
				for (int k=0; k<weightL.size(); ++k) {
					uL[j]+=weightL[k][j]*valueL[k][j];
				}
			}
			flux = hll_flux(uL, uR, gamma); //find flux
			

			// residual upadate now
			for (int j=0; j<3; ++j) {
				cells[i].res[j] += flux[j];
				cells[i+1].res[j] -= flux[j];
			}
		}

		// Left most face: left face of cell i=1
		std::vector<double> uR = calc_face_value(0, "R", cells, task_size)[0];
		flux = hll_flux(uR, uR, gamma);
		for (int i=0; i<3; ++i) cells[1].res[i] -= flux[i];
		
		// Right most face: right face of cell i=ncells.
		std::vector<double> uL = calc_face_value(ncells, "L", cells, task_size)[0];
		flux = hll_flux(uL, uL, gamma);
		for (int i=0; i<3; ++i) cells[ncells].res[i] += flux[i];
		
		// Solution update
		for (int i=1; i<ncells+1; ++i) {
			for (int j=0; j<3; ++j) cells[i].u[j] -= dt/dx * cells[i].res[j];
			cells[i].w = u2w(cells[i].u, gamma);
			
		}
		// Copy the solutions to the ghost cells.
		for (int i=0; i<3; ++i) {
			cells[0].w[i] = cells[1].w[i];
			cells[ncells+1].w[i] = cells[ncells].w[i];
		}
	}

	std::ofstream f(out);
	for (int i=1; i<ncells+1; ++i) f << cells[i].xc << ", " << cells[i].w[0] << ", " 
									 << cells[i].w[1] << ", " << cells[i].w[2] << '\n';
	f.close();
}