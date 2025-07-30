#ifndef __ERF2D4_H__
#define __ERF2D4_H__

#include "includes.h"
#include "unsteady_heat.h"

class Erf2D4Problem {
	public:
		static constexpr double Lx = 1.0;
		static constexpr double Ly = 1.0;
		uint Nx1;
		uint Nx2;
		uint Nx;
		uint Ny;
		double dt;
		uint t_steps;
		double t_shift;
		static constexpr double vel = -0.5;

		// double *flux_pt = (double *) malloc(sizeof(double));
		static constexpr double *flux_pt = nullptr;

		unordered_map<uint, FluxFctPt> flux_boundaries = {
			// {0, get_flux},
			// {2, get_flux}
			{4, get_flux}
		};

		// pin center to 0.0
		unordered_map<uint, double> pinned_boundaries = {
			// {4, 0.0}
		};

		// use analytical solution on substrate and liquid boundaries
		Vector<uint> analytical_boundaries = {
			// 0,
			1,
			// 2,
			3,
			// 4
		};

		static constexpr double kappa1 = 1.0;
		static constexpr double kappa2 = 0.5;
		static constexpr double Ts = -1.0;
		static constexpr double Tfr = 1.0;

		Erf2D4Problem(uint Nx1, uint Nx2, uint Ny, uint t_steps, double dt, double t_shift)
			: Nx1(Nx1), Nx2(Nx2), Nx(Nx1+Nx2), Ny(Ny), dt(dt), t_steps(t_steps), t_shift(t_shift) {
				// flux_pt = malloc(sizeof(double));
			}
		
		// x: 2D vector, u: 1D vector
		static void get_exact_u(const double &t, const Vector<double> &x, Vector<double> &u) {
			double t0 = 0.0;
			// double xi = x[0] / sqrt(t);
			double kappa = (x[0] < 0) ? kappa1 : kappa2;
			double e = (x[0] < 0) ? 1 : sqrt(kappa1 / kappa2);
			// double e = 1;

			u[0] = t0 - (t0 - Ts) * e * erf(x[0] / (2 * sqrt(kappa * t)));
		}

		static void get_source(const double &t, const Vector<double> &x, double &source) {
			source = 0.0;
		}

		static void get_initial_condition(const double &t, const Vector<double> &x, Vector<double> &u) {
			get_exact_u(t, x, u);
		}

		static constexpr void get_flux(const double &t, const Vector<double> &x, double &flux) {
			// flux = -0.1;

			flux = 0.2 * vel;

			// flux = 0.0;

		// 	if (flux_pt)
		// 		flux = *flux_pt;
		}
};

#endif