#ifndef __ERF2D1_H__
#define __ERF2D1_H__

#include "includes.h"
#include "unsteady_heat.h"

class Erf2D1Problem {
	public:
		static constexpr double Lx = 1.0;
		static constexpr double Ly = 1.0;
		uint Nx;
		uint Ny;
		double dt;
		uint t_steps;
		double t_shift;

		unordered_map<uint, FluxFctPt> flux_boundaries = {
			// {0, get_flux},
			// {2, get_flux}
		};

		// pin center to 0.0
		unordered_map<uint, double> pinned_boundaries = {
			// {2, 0.0}
		};

		// use analytical solution on substrate and liquid boundaries
		Vector<uint> analytical_boundaries = {
			0,
			1,
			2,
			3
		};

		static constexpr double kappa1 = 1.0;
		static constexpr double kappa2 = 0.5;
		static constexpr double Ts = -1.0;
		static constexpr double Tfr = 1.0;

		Erf2D1Problem(uint Nx, uint Nx2, uint Ny, uint t_steps, double dt, double t_shift)
			: Nx(Nx), Ny(Ny), dt(dt), t_steps(t_steps), t_shift(t_shift) {}
		
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

		static void get_flux(const double &t, const Vector<double> &x, double &flux) {
			flux = 0.0;
		}
};

#endif