#ifndef __ERF2D5_H__
#define __ERF2D5_H__

#include "includes.h"
#include "unsteady_heat.h"

class Erf2D5Problem {
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
		static constexpr double vel = 0.5;

		unordered_map<uint, FluxFctPt> flux_boundaries = {
			// {0, get_flux},
			// {2, get_flux}
			// {4, get_flux}
		};

		// pin center to 0.0
		unordered_map<uint, double> pinned_boundaries = {
			// {4, 0.0}
			{1, Tm},
			{3, Ts},
			{4, Tm}
		};

		// use analytical solution on substrate and liquid boundaries
		Vector<uint> analytical_boundaries = {
			// 0,
			// 1,
			// 2,
			// 3,
			// 4
		};

		static constexpr double kappa1 = 1.0;
		static constexpr double kappa2 = 0.5;
		static constexpr double Ts = -1.0;
		static constexpr double Tfr = 1.0;
		static constexpr double Tm = 0.0;
		static constexpr double alpha = 1.0;

		Erf2D5Problem(uint Nx1, uint Nx2, uint Ny, uint t_steps, double dt, double t_shift)
			: Nx1(Nx1), Nx2(Nx2), Nx(Nx1+Nx2), Ny(Ny), dt(dt), t_steps(t_steps), t_shift(t_shift) {
				// flux_pt = malloc(sizeof(double));
			}
		
		// x: 2D vector, u: 1D vector
		static void get_exact_u(const double &t, const Vector<double> &x, Vector<double> &u) {
			if (x[0] <= vel * t) {
				// u[0] = sqrt(Pi) * alpha * vel * exp(0.25 * vel * vel * t / kappa1) * erf(x[0] / (2 * sqrt(kappa1 * t))) + Ts;
				u[0] = Ts + (Tm - Ts) * erf(0.5 * x[0]/sqrt(kappa1 * t)) / erf(0.5 * vel * t / sqrt(kappa1 * t));
			} else {
				u[0] = Tm;
			}
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