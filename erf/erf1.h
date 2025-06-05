#include "includes.h"
#include "unsteady_heat.h"

// #define T_ORDER 2
// #define X_ORDER 

typedef void (*FluxFctPt)(const double &, const Vector<double>&, double &) ;

class Erf1Problem {
	public:
		static constexpr double Lx = 1.0;
		uint Nx;
		uint t_steps;
		double dt;
		double t_shift;
		
		unordered_map<uint, FluxFctPt> flux_boundaries = {
			// {1, get_flux}
		};

		unordered_map<uint, double> pinned_boundaries = {
			{0, 0.0}
		};

		Vector<uint> analytical_boundaries = {
			1
		};
		
		static constexpr double kappa = 1.0;

		Erf1Problem(uint Nx, uint t_steps, double dt, double t_shift)
			: Nx(Nx), t_steps(t_steps), dt(dt), t_shift(t_shift) {}

		static void get_exact_u(const double &t, const Vector<double> &x, Vector<double> &u) {
			double xi = x[0] / sqrt(t);

			u[0] = erf(xi / (2 * kappa));
		}

		static void get_initial_condition(const double &t, const Vector<double> &x, Vector<double> &u) {
			get_exact_u(t, x, u);
		}

		static void get_source(const double &t, const Vector<double> &x, double &source) {
			source = 0.0;
		}

		static void get_flux(const double &t, const Vector<double> &x, double &flux) {
			flux = 0.0;
		}
};