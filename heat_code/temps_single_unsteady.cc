#include <cmath>

#include "generic.h"
#include "meshes/one_d_mesh.h"
#include "unsteady_heat.h"

using namespace std;
using namespace oomph;
using namespace MathematicalConstants;

#define bnode_pt(b, n) mesh_pt()->boundary_node_pt(b, n)

typedef Vector<double> vecd;

namespace ProblemParams {
	const double alpha = 1.0;
	const double k = 0.3;

	void get_exact_u(const double &t, const vecd &x, vecd &u) {
		// u[0] = erf((2*x[0])/(alpha*sqrt(t)));
		// u[0] = exp(-k*k*t) * sin(k*x[0]) / sin(k);
		double val = 2.0/(alpha*sqrt(t));
		u[0] = erf(val * x[0])/erf(val);
	}

	void get_source(const double &t, const vecd &x, double &source) {
		source = 0.0;
	}
}

template<class EL> class OneDUnsteadyHeatProblem : public Problem {
	public:
		OneDUnsteadyHeatProblem(const uint Nx, const double Lx, UnsteadyHeatEquations<1>::UnsteadyHeatSourceFctPt source_inp) : source_pt(source_inp), nelem(Nx) {
			add_time_stepper_pt(new BDF<2>);

			mesh_pt() = new OneDMesh<EL>(nelem, Lx, time_stepper_pt());

			nbound = mesh_pt()->nboundary();
			for (uint b = 0; b < nbound; b++) {
				int nnode = mesh_pt()->nboundary_node(b);
				for (int n = 0; n < nnode; n++)
					mesh_pt()->boundary_node_pt(b, n)->pin(0);
			}

			for (uint i = 0; i < Nx; i++)
				dynamic_cast<EL *>(mesh_pt()->element_pt(i))->source_fct_pt() = source_pt;

			cout << "Number of equations: " << assign_eqn_numbers() << endl;
		}

		~OneDUnsteadyHeatProblem() {
			delete mesh_pt();
		}

		void actions_before_newton_solve() {}
		void actions_after_newton_solve() {}
		void actions_before_implicit_timestep() {
			double t = time_pt()->time();

			for (uint b = 0; b < nbound; b++) {
				int nnode = mesh_pt()->nboundary_node(b);
				for (int n = 0; n < nnode; n++) {
					vecd x(1), u(1);
					x[0] = mesh_pt()->boundary_node_pt(b, n)->x(0);
					ProblemParams::get_exact_u(t, x, u);
					bnode_pt(b, n)->set_value(0, u[0]);
				}
			}
		}
		
		void actions_after_implicit_timestep() {}

		void set_initial_conditions() {
			double t_bak = time_pt()->time();

			vecd x(1);

			int nnode = mesh_pt()->nnode();
			int nprev_steps = time_stepper_pt()->nprev_values();

			for (int t = nprev_steps; t >= 0; t--) {
				// double time = prev_time[t];
				cout << "Setting Initial Condition at time = " << time_pt()->time((uint)t) << endl;

				for (int n = 2; n < nnode; n++) {
					x[0] = mesh_pt()->node_pt(n)->x(0);
					mesh_pt()->node_pt(n)->set_value(t, 0, 1.0);
				}

				mesh_pt()->node_pt(0)->set_value(t, 0, 0.0);
				mesh_pt()->node_pt(1)->set_value(t, 0, 0.0);
			}
			cout << nnode << endl;

			time_pt()->time() = t_bak;

		}

		void doc_solution(DocInfo &info) {
			const uint npts = 5;

			char filename[100];
			ofstream outfile;

			sprintf(filename, "%s/soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			mesh_pt()->output(outfile, npts);
			outfile.close();

			sprintf(filename, "%s/exact_soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			mesh_pt()->output_fct(outfile, npts, time_pt()->time(), ProblemParams::get_exact_u);
			outfile.close();

			double error, norm;
			sprintf(filename, "%s/error%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			mesh_pt()->compute_error(outfile, ProblemParams::get_exact_u, time_pt()->time(), error, norm);
			outfile.close();

			info.number()++;
		}
	
	private:
		UnsteadyHeatEquations<1>::UnsteadyHeatSourceFctPt source_pt;
		uint nbound;
		uint nelem;
};

int main() {
	OneDUnsteadyHeatProblem<QUnsteadyHeatElement<1,2>> problem(1000, 1.0, ProblemParams::get_source);

	DocInfo info;
	info.set_directory("RESLT");
	info.number() = 0;

	const double t_max = 5.0;
	const double dt = 0.05;

	problem.initialise_dt(dt);
	problem.set_initial_conditions();
	problem.doc_solution(info);


	uint nsteps = (uint) t_max / dt;
	for (uint t = 0; t < nsteps; t++) {
		cout << "Timestep " << t << endl;

		problem.unsteady_newton_solve(dt);
		problem.doc_solution(info);
	}
}