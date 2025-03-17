#include <cmath>

#include "generic.h"
#include "solid.h"
#include "meshes/one_d_mesh.h"
#include "unsteady_heat.h"

using namespace std;
using namespace oomph;
using namespace MathematicalConstants;

namespace ExactSoln {
	double Ts = -10.0; // substrate temp
	double T_fr = 10.0; // temp at +inf

	double lambda_s = 1.0;
	double lambda_i = 2.0;

	double Ds = 10.0;
	double Di = 20.0;

	double Ds_r = sqrt(Ds);
	double Di_r = sqrt(Di);

	double e_rat = sqrt(lambda_s/lambda_i); // effusivity ratio

	double eps = 1e-6;

	double T0(double time) {
		return Ts * erf(eps/(2*sqrt(Ds*time)));
	}

	void get_exact_u(const double &time, const Vector<double> &x, Vector<double> &u) {
		double t0 = T0(time);
		double D = (x[0] < 0) ? Ds : Di;
		double e = (x[0] < 0) ? 1 : e_rat;

		u[0] = t0 + (t0 - Ts) * erf(x[0]/(2*sqrt(D*time))) * e;
	}

	void get_source(const double &time, const Vector<double> &x, double &source) {
		source = 0.0;
	}
}

template<class EL> class OneDUnsteadyHeatProblem : public Problem {
	public:
		OneDUnsteadyHeatProblem(UnsteadyHeatEquations<1> :: UnsteadyHeatSourceFctPt source_ptr_inp) : source_ptr(source_ptr_inp) {
			add_time_stepper_pt(new BDF<2>);

			mesh_pt() = new OneDMesh<EL>(Nx, Lx, time_stepper_pt());

			int n_bound = mesh_pt()->nboundary();
			for (int b = 0; b < n_bound; b++) {
				int n_node = mesh_pt()->nboundary_node(b);
				for (int n = 0; n < n_node; n++)
					mesh_pt()->boundary_node_pt(b, n)->pin(0);
			}

			unsigned int n_elem = mesh_pt()->nelement();
			for (unsigned int i = 0; i < n_elem; i++) {
				EL *el_pt = dynamic_cast<EL*>(mesh_pt()->element_pt(i));
				el_pt->source_fct_pt() = source_ptr;
			}

			uint control_el = (unsigned int) n_elem / 2;
			control_node_ptr = mesh_pt()->finite_element_pt(control_el)->node_pt(0);

			for (uint i = 0; i < control_el; i++) {
				EL *el_pt = dynamic_cast<EL*>(mesh_pt()->element_pt(i));
				el_pt->alpha_pt() = &ExactSoln::Ds_r;
				el_pt->beta_pt() = &ExactSoln::Ds;
			}

			for (uint i = control_el; i < n_elem; i++) {
				EL *el_pt = dynamic_cast<EL*>(mesh_pt()->element_pt(i));
				el_pt->alpha_pt() = &ExactSoln::Di_r;
				el_pt->beta_pt() = &ExactSoln::Di;
			}

			cout << "Number of equations: " << assign_eqn_numbers() << endl;
		}

		~OneDUnsteadyHeatProblem() {}

		void actions_after_newton_solve() {}
		void actions_before_newton_solve() {}
		void actions_after_implicit_timestep() {}
		void actions_before_implicit_timestep() {
			double time = time_pt()->time();

			int n_bound = mesh_pt()->nboundary();
			for (int b = 0; b < n_bound; b++) {
				int n_node = mesh_pt()->nboundary_node(b);
				for (int n = 0; n < n_node; n++) {
					Node *node_ptr = mesh_pt()->boundary_node_pt(b, n);

					Vector<double> x(1), u(1);
					x[0] = node_ptr->x(0);

					ExactSoln::get_exact_u(time, x, u);
					node_ptr->set_value(0, u[0]);
				}
			}
		}

		void set_initial_conditions() {
			double t_bak = time_pt()->time();

			Vector<double> x(1);

			int n_node = mesh_pt()->nnode();

			int nprev_steps = time_stepper_pt()->nprev_values();
			Vector<double> prev_time(nprev_steps+1);

			for (int t = nprev_steps; t >= 0; t--)
				prev_time[t] = time_pt()->time((unsigned int) t);

			for (int t = nprev_steps; t > -1; t--) {
				double time = prev_time[t];
				cout << "setting IC at time = " << time << endl;

				for (int n = 0; n < n_node; n++) {
					x[0] = mesh_pt()->node_pt(n)->x(0);

					Vector<double> prev_soln(1);
					Vector<double> real_x(1);
					real_x[0] = x[0] - Lx/2.0;

					prev_soln[0] = (real_x[0] < 0.0) ? ExactSoln::Ts : ExactSoln::T_fr;
					mesh_pt()->node_pt(n)->set_value(t, 0, prev_soln[0]);
					mesh_pt()->node_pt(n)->x(t, 0) = real_x[0];
				}
			}

			time_pt()->time() = t_bak;
		}

		void doc_solution(DocInfo &doc_info) {
			const unsigned int npts = 5;
			
			char filename[100];
			ofstream outfile;

			sprintf(filename, "%s/soln%i.dat", doc_info.directory().c_str(), doc_info.number());
			outfile.open(filename);
			mesh_pt()->output(outfile, npts);
			outfile.close();

			sprintf(filename, "%s/exact_soln%i.dat", doc_info.directory().c_str(), doc_info.number());
			outfile.open(filename);
			mesh_pt()->output_fct(outfile, npts, time_pt()->time(), ExactSoln::get_exact_u);
			outfile.close();

			double error, norm;
			sprintf(filename, "%s/error%i.dat", doc_info.directory().c_str(), doc_info.number());
			outfile.open(filename);
			mesh_pt()->compute_error(outfile, ExactSoln::get_exact_u, time_pt()->time(), error, norm);
			outfile.close();
		}

	private:
		const unsigned int Nx = 1000;
		const double Lx = 20.0;

		UnsteadyHeatEquations<1> :: UnsteadyHeatSourceFctPt source_ptr;
		Node *control_node_ptr;
};

int main() {
	OneDUnsteadyHeatProblem<QUnsteadyHeatElement<1,2>> problem(ExactSoln::get_source);

	DocInfo info;
	info.set_directory("RESLT");
	info.number() = 0;

	const double t_max = 10.0;
	const double dt = 0.1;

	problem.initialise_dt(dt);
	problem.set_initial_conditions();
	problem.doc_solution(info);
	info.number()++;

	int n_steps = (unsigned int) t_max / dt;
	for (int t = 0; t < n_steps; t++) {
		cout << "Timestep " << t << endl;

		problem.unsteady_newton_solve(dt);
		problem.doc_solution(info);
		info.number()++;
	}
}