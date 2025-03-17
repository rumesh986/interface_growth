#include <cmath>

#include "generic.h"
#include "meshes/one_d_mesh.h"
#include "unsteady_heat.h"

using namespace std;
using namespace oomph;
using namespace MathematicalConstants;

typedef Vector<double> vecd;

namespace ProblemParams {
	double D1 = 1.0;
	double D2 = 2.0;

	const double Ts = -1.0;
	const double Tfr = 1.0;

	const double lambda_1 = 1.0;
	const double lambda_2 = 2.0;

	const double e_rat = sqrt(lambda_1/lambda_2);
	// const double e_rat = lambda_1/lambda_2;

	const double eps = 1e-6;
	const double ic_alpha = 2.5;

	const double t_shift = 0.015;

	double T0(double t) {
		return Ts * erf(eps/(2*sqrt(D1*t)));
	}

	void get_exact_u(const double &t, const vecd &x, vecd &u) {
		double t0 = T0(t);
		double D = (x[0] < 0) ? D1 : D2;
		double e = (x[0] < 0) ? 1 : e_rat;

		u[0] = t0 + (t0 - Ts) * erf(x[0]/(2*sqrt(D*t))) * e;
	}

	void get_source(const double &t, const vecd &x, double &source) {
		source = 0.0;
	}

	void get_initial_condition(const double &t, const vecd &x, vecd &u) {
		// double T = (x[0] < 0) ? -Ts : Tfr;

		// u[0] = T * erf(ic_alpha * x[0]);
		get_exact_u(t, x, u);
	}
}

template<class EL> class TwoLayerMesh : public virtual OneDMesh<EL> {
	public:
		TwoLayerMesh(const uint nx1, const uint nx2, const double x0, const double x1, const double x2, TimeStepper *ts_pt = &Mesh::Default_TimeStepper)
			: OneDMesh<EL>(nx1+nx2, x0, x2, ts_pt) {
				this->set_nboundary(3);

				FiniteElement *el_pt = this->finite_element_pt(nx1);
				// const uint nnode = el_pt->nnode();
				Node *node_pt = el_pt->node_pt(0);
				this->convert_to_boundary_node(node_pt);
				this->add_boundary_node(2, node_pt);

				this->setup_boundary_element_info();
			}
};

template<class EL> class TwoMeshUnsteadyHeatProblem : public Problem {
	public:
		TwoMeshUnsteadyHeatProblem(
			const uint nx1, const uint nx2, 
			const double x_0, const double x_1, const double x_2, 
			UnsteadyHeatEquations<1>::UnsteadyHeatSourceFctPt source_inp)
				: Nx1(nx1), Nx2(nx2), x0(x_0), x1(x_1), x2(x_2), source_pt(source_inp) {
					add_time_stepper_pt(new BDF<2>);

					bulk_mesh_pt = new TwoLayerMesh<EL>(Nx1, Nx2, x0, x1, x2, time_stepper_pt());
					surf_mesh_pt = new Mesh;

					// create_interface_elements();

					add_sub_mesh(bulk_mesh_pt);
					add_sub_mesh(surf_mesh_pt);

					build_global_mesh();

					nbound = bulk_mesh_pt->nboundary();

					// set dirichlet conditions
					for (uint b = 0; b < nbound; b++) {
						uint nnode = bulk_mesh_pt->nboundary_node(b);
						for (uint n = 0; n < nnode; n++) {
							bulk_mesh_pt->boundary_node_pt(b, n)->pin(0);
							bulk_mesh_pt->boundary_node_pt(b, n)->set_value(0, (double) ((int)b)-1);

							cout << "Value at b = " << b << " : " << bulk_mesh_pt->boundary_node_pt(b, n)[0] << endl;
						}
					}

					// set coefficient values
					for (uint e = 0; e < Nx1; e++) {
						EL *el_pt = dynamic_cast<EL*>(bulk_mesh_pt->element_pt(e));
						el_pt->beta_pt() = &ProblemParams::D1;
						el_pt->source_fct_pt() = source_pt;
					}

					for (uint e = Nx1; e < Nx1+Nx2; e++) {
						EL *el_pt = dynamic_cast<EL*>(bulk_mesh_pt->element_pt(e));
						el_pt->beta_pt() = &ProblemParams::D2;
						el_pt->source_fct_pt() = source_pt;
					}

					cout << "Number of equations : " << assign_eqn_numbers() << endl;
				}

		void actions_before_newton_solve() {}
		void actions_after_newton_solve() {}
		
		void actions_before_implicit_timestep() {
			double t = time_pt()->time() + ProblemParams::t_shift;

			vecd x(1);
			vecd u(1);

			for (uint b = 0; b < nbound; b++) {
				int nnode = bulk_mesh_pt->nboundary_node(b);
				for (int n = 0; n < nnode; n++) {
					Node *node_pt = bulk_mesh_pt->boundary_node_pt(b, n);
					x[0] = node_pt->x(0);

					ProblemParams::get_exact_u(t, x, u);
					node_pt->set_value(0, 0, u[0]);
				}
			}
		}

		void actions_after_implicit_timestep() {}

		void set_initial_condition() {
			double t_bak = time_pt()->time();

			vecd x(1);
			vecd u(1);

			int nprev_steps = time_stepper_pt()->nprev_values();

			for (int t = nprev_steps; t >= 0; t--) {
				double t_prev = time_pt()->time((uint)t) + ProblemParams::t_shift;
				cout << "Setting Initial Condition at time = " << t_prev << endl;

				// bulk_mesh_pt->node_pt(0)->set_value(t, 0, ProblemParams::Ts);
				// bulk_mesh_pt->node_pt(1)->set_value(t, 0, ProblemParams::Ts);

				// for (int n = 2; n < nnode; n++) {
				// 	x[0] = bulk_mesh_pt->node_pt(n)->x(0);
				// 	bulk_mesh_pt->node_pt(n)->set_value(t, 0, ProblemParams::Tfr);
				// }

				// for (uint e = 0; e < Nx1; e++) {
				// 	EL *el_pt = dynamic_cast<EL*>(bulk_mesh_pt->element_pt(e));
				// 	uint nnode = el_pt->nnode();
				// 	for (uint n = 0; n < nnode; n++)
				// 		el_pt->node_pt(n)->set_value(0, ProblemParams::Ts);
				// 		// el_pt->node_pt(n)->set_value(t, 0, ProblemParams::Ts);
				// }

				// for (uint e = Nx1; e < Nx1+Nx2; e++) {
				// 	EL *el_pt = dynamic_cast<EL*>(bulk_mesh_pt->element_pt(e));
				// 	uint nnode = el_pt->nnode();
				// 	for (uint n = 0; n < nnode; n++)
				// 		el_pt->node_pt(n)->set_value(0, ProblemParams::Tfr);//*ProblemParams::e_rat);
				// 		// el_pt->node_pt(n)->set_value(t, 0, ProblemParams::Tfr);
				// }
				for (uint e = 0; e < Nx1+Nx2; e++) {
					EL *el_pt = dynamic_cast<EL*>(bulk_mesh_pt->element_pt(e));
					uint nnode = el_pt->nnode();
					for (uint n = 0; n < nnode; n++) {
						x[0] = el_pt->node_pt(n)->x(0);
						ProblemParams::get_initial_condition(t_prev, x, u);
						el_pt->node_pt(n)->set_value(t, 0, u[0]);
					}
				}

			}

			time_pt()->time() = t_bak;

			// assign_initial_values_impulsive();
		}

		void doc_solution(DocInfo &info) {
			const uint npts = 5;

			char filename[100];
			ofstream outfile;
			// char timestring[100];

			sprintf(filename, "%s/soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			mesh_pt()->output(outfile, npts);
			outfile.close();

			double time = time_pt()->time()+ProblemParams::t_shift;

			sprintf(filename, "%s/exact_soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			mesh_pt()->output_fct(outfile, npts, time, ProblemParams::get_exact_u);
			outfile.close();

			double error, norm;
			sprintf(filename, "%s/error%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			mesh_pt()->compute_error(outfile, ProblemParams::get_exact_u, time, error, norm);
			outfile.close();

			sprintf(filename, "%s/times.dat", info.directory().c_str());
			outfile.open(filename, ios::app);
			// sprintf(timestring, "%10.5f\n", time);
			// fwrite(timestring, sizeof(char), 100, outfile);
			outfile << time << endl;
			outfile.close();

			info.number()++;
		}

	private:
		uint Nx1;
		uint Nx2;
		double x0;
		double x1;
		double x2;
		UnsteadyHeatEquations<1> :: UnsteadyHeatSourceFctPt source_pt;

		uint nbound;
		TwoLayerMesh<EL> *bulk_mesh_pt;
		Mesh *surf_mesh_pt;
};

int main() {
	TwoMeshUnsteadyHeatProblem<QUnsteadyHeatElement<1,2>> problem(1000, 1000, -1.0, 0.0, 1.0, ProblemParams::get_source);

	DocInfo info;
	info.set_directory("RESLT");
	info.number() = 0;

	const double t_max = 1.0;
	const double dt = 0.005;

	problem.initialise_dt(dt);
	problem.set_initial_condition();
	problem.doc_solution(info);

	uint nsteps = (uint) t_max / dt;
	for (uint t = 0; t < nsteps; t++) {
		cout << "Timestep " << t << endl;

		problem.unsteady_newton_solve(dt);
		problem.doc_solution(info);
	}
}