#include <cmath>

#include "generic.h"
#include "meshes/one_d_mesh.h"
#include "solid.h"
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

	const double eps = 1e-6;

	const double nu = 0.1; // for pseudo solide node thing

	double T0(double t) {
		return Ts * erf(0.5 * eps/sqrt(D1 * t));
	}

	void get_exact_u(const double &t, const vecd &x, vecd &u) {
		double t0 = T0(t);
		double D = (x[0] < 0) ? D1 : D2;
		double e = (x[0] < 0) ? 1 : e_rat;

		u[0] = t0 + (t0 - Ts) * e * erf(0.5 * x[0]/sqrt(D * t));
	}

	void get_source(const double &t, const vecd &u, double &source) {
		source = 0.0;
	}
}

template<class EL> class RefineableTwoLayerMesh : public virtual RefineableOneDMesh<EL> {
	public:
		RefineableTwoLayerMesh(const uint nxs[2], const double xs[3], TimeStepper *ts_pt = &Mesh::Default_TimeStepper)
			: OneDMesh<EL>(nxs[0] + nxs[1], xs[0], xs[2], ts_pt),
			RefineableOneDMesh<EL>(nxs[0] + nxs[1], xs[0], xs[2], ts_pt) {
				this->set_nboundary(3);

				FiniteElement *el_pt = this->finite_element_pt(nxs[0]);
				Node *node_pt = el_pt->node_pt(0);

				this->convert_to_boundary_node(node_pt);
				this->add_boundary_node(2, node_pt);

				this->setup_boundary_element_info();
			}
};

template<class EL> class TwoMeshRefineableUnsteadyHeatProblem : public Problem {
	private:
		const uint Nx1;
		const uint Nx2;
		const double x0;
		const double x1;
		const double x2;
		UnsteadyHeatEquations<1>::UnsteadyHeatSourceFctPt source_pt;

		uint nbound;
		RefineableTwoLayerMesh<EL> *bulk_mesh_pt;
		Mesh *surf_mesh_pt;

		ConstitutiveLaw *const_law_pt;

	public:
		TwoMeshRefineableUnsteadyHeatProblem(const uint nxs[2], const double xs[3], UnsteadyHeatEquations<1>::UnsteadyHeatSourceFctPt source_inp)
		 : Nx1(nxs[0]), Nx2(nxs[1]), x0(xs[0]), x1(xs[1]), x2(xs[2]), source_pt(source_inp)
		{
			add_time_stepper_pt(new BDF<2>);

			bulk_mesh_pt = new RefineableTwoLayerMesh<EL>(nxs, xs, time_stepper_pt());
			surf_mesh_pt = new Mesh;

			add_sub_mesh(bulk_mesh_pt);
			add_sub_mesh(surf_mesh_pt);

			build_global_mesh();

			nbound = bulk_mesh_pt->nboundary();
			for (uint b = 0; b < nbound-1; b++) {
				uint nnode = bulk_mesh_pt->nboundary_node(b);
				for (uint n = 0; n < nnode; n++) {
					bulk_mesh_pt->boundary_node_pt(b, n)->pin(0);
					bulk_mesh_pt->boundary_node_pt(b, n)->set_value(0, (double)((int)b)-1);

					cout << "Value at b = " << b << " : " << bulk_mesh_pt->boundary_node_pt(b, n)[0] << endl;
				}
			}

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

			// const_law_pt = new GeneralisedHookean(&ProblemParams::nu);

			// uint nelem = Nx1 + Nx2;

			// for (uint e = 0; e < nelem; e++)
				// EL *el_pt = dynamic_cast<EL*>(bulk_mesh_pt->element_pt(e))->constitutive_law_pt() = const_law_pt;

			cout << "Number of equations : " << assign_eqn_numbers() << endl;
		}

		void actions_before_newton_solve() {}
		void actions_after_newton_solve() {}

		void actions_before_implicit_timestep() {
			double t = time_pt()->time();

			vecd x(1);
			vecd u(1);

			for (uint b = 0; b < nbound; b++) {
				int nnode = bulk_mesh_pt->nboundary_node(b);
				for (int n = 0; n < nnode; n++) {
					Node *node_pt = bulk_mesh_pt->boundary_node_pt(b, n);
					x[0] = node_pt->x(0);

					ProblemParams::get_exact_u(t, x, u);
					node_pt->set_value(0, u[0]);
				}
			}
		}

		void actions_after_implicit_timestep() {}

		void actions_before_adapt() {

		}

		void actions_after_adapt() {
			
		}

		void set_initial_condition() {
			vecd x(1);
			vecd u(1);

			for (uint e = 0; e < Nx1; e++) {
				EL *el_pt = dynamic_cast<EL*>(bulk_mesh_pt->element_pt(e));
				uint nnode = el_pt->nnode();
				for (uint n = 0; n < nnode; n++)
					el_pt->node_pt(n)->set_value(0, ProblemParams::Ts);
			}

			for (uint e = Nx1; e < Nx1+Nx2; e++) {
				EL *el_pt = dynamic_cast<EL*>(bulk_mesh_pt->element_pt(e));
				uint nnode = el_pt->nnode();
				for (uint n = 0; n < nnode; n++)
					el_pt->node_pt(n)->set_value(0, ProblemParams::Tfr);
			}

			assign_initial_values_impulsive();
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
};

int main() {
	const uint Nxs[2] = {100, 100};
	const double Xs[3] = {-1.0, 0.0, 1.0};

	TwoMeshRefineableUnsteadyHeatProblem<RefineableQUnsteadyHeatElement<1,2>> problem(Nxs, Xs, ProblemParams::get_source);
	//  QUnsteadyHeatElement<1,2>> 

	DocInfo info;
	info.set_directory("RESLT");
	info.number() = 0;

	const double t_max = 1.0;
	const double dt = 0.01;

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