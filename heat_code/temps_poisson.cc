#include <cmath>

#include "generic.h"
#include "solid.h"
#include "poisson.h"
#include "meshes/one_d_mesh.h"
// #include "unsteady_heat.h"

using namespace std;
using namespace oomph;
using namespace MathematicalConstants;

typedef Vector<double> vecd;

namespace SteadyParams {
	const double alpha = 0.1;

	void get_exact_u(const vecd& x, vecd& u) {
		u[0] = alpha*x[0]*x[0] + (1-alpha)*x[0];
	}

	// void source_func(const double &x, double &source) {
	void source_func(const vecd &x, double &source) {
		source = alpha;
	}	
}

template <class EL> class PoissonProblem : public Problem {
	public:
		PoissonProblem(const double Lx, const uint& n_elem, PoissonEquations<1>::PoissonSourceFctPt source_pt_inp) : source_pt(source_pt_inp) {
			mesh_pt() = new OneDMesh<EL>(n_elem, Lx);
			
			uint nbound = mesh_pt()->nboundary();
			for (uint i = 0; i < nbound; i++) 
				mesh_pt()->boundary_node_pt(i, 0)->pin(0);
			
			// for (uint i = 0; i < n_elem; i++) {
				// EL *elem_pt = dynamic_cast<EL*>(mesh_pt()->element_pt(i));
				// elem_pt->source_fct_pt() = source_pt;
			// }
			for (uint i = 0; i < n_elem; i++)
				dynamic_cast<EL*>(mesh_pt()->element_pt(i))->source_fct_pt() = source_pt;

			assign_eqn_numbers();
		}

		~PoissonProblem() {
			delete mesh_pt();
		}

		void actions_before_newton_solve() {
			
			// double x = bound_node_pt->x(0);
			// double u;
			vecd x(1), u(1);

			// uint nbound = mesh_pt()->nboundary();
			// for (uint i = 0; i < nbound; i++) {
			// 	x[0] = mesh_pt()->boundary_node_pt(i, 0)->x(0);
			// 	SteadyParams::get_exact_u(x, u);
			// 	mesh_pt()->boundary_node_pt(i, 0)->set_value(0, u[0]);
			// }

			Node *bound_node_pt = mesh_pt()->node_pt(0);
			x[0] = bound_node_pt->x(0);

			SteadyParams::get_exact_u(x, u);
			bound_node_pt->set_value(0, u[0]);

			bound_node_pt = mesh_pt()->node_pt(mesh_pt()->nnode()-1);
			x[0] = bound_node_pt->x(0);
			
			SteadyParams::get_exact_u(x, u);
			bound_node_pt->set_value(0, u[0]);
		}
		
		void actions_after_newton_solve() {}

		void doc_solution(DocInfo &info) {
			using namespace StringConversion;

			const uint npts = 5;

			char filename[100];
			ofstream outfile;

			sprintf(filename, "%s/soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			mesh_pt()->output(outfile, npts);
			outfile.close();

			sprintf(filename, "%s/exact_soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			mesh_pt()->output_fct(outfile, npts, SteadyParams::get_exact_u);
			outfile.close();

			double error, norm;
			sprintf(filename, "%s/error%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			mesh_pt()->compute_error(outfile, SteadyParams::get_exact_u, error, norm);
			outfile.close();

			cout << "\nNorm of error	: " << sqrt(error) << endl;
			cout << "Norm of solution 	: " << sqrt(norm) << endl;
			cout << endl;

			info.number()++;
		}
	
	private:
		PoissonEquations<1>::PoissonSourceFctPt source_pt;
};

int main() {
	const double Lx = 1.0;
	const uint n_elem = 5000;
	PoissonProblem<QPoissonElement<1, 4>> problem(Lx, n_elem, SteadyParams::source_func);

	cout << "\n\nProblem self-test";
	if (problem.self_test() == 0)
		cout << "Passed: Problem can be solved" << endl;
	else
		throw OomphLibError("Failed!", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
	
	problem.newton_solve();
	DocInfo info;
	info.set_directory("RESLT");
	info.number() = 0;

	problem.doc_solution(info);
}