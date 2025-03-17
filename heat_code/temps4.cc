#include <cmath>

#include "generic.h"
#include "solid.h"
#include "meshes/one_d_mesh.h"
#include "unsteady_heat.h"

using namespace std;
using namespace oomph;
using namespace MathematicalConstants;

typedef Vector<double> vecd;

namespace Params {
	const double Ts = -10.0;
	const double T_fr = 10.0;

	const double lambda_s = 1.0;
	const double lambda_i = 2.0;

	const double Ds = 10.0;
	const double Di = 20.0;

	double e_rat = sqrt(lambda_s/lambda_i);


	double T0(double t) {
		const double eps = 1e-6;
		return Ts * erf(eps/(2*sqrt(Ds*t)));
	}

	void get_exact_u(const double &t, const vecd &x, vecd &u) {
		double t0 = T0(t);
		double D = (x[0] < 0) ? Ds : Di;
		double e = (x[0] < 0) ? 1 : e_rat;

		u[0] = t0 + (t0 - Ts) * erf(x[0]/(2*sqrt(D*t))) * e;
	}

	void get_source(const double &t, const vecd &x, double &source) {
		source = 0.0;
	}
}

template<class EL> class ElasticRefineableTwoLayerMesh : public virtual Elas

template<class EL> class UnsteadyHeat1D2Mesh : public Problem {
	public:
		UnsteadyHeat1D2Mesh(const uint Nx1, const uint Nx2, const double Lx1, const double Lx2, UnsteadyHeatEquations<1>::UnsteadyHeatSourceFctPt inp) : source_pt(inp) {
			add_time_stepper_pt(new BDF<2>);

			mesh_pt_neg = new OneDMesh<EL>(Nx1, Lx1);
			mesh_pt_pos = new OneDMesh<EL>(Nx2, Lx2);
		}

	private:
		UnsteadyHeatEquations<1>::UnsteadyHeatSourceFctPt source_pt;
}

int main() {
	cout << "Hello World!" << endl;
}