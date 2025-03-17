#include "generic.h"
#include "fluid_interface.h"
#include "solid.h"
#include "constitutive.h"
#include "navier_stokes.h"

#include "meshes/rectangular_quadmesh.h"

using namespace std;
using namespace oomph;

namespace Global_Params {
	double Re = 5.0;
	double St = 1.0;

	double ReSt;

	double ReInvFr = 5.0;

	double Ca = 0.01;

	Vector<double> G(2);

	double Nu = 0.1;
}

template<class ELEMENT, class TIMESTEPPER> class InterfaceProblem : public Problem {
	public:
		InterfaceProblem(const unsigned int &Nx, const unsigned int &Ny, const double &Lx, const double &h);
		~InterfaceProblem() {}

		void set_initial_condition();
		void set_boundary_condition();
		void doc_solution(DocInfo &doc_info);

		void unsteady_run(const double &t_max, const double &dt);

	private:
		void actions_before_newton_solve() {}
		void actions_after_newton_solve() {}

		void actions_before_implicit_timestep() {
			bulk_mesh_ptr->set_lagrangian_nodal_coordinates();
		}

		void deform_free_surface(const double &epsilon, const unsigned int &n_periods);

		ElasticRectangularQuadMesh<ELEMENT> *bulk_mesh_ptr;
		Mesh *surf_mesh_ptr;

		ConstitutiveLaw *const_law_ptr;

		double Lx;

		ofstream tracef;
};

template<class ELEMENT, class TIMESTEPPER> InterfaceProblem<ELEMENT, TIMESTEPPER>::InterfaceProblem(
	const unsigned int &Nx, const unsigned int &Ny, const double &_Lx, const double &h) : Lx(_Lx) {
		add_time_stepper_pt(new TIMESTEPPER);

		bulk_mesh_ptr = new ElasticRectangularQuadMesh<ELEMENT>(Nx, Ny, Lx, h, true, time_stepper_pt());
		surf_mesh_ptr = new Mesh;

		for (unsigned int e = 0; e < Nx; e++) {
			FiniteElement *bulk_elem_ptr = bulk_mesh_ptr->finite_element_pt(Nx*(Ny-1)+e);
			FiniteElement *if_elem_ptr = new ElasticLineFluidInterfaceElement<ELEMENT>(bulk_elem_ptr, 2);

			this->surf_mesh_ptr->add_element_pt(if_elem_ptr);
		}

		add_sub_mesh(bulk_mesh_ptr);
		add_sub_mesh(surf_mesh_ptr);

		build_global_mesh();

		const unsigned int n_bound = bulk_mesh_ptr->nboundary();

		for (unsigned int b= 0; b < n_bound; b++) {
			const unsigned int n_node = bulk_mesh_ptr->nboundary_node(b);

			for (unsigned int n = 0; n < n_node; n++) {
				if (b == 0) {
					bulk_mesh_ptr->boundary_node_pt(b, n)->pin(0);
					bulk_mesh_ptr->boundary_node_pt(b, n)->pin(1);

					bulk_mesh_ptr->boundary_node_pt(b,n)->pin_position(1);
				} else if (b == 1 || b == 3) {
					bulk_mesh_ptr->boundary_node_pt(b, n)->pin(0);
				}

			}
		}

		const unsigned int n_node = bulk_mesh_ptr->nnode();
		for (unsigned int n = 0; n < n_node; n++)
			bulk_mesh_ptr->node_pt(n)->pin_position(0);
		
		const_law_ptr = new GeneralisedHookean(&Global_Params::Nu);

		const unsigned int n_bulk_elem = bulk_mesh_ptr->nelement();

		for (unsigned int e = 0; e < n_bulk_elem; e++) {
			ELEMENT *elem_ptr = dynamic_cast<ELEMENT*>(bulk_mesh_ptr->element_pt(e));

			elem_ptr->re_pt() = &Global_Params::Re;
			elem_ptr->re_st_pt() = &Global_Params::ReSt;
			elem_ptr->re_invfr_pt() = &Global_Params::ReInvFr;
			elem_ptr->g_pt() = &Global_Params::G;
			
			elem_ptr->constitutive_law_pt() = const_law_ptr;

		}

		Data *P_ext_ptr = new Data(1);
		P_ext_ptr->pin(0);
		P_ext_ptr->set_value(0, 1.31);

		const unsigned int n_surf_elem = surf_mesh_ptr->nelement();

		for (unsigned int e = 0; e < n_surf_elem; e++) {
			ElasticLineFluidInterfaceElement<ELEMENT>* elem_ptr = dynamic_cast<ElasticLineFluidInterfaceElement<ELEMENT>*>(surf_mesh_ptr->element_pt(e));

			elem_ptr->st_pt() = &Global_Params::St;
			elem_ptr->ca_pt() = &Global_Params::Ca;

			elem_ptr->set_external_pressure_data(P_ext_ptr);
		}

		set_boundary_condition();

		cout << "Number of equations: " << assign_eqn_numbers() << endl;		
}

template <class ELEMENT, class TIMESTEPPER> void InterfaceProblem<ELEMENT, TIMESTEPPER>::set_initial_condition() {
	const unsigned int n_node = mesh_pt()->nnode();

	for (unsigned int n = 0; n < n_node; n++) {
		for (unsigned int i = 0; i < 2; i++) {
			mesh_pt()->node_pt(n)->set_value(i, 0.0);
		}
	}

	assign_initial_values_impulsive();
}

template<class EL, class TS> void InterfaceProblem<EL,TS>::set_boundary_condition() {
	const unsigned int n_bound = bulk_mesh_ptr->nboundary();

	for(unsigned int b = 0; b < n_bound; b++) {
		const unsigned int n_node = bulk_mesh_ptr->nboundary_node(b);
		for (unsigned int n = 0; n < n_node; n++) {
			if (b == 0)
				bulk_mesh_ptr->boundary_node_pt(b,n)->set_value(1, 0.0);
			
			if (b != 2)
				bulk_mesh_ptr->boundary_node_pt(b,n)->set_value(0, 0.0);
		}
	}
}

template<class EL, class TS> void InterfaceProblem<EL,TS>::deform_free_surface(const double &epsilon, const unsigned int &n_periods) {
	const unsigned int n_node = bulk_mesh_ptr->nnode();
	for (unsigned int n = 0; n < n_node; n++) {
		const double current_x_pos =bulk_mesh_ptr->node_pt(n)->x(0);
		const double current_y_pos =bulk_mesh_ptr->node_pt(n)->x(1);

		const double new_y_pos = current_y_pos + (1.0 - fabs(1.0-current_y_pos))*epsilon*cos(2.0*n_periods*MathematicalConstants::Pi*current_x_pos/Lx);

		bulk_mesh_ptr->node_pt(n)->x(1) = new_y_pos;
	}
}

template<class EL, class TS> void InterfaceProblem<EL,TS>::doc_solution(DocInfo &doc_info) {
	cout << "Time is now " << time_pt()->time() << endl;

	ElasticLineFluidInterfaceElement<EL> *elem_ptr = dynamic_cast<ElasticLineFluidInterfaceElement<EL>*>(surf_mesh_ptr->element_pt(0));

	tracef << time_pt()->time() << " " << elem_ptr->node_pt(0)->x(1) << endl;

	ofstream outfile;
	char filename[100];
	const unsigned int npts = 5;

	sprintf(filename, "%s/soln%i.dat", doc_info.directory().c_str(), doc_info.number());
	outfile.open(filename);
	bulk_mesh_ptr->output(outfile, npts);
	outfile.close();

	sprintf(filename, "%s/interface_soln%i.dat", doc_info.directory().c_str(), doc_info.number());
	outfile.open(filename);
	surf_mesh_ptr->output(outfile, npts);
	outfile.close();
}

template <class EL, class TS> void InterfaceProblem<EL,TS>::unsteady_run(const double &t_max, const double &dt) {
	const double epsilon = 0.1;
	const unsigned int n_periods = 1;

	deform_free_surface(epsilon, n_periods);

	DocInfo info;
	info.set_directory("RESLT");
	info.number() = 0;

	char filename[100];
	sprintf(filename, "%s/trace.dat", info.directory().c_str());
	tracef.open(filename);

	tracef << "time, free surface height" << endl;

	initialise_dt(dt);
	set_initial_condition();

	const unsigned n_timestep = unsigned(t_max / dt);
	doc_solution(info);
	info.number()++;

	for (unsigned int t = 1; t < n_timestep; t++) {
		cout <<"\nTimestep " << t << "/" << n_timestep << endl;

		unsteady_newton_solve(dt);
		doc_solution(info);
		info.number()++;
	}
}

int main(int argc, char **argv) {
	CommandLineArgs::setup(argc, argv);

	Global_Params::ReSt = Global_Params::Re * Global_Params::St;

	double t_max = 0.6;

	const double dt = 0.0025;

	const unsigned int Nx = 12;
	const unsigned int Ny = 12;

	const double Lx = 1.0;
	const double h = 1.0;

	Global_Params::G[0] = 0.0;
	Global_Params::G[1] = -1.0;

	InterfaceProblem<PseudoSolidNodeUpdateElement<QCrouzeixRaviartElement<2>, QPVDElement<2,3>>, BDF<2>> problem(Nx, Ny, Lx, h);

	problem.unsteady_run(t_max, dt);
}