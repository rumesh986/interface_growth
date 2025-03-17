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

	double mu_rat = 0.1;
	double rho_rat = 0.5;

	double Ca = 0.01;

	Vector<double> G(2);

	double Nu = 0.1;
}

template<class EL> class ElasticRefineableTwoLayerMesh : public virtual ElasticRefineableRectangularQuadMesh<EL> {
	public:
		ElasticRefineableTwoLayerMesh(
			const unsigned int &nx, const unsigned int &ny1, const unsigned int &ny2, 
			const double &lx, const double &h1, const double &h2, 
			const bool &periodicity, 
			TimeStepper *ts_ptr = &Mesh::Default_TimeStepper) :
				RectangularQuadMesh<EL>(nx, ny1+ny2, lx, h1+h2, periodicity, ts_ptr),
				ElasticRectangularQuadMesh<EL>(nx, ny1+ny2, lx, h1+h2, periodicity, ts_ptr),
				ElasticRefineableRectangularQuadMesh<EL>(nx, ny1+ny2, lx, h1+h2, periodicity, ts_ptr) {
					this->set_nboundary(5);

					for (unsigned int e = 0; e < nx; e++) {
						FiniteElement *el_ptr = this->finite_element_pt(nx*(ny1-1)+e);

						const unsigned n_node = el_ptr->nnode();

						for (unsigned int n = n_node-1; n > n_node-4; n--) {
							Node *node_ptr = el_ptr->node_pt(n);
							this->convert_to_boundary_node(node_ptr);
							this->add_boundary_node(4, node_ptr);
						}
					}

					this->setup_boundary_element_info();
				}
};

template<class ELEMENT, class TIMESTEPPER> class InterfaceProblem : public Problem {
	public:
		InterfaceProblem(const unsigned int &Nx, const unsigned int &Ny1, const unsigned int &Ny2, const double &Lx, const double &h1, const double &h2);
		~InterfaceProblem() {}

		void set_initial_condition();
		void set_boundary_condition();
		void doc_solution(DocInfo &doc_info);

		void unsteady_run(const double &t_max, const double &dt);

	private:
		void actions_before_newton_solve() {}
		void actions_after_newton_solve() {}

		void actions_before_adapt();
		void actions_after_adapt();

		void create_interface_elements();
		void delete_interface_elements();

		void fix_pressure(const unsigned int &e, const unsigned int &dof_ptr, const double &value_ptr) {
			dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->fix_pressure(dof_ptr, value_ptr);
		}

		void actions_before_implicit_timestep() {
			bulk_mesh_ptr->set_lagrangian_nodal_coordinates();
		}

		void deform_free_surface(const double &epsilon, const unsigned int &n_periods);

		ElasticRefineableRectangularQuadMesh<ELEMENT> *bulk_mesh_ptr;
		Mesh *surf_mesh_ptr;

		ConstitutiveLaw *const_law_ptr;

		double Lx;

		ofstream tracef;
};

template<class ELEMENT, class TIMESTEPPER> InterfaceProblem<ELEMENT, TIMESTEPPER>::InterfaceProblem(
	const unsigned int &Nx, const unsigned int &Ny1, const unsigned int &Ny2, const double &_Lx, const double &h1, const double &h2) : Lx(_Lx) {
		add_time_stepper_pt(new TIMESTEPPER);

		bulk_mesh_ptr = new ElasticRefineableTwoLayerMesh<ELEMENT>(Nx, Ny1, Ny2, Lx, h1, h2, true, time_stepper_pt());
		bulk_mesh_ptr->spatial_error_estimator_pt() = new Z2ErrorEstimator;
		bulk_mesh_ptr->max_refinement_level() = 4;

		surf_mesh_ptr = new Mesh;

		create_interface_elements();

		add_sub_mesh(bulk_mesh_ptr);
		add_sub_mesh(surf_mesh_ptr);

		build_global_mesh();

		const unsigned int n_bound = bulk_mesh_ptr->nboundary();

		for (unsigned int b = 0; b < n_bound; b++) {
			const unsigned int n_node = bulk_mesh_ptr->nboundary_node(b);

			for (unsigned int n = 0; n < n_node; n++) {
				if (b == 0 || b == 2) {
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

		const unsigned int n_lower = Nx * Ny1;
		const unsigned int n_bulk_elem = bulk_mesh_ptr->nelement();

		for (unsigned int e = 0; e < n_lower; e++) {
			ELEMENT *elem_ptr = dynamic_cast<ELEMENT*>(bulk_mesh_ptr->element_pt(e));

			elem_ptr->re_pt() = &Global_Params::Re;
			elem_ptr->re_st_pt() = &Global_Params::ReSt;
			elem_ptr->re_invfr_pt() = &Global_Params::ReInvFr;
			elem_ptr->g_pt() = &Global_Params::G;
			
			elem_ptr->constitutive_law_pt() = const_law_ptr;
		}

		for (unsigned int e = n_lower; e < n_bulk_elem; e++) {
			ELEMENT *elem_ptr = dynamic_cast<ELEMENT*>(bulk_mesh_ptr->element_pt(e));

			elem_ptr->re_pt() = &Global_Params::Re;
			elem_ptr->re_st_pt() = &Global_Params::ReSt;
			elem_ptr->re_invfr_pt() = &Global_Params::ReInvFr;
			elem_ptr->g_pt() = &Global_Params::G;

			elem_ptr->viscosity_ratio_pt() = &Global_Params::mu_rat;
			elem_ptr->density_ratio_pt() = &Global_Params::rho_rat;
			
			elem_ptr->constitutive_law_pt() = const_law_ptr;
		}

		fix_pressure(0, 0, 0.0);

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

template<class EL, class TS> void InterfaceProblem<EL,TS>::actions_before_adapt() {
	delete_interface_elements();
	rebuild_global_mesh();
}

template<class EL, class TS> void InterfaceProblem<EL,TS>::actions_after_adapt() {
	create_interface_elements();
	rebuild_global_mesh();

	const unsigned int n_node = bulk_mesh_ptr->nnode();
	for (unsigned int n = 0; n < n_node; n++)
		bulk_mesh_ptr->node_pt(n)->pin_position(0);
	
	RefineableNavierStokesEquations<2>::unpin_all_pressure_dofs(bulk_mesh_ptr->element_pt());
	RefineableNavierStokesEquations<2>::pin_redundant_nodal_pressures(bulk_mesh_ptr->element_pt());

	fix_pressure(0, 0, 0.0);

	PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(bulk_mesh_ptr->element_pt());

	set_boundary_condition();
}

template<class EL, class TS> void InterfaceProblem<EL, TS>::create_interface_elements() {
	const unsigned int n_elem = bulk_mesh_ptr->nboundary_element(4);
	for (unsigned int e = 0; e < n_elem; e++) {
		EL *bulk_elem_ptr = dynamic_cast<EL*>(bulk_mesh_ptr->boundary_element_pt(4, e));

		if (bulk_elem_ptr->viscosity_ratio_pt() != &Global_Params::mu_rat) {
			const int face_i = bulk_mesh_ptr->face_index_at_boundary(4, e);

			FiniteElement *if_elem_ptr = new ElasticLineFluidInterfaceElement<EL>(bulk_elem_ptr, face_i);
			dynamic_cast<ElasticLineFluidInterfaceElement<EL>*>(if_elem_ptr)->st_pt() = &Global_Params::St;
			dynamic_cast<ElasticLineFluidInterfaceElement<EL>*>(if_elem_ptr)->ca_pt() = &Global_Params::Ca;
			surf_mesh_ptr->add_element_pt(if_elem_ptr);
		}
	}
}

template<class EL, class TS> void InterfaceProblem<EL,TS>::delete_interface_elements() {
	const unsigned int n_if_elem = surf_mesh_ptr->nelement();
	for (unsigned int e = 0; e < n_if_elem; e++)
		delete surf_mesh_ptr->element_pt(e);
	
	surf_mesh_ptr->flush_element_and_node_storage();
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
	bulk_mesh_ptr->output_paraview(outfile, npts);
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

	unsigned int max_adapt = 2;
	for (unsigned int i = 0; i < 2; i++)
		refine_uniformly();
	
	doc_solution(info);
	info.number()++;

	const unsigned n_timestep = unsigned(t_max / dt);
	bool first_timestep = true;

	for (unsigned int t = 1; t < n_timestep; t++) {
		cout <<"\nTimestep " << t << "/" << n_timestep << endl;

		unsteady_newton_solve(dt, max_adapt, first_timestep);

		first_timestep = false;
		max_adapt = 1;

		doc_solution(info);
		info.number()++;
	}
}

int main(int argc, char **argv) {
	CommandLineArgs::setup(argc, argv);

	Global_Params::ReSt = Global_Params::Re * Global_Params::St;

	double t_max = 0.6;

	const double dt = 0.0025;

	const unsigned int Nx = 3;
	const unsigned int Ny1 = 3;
	const unsigned int Ny2 = 3;

	const double Lx = 1.0;
	const double h1 = 1.0;
	const double h2 = 1.0;

	Global_Params::G[0] = 0.0;
	Global_Params::G[1] = -1.0;

	InterfaceProblem<RefineablePseudoSolidNodeUpdateElement<RefineableQCrouzeixRaviartElement<2>, RefineableQPVDElement<2,3>>, BDF<2>> problem(Nx, Ny1, Ny2, Lx, h1, h2);

	problem.unsteady_run(t_max, dt);
}