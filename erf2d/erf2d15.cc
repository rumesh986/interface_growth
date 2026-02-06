// moving interface
// calculate as part of solution in newton iterations (hopefully)
// move with spine meshes
// using proper non-dimensionalized form with non-dimensionalized constants
// temperature scaling is changed to non-dim temps have same sign as dimensional temps

#include <filesystem>

#include "includes.h"
#include "two_phase_free_boundary_mesh.h"

#ifndef RUN_SCRIPT
#define X_ORDER 2
#define T_ORDER 1
#endif

// dimensionless parameters
double k = 1.0;
double D = 1.0;
double De = 1.0;
double St = 1.0;

bool ic_set = false;

// boundary condition
double Tl = 1.0;

// simulation domain
double xs[2] = {0.0, 3.0};
double ys[2] = {0.0, 0.0001};

FreeBoundaryElement *geom_element;

void no_flux_fct(const double &t, const Vector<double> &x, double &flux) {
	flux = 0.0;
}

void get_source(const double &t, const Vector<double> &x, double &source) {
	source = 0.0;
}

// analytical solution
// u[0] -> solution with interface in current position (from numerical solution)
// u[1] -> solution with interface position from analytical solution
void get_exact_u(const double &t, const Vector<double> &x, Vector<double> &u) {
	double h = geom_element->get_interface();
	double h_ana = sqrt(De * t);
	
	// get values based on numerical position of interface
	if (x[0] < h) {
		double denom = 0.5 / sqrt(t);
		u[0] = erf(x[0] * denom) / erf(h * denom) - 1.0;
	} else {
		double denom = 0.5 / sqrt(D * t);
		u[0] = Tl * (erf(x[0] * denom) - erf(h * denom)) / (1.0 - erf(h * denom));
	}

	// get values based on analytical position of interface
	if (x[0] < h_ana) {
		double denom = 0.5 / sqrt(t);
		u[1] = erf(x[0] * denom) / erf(h_ana * denom) - 1.0;
	} else {
		double denom = 0.5 / sqrt(D * t);
		u[1] = Tl * (erf(x[0] * denom) - erf(h_ana * denom)) / (1.0 - erf(h_ana * denom));
	}
}

Vector<unsigned int> analytical_boundaries = {
	1, // right
	3, // left
};

Vector<unsigned int> pinned_boundaries = {
	// 3, // left
	4 // interface
};

map<unsigned int, FluxFctPt> flux_boundaries = {
	// {0, no_flux_fct},
	// {2, no_flux_fct},
};

template<class EL>
class Erf2DProblem : public Problem {
	private:
		unsigned int nx1, nx2, nx, ny;
		unsigned int t_steps;
		double dt, t_shift;

		DocInfo info;

		TwoPhaseFreeBoundarySpineMesh<SpineElement<EL>> *bulk_mesh_pt;
		Mesh *geometry_mesh_pt;

	public:
		Erf2DProblem(
			unsigned int nx1_,
			unsigned int nx2_,
			unsigned int ny_,
			unsigned int tsteps_,
			double dt_,
			double tshift_,
			DocInfo info_
		) : nx1(nx1_), nx2(nx2_), nx(nx1_+nx2_), ny(ny_), t_steps(tsteps_), dt(dt_), t_shift(tshift_), info(info_) {

			add_time_stepper_pt(new BDF<T_ORDER>);

			geom_element = new FreeBoundaryElement(xs[0], sqrt(De * t_shift), xs[1], ys[0], ys[1], St, time_stepper_pt());
			printf("x0=%8.6f x1=%8.6f x2=%8.6f\n", geom_element->x0(), geom_element->x1(), geom_element->x2());
		
			bulk_mesh_pt = new TwoPhaseFreeBoundarySpineMesh<SpineElement<EL>>(nx1, nx2, ny, geom_element, time_stepper_pt());
			bulk_mesh_pt->setup_boundary_element_info();
			add_sub_mesh(bulk_mesh_pt);

			geometry_mesh_pt = new Mesh;
			geometry_mesh_pt->add_element_pt(geom_element);
			add_sub_mesh(geometry_mesh_pt);

			build_global_mesh();

			for (unsigned int e = 0; e < bulk_mesh_pt->nelement(); e++)
				dynamic_cast<EL *>(bulk_mesh_pt->element_pt(e))->source_fct_pt() = get_source;
			
			for (unsigned int b : pinned_boundaries) {
				unsigned long int nnode = bulk_mesh_pt->nboundary_node(b);
				for (unsigned long int n = 0; n < nnode; n++)
					bulk_mesh_pt->boundary_node_pt(b, n)->pin_all();
			}

			for (unsigned int b : analytical_boundaries) {
				unsigned long int nnode = bulk_mesh_pt->nboundary_node(b);
				for (unsigned long int n = 0; n < nnode; n++)
					bulk_mesh_pt->boundary_node_pt(b, n)->pin_all();
			}

			for (auto iter : flux_boundaries) {
				create_flux_elements(iter.first, iter.second);
			}

			for (unsigned int yi = 0; yi < ny; yi++) {
				unsigned int base = yi * nx;
				for (unsigned int e = nx1; e < nx; e++)
					dynamic_cast<EL *>(bulk_mesh_pt->element_pt(base + e))->beta_pt() = &D;
			}

			printf("Total number of equations: %lu\n", assign_eqn_numbers());
			printf("NDOF: %lu\n", ndof());

			linear_solver_pt()->disable_doc_time();
			disable_info_in_newton_solve();
			// newton_solver_tolerance() = 1e-9;
			max_newton_iterations() = 1e2;
			max_residuals() = 1e3;
		}

		~Erf2DProblem() {
			delete bulk_mesh_pt;
			delete geometry_mesh_pt;
		}

		// Neumann boundaries need special flux elements associated at the boundaries
		// this method will create the required elements and attach the correct flux functions to them
		void create_flux_elements(unsigned int b, FluxFctPt &flux_pt) {
			unsigned int nelems = bulk_mesh_pt->nboundary_element(b);
			for (unsigned int e = 0; e < nelems; e++) {
				EL *elem = dynamic_cast<EL *>(bulk_mesh_pt->boundary_element_pt(b, e));
				int face_i = bulk_mesh_pt->face_index_at_boundary(b, e);
				UnsteadyHeatFluxElement<EL> *flux_elem = new UnsteadyHeatFluxElement<EL>(elem, face_i);
				flux_elem->flux_fct_pt() = flux_pt;
				bulk_mesh_pt->add_element_pt(flux_elem);
			}
		}

		// update boundary values to match (semi) infinite domain
		void actions_before_implicit_timestep() {
			Vector<double> x(2), u(2);

			double time = time_pt()->time() + dt;

			for (unsigned int b : analytical_boundaries) {
				unsigned long int nnode = bulk_mesh_pt->nboundary_node(b);
				for (unsigned long int n = 0; n < nnode; n++) {
					bulk_mesh_pt->boundary_node_pt(b, n)->position(x);
					get_exact_u(time, x, u);
					bulk_mesh_pt->boundary_node_pt(b,n)->set_value(0, u[1]);
				}
			}
		}
		void actions_after_implicit_timestep() {};

		// update flux at each newton step
		void actions_before_newton_step() {
			Vector<double> flux(2), flux_fd(2);

			double tot_flux = 0.0;
			unsigned long int nelems = bulk_mesh_pt->nboundary_element(4);

			Vector<double> s(2);
			s[1] = 0.0;
			
			for (unsigned long int e = 0; e < nelems; e++) {
				int face_index = bulk_mesh_pt->face_index_at_boundary(4, e);
				EL *elem = dynamic_cast<EL *>(bulk_mesh_pt->boundary_element_pt(4, e));
				
				double factor = 0.0;
				if (face_index == 1) {
					s[0] = 1.0;
					factor = 1.0;
				} else if (face_index == -1) {
					s[0] = -1.0;
					factor = -k;
				}

				elem->get_flux(s, flux);
				tot_flux += factor * flux[0];
			}
			
			printf("Setting total flux to %8.6f interface at %16.14f\n", tot_flux / ny, geom_element->get_interface());
			geom_element->set_flux(tot_flux / ny);
		}

		// set new interface position after each newton step
		void actions_after_newton_step() {
			for (unsigned s = 0; s < bulk_mesh_pt->nspine(); s++) {
				bulk_mesh_pt->spine_pt(s)->height() = geom_element->x1();
			}

			bulk_mesh_pt->node_update();
		}

		void set_initial_condition() {
			time_pt()->time() = t_shift;
			unsigned long int nnode = bulk_mesh_pt->nnode();

			Vector<double> x(2);
			Vector<double> u(1);

			unsigned int tsteps = time_stepper_pt()->nprev_values();

			for (unsigned long int n = 0; n < nnode; n++)
				time_stepper_pt()->assign_initial_positions_impulsive(bulk_mesh_pt->node_pt(n));

			for (unsigned long int n = 0; n < nnode; n++) {
				bulk_mesh_pt->node_pt(n)->position(tsteps, x);
				get_exact_u(time_pt()->time(tsteps), x, u);
				bulk_mesh_pt->node_pt(n)->set_value(tsteps, 0, u[1]);
			}

			unsigned int step = 0;
			doc_step(step, tsteps);
			step++;

			unsigned long int nelems = bulk_mesh_pt->nboundary_element(4);

			for (int t = tsteps-1; t >= 0; t--) {
				double time = time_pt()->time((unsigned int) t);

				// analytical interface position
				double new_h = sqrt(De * time);

				printf("IC %u: x_old=%8.6f x_new=%8.6f\n", t, geom_element->get_interface(), new_h);
				geom_element->set_interface(new_h);
				bulk_mesh_pt->node_update();

				for (unsigned long int n = 0; n < nnode; n++) {
					bulk_mesh_pt->node_pt(n)->position(t, x);
					get_exact_u(time, x, u);
					bulk_mesh_pt->node_pt(n)->set_value(t, 0, u[1]);
				}

				printf("[% 4d] Setting initial condition at t=%8.6f\n", step, time);
				doc_step(step, t);
				step++;
			}

			time_pt()->time() = t_shift;
			ic_set = true;
		}

		void doc_step(const unsigned int &timestep, const unsigned int &t = 0) {
			double time = time_pt()->time(t);
			unsigned long int nnode = bulk_mesh_pt->nnode();

			Vector<double> x(2);
			Vector<double> exact_u(2);
			Vector<double> numerical_u(1);

			double tot_error = 0.0; // L2 norm of errors

			char fname[256];
			sprintf(fname, "%s/steps/step%u.dat", info.directory().c_str(), info.number());
			FILE *file = fopen(fname, "w");
			for (unsigned long int n = 0; n < nnode; n++) {
				bulk_mesh_pt->node_pt(n)->position(t, x);
				bulk_mesh_pt->node_pt(n)->value(t, numerical_u);
				get_exact_u(time, x, exact_u);

				double error = numerical_u[0] - exact_u[1]; // error in each node
				tot_error += error * error;

				fprintf(file, "%16.14f %16.14f %16.14f %16.14f %16.14f\n", x[0], x[1], exact_u[1], numerical_u[0], error);
			}
			fclose(file);

			tot_error = sqrt(tot_error) / nnode;

			// get element size for error analysis graphs
			double max_elem_size, min_elem_size;
			bulk_mesh_pt->max_and_min_element_size(max_elem_size, min_elem_size);

			sprintf(fname, "%s/results.dat", info.directory().c_str());
			file = fopen(fname, "a");
			fprintf(file, "%16.14f %16.14f %16.14f %16.14f %16.14f\n", time, tot_error, geom_element->get_interface(), sqrt(De*time), max_elem_size/ys[1]);
			fclose(file);

			printf("[%4u] time=%8.6f error=%e iface_err=%e interface=%8.6f expected=%8.6f\n", timestep, time, tot_error, fabs(geom_element->get_interface() - sqrt(De * time)), geom_element->get_interface(), sqrt(De*time));
			info.number()++;
		}
};

int main(int argc, char **argv) {
	CommandLineArgs::setup(argc, argv);

	unsigned int nx1 = 0;
	unsigned int nx2 = 0;
	unsigned int ny = 1;
	unsigned int t_steps = 100;
	double dt = 0.0;
	double t_shift = 0.0;
	unsigned int wf = 1;
	std::string dname;

	CommandLineArgs::specify_command_line_flag("--nx", &nx1);
	CommandLineArgs::specify_command_line_flag("--ny", &ny);
	CommandLineArgs::specify_command_line_flag("--tsteps", &t_steps);
	CommandLineArgs::specify_command_line_flag("--dt", &dt);
	CommandLineArgs::specify_command_line_flag("--tshift", &t_shift);
	CommandLineArgs::specify_command_line_flag("--vardt");
	CommandLineArgs::specify_command_line_flag("--write-freq", &wf);
	CommandLineArgs::specify_command_line_flag("--k", &k);
	CommandLineArgs::specify_command_line_flag("--D", &D);
	CommandLineArgs::specify_command_line_flag("--De", &De);
	CommandLineArgs::specify_command_line_flag("--St", &St);
	CommandLineArgs::specify_command_line_flag("--dname", &dname, "doc");

	CommandLineArgs::parse_and_assign();

	if (nx1 == 0) {
		cout << "Error: Nx not specified" << endl;
		exit(1);
	}

	if (dt == 0.0) {
		cout << "Error: dt not specified" << endl;
		exit(1);
	}

	if (!CommandLineArgs::command_line_flag_has_been_set("--dname")) {
		char temp[256];
		sprintf(temp, "RESLT/%dn%u_%dt%.2e", X_ORDER, nx1, T_ORDER, dt);
		dname.assign(temp);
	}

	printf("Saving results to %s\n", dname.c_str());
	
	if (std::filesystem::exists(dname.c_str())) {
		cout << dname << " exists" << endl;
	} else {
		std::filesystem::create_directories(dname.c_str());
	}

	char sub_dname[1024];
	sprintf(sub_dname, "%s/steps", dname.c_str());
	if (std::filesystem::exists(sub_dname)) {
		cout << sub_dname << " exists" << endl;
	} else {
		std::filesystem::create_directories(sub_dname);
	}

	DocInfo info;
	info.set_directory(dname.c_str());
	info.number() = 0;

	cout << "Output directory: " << info.directory() << endl;

	if (t_shift == 0.0) {
		t_shift = T_ORDER * dt;
	}

	printf("Problem Def:\n");
	printf("\tk=%8.6f D=%8.6f, St=%8.6f\n", k, D, St);
	printf("\tDe=%e\n", De);
	printf("\tx0=%8.6f x2=%8.6f\n", xs[0], xs[1]);

	nx2 = nx1;
	auto problem = Erf2DProblem<QUnsteadyHeatElement<2, X_ORDER>>(nx1, nx2, ny, t_steps, dt, t_shift, info);

	problem.initialise_dt(dt);
	problem.set_initial_condition();

	char config_fname[256];
	sprintf(config_fname, "%s/config", dname.c_str());

	FILE *file = fopen(config_fname, "w");
	fprintf(file, "nx1=%d\n", nx1);
	fprintf(file, "nx2=%d\n", nx2);
	fprintf(file, "nx=%d\n", nx1+nx2);
	fprintf(file, "x_order=%d\n", X_ORDER);
	fprintf(file, "dt=%e\n", dt);
	fprintf(file, "t_shift=%e\n", t_shift);
	fprintf(file, "t_order=%d\n", T_ORDER);
	fprintf(file, "t_steps=%d\n", t_steps);
	fprintf(file, "write_freq=%u\n", wf);
	fprintf(file, "k=%10.8f\n", k);
	fprintf(file, "D=%e\n", D);
	fprintf(file, "De=%e\n", De);
	fclose(file);

	int prev_steps = problem.time_stepper_pt()->nprev_values() + 1;

	for (unsigned int t = 0; t < t_steps; t++) {
		problem.unsteady_newton_solve(dt);

		// write data to file
		if (t % wf == 0 || t == t_steps-1) 
			problem.doc_step(t+prev_steps);
		
		double h = geom_element->get_interface();

		// exit if interface is leaving domain or has problems
		if (h >= xs[1] || isnan(h)) {
			printf("Interface reached right boundary, exiting...\n");
			break;
		}
	}

}