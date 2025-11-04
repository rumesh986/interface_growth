
// moving interface
// calculate as part of solution in newton iterations (hopefully)
// move with spine meshes
// using proper non-dimensionalized form (without removing constants)

#include <filesystem>

#include "includes.h"
#include "two_phase_free_boundary_mesh.h"

#ifndef RUN_SCRIPT
#define X_ORDER 2
#define T_ORDER 1
#endif

double k[2] = {1.0, 0.25};
double rho[2] = {1.0, 2.0};
double Cp[2] = {1.0, 0.25};

double D[3] {k[0]/(Cp[0]*rho[0]), k[1]/(Cp[1]*rho[1]), 0.0};
double L = 1.0;
double St = 0.0;

bool ic_set = false;

static const double T_s = -1.0;
static const double T_m = 0.0;
static const double T_l = 1.0;
static const double T_p = (T_l - T_m) / (T_l - T_s);

double xs[3] = {0.0, 0.025, 3.0};
static const double ys[2] = {0.0, 0.0001};

FreeBoundaryElement *geometry;

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
	double h = 0.0;
	double h_ana = 0.0;
	if (ic_set) {
		h = (D[2] == 0.0) ? geometry->get_interface() : sqrt(D[2] * t);
		h_ana = sqrt(D[2] * t);
	} else {
		h = geometry->get_interface();
		h_ana = sqrt(D[2] * t);
	}

	if (x[0] < h) {
		double erf_zeta = erf(0.5 * x[0] / sqrt(D[0] * t));
		double erf_iface = erf(0.5 * h / sqrt(D[0]*t));
		u[0] = 1.0 + (T_p - 1.0) * erf_zeta/erf_iface;
	} else {
		double erf_zeta = erf(0.5  * x[0] / sqrt(D[1] * t));
		double erf_iface = erf(0.5 * h / sqrt(D[1]*t));
		u[0] = T_p * (1.0 - erf_zeta) / (1.0 - erf_iface);
	}

	if (x[0] < h_ana) {
		double erf_zeta = erf(0.5 * x[0] / sqrt(D[0] * t));
		double erf_iface = erf(0.5 * h_ana / sqrt(D[0]*t));
		u[1] = 1.0 + (T_p - 1.0) * erf_zeta/erf_iface;
	} else {
		double erf_zeta = erf(0.5  * x[0] / sqrt(D[1] * t));
		double erf_iface = erf(0.5 * h_ana / sqrt(D[1]*t));
		u[1] = T_p * (1.0 - erf_zeta) / (1.0 - erf_iface);
	}
}

Vector<unsigned int> analytical_boundaries = {
	1, // right
	3, // left
};

Vector<unsigned int> pinned_boundaries = {
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

		TimeStepper *iface_ts_pt;

		public:
		Erf2DProblem(
			unsigned int nx1_,
			unsigned int nx2_,
			unsigned int ny_,
			unsigned int tsteps,
			double dt_,
			double tshift,
			DocInfo info_
		) : nx1(nx1_), nx2(nx2_), nx(nx1_+nx2_), ny(ny_), t_steps(tsteps), dt(dt_), t_shift(tshift), info(info_) {

			add_time_stepper_pt(new BDF<T_ORDER>);

			iface_ts_pt = new BDF<1>;
			add_time_stepper_pt(iface_ts_pt);

			geometry = new FreeBoundaryElement(xs[0], xs[1], xs[2], ys[0], ys[1], St, iface_ts_pt);

			printf("x0=%8.6f x1=%8.6f x2=%8.6f\n", geometry->x0(), geometry->x1(), geometry->x2());

			bulk_mesh_pt = new TwoPhaseFreeBoundarySpineMesh<SpineElement<EL>>(nx1, nx2, ny, geometry, time_stepper_pt());
			bulk_mesh_pt->setup_boundary_element_info();
			add_sub_mesh(bulk_mesh_pt);

			geometry_mesh_pt = new Mesh;
			geometry_mesh_pt->add_element_pt(geometry);
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
				for (unsigned int e = 0; e < nx1; e++)
					dynamic_cast<EL *>(bulk_mesh_pt->element_pt(base + e))->beta_pt() = &D[0];
				for (unsigned int e = nx1; e < nx; e++)
					dynamic_cast<EL *>(bulk_mesh_pt->element_pt(base + e))->beta_pt() = &D[1];
			}

			printf("Total number of equations: %lu\n", assign_eqn_numbers());
			printf("NDOF: %lu\n", ndof());

			// for (unsigned int e = 0; e < nx*ny; e++) {
			// 	printf("[%u] ngeom=%u\n", e, dynamic_cast<SpineElement<EL>*>(bulk_mesh_pt->element_pt(e))->ngeom_data());
			// }

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

		// void actions_before_newton_solve() {};
		void actions_after_newton_solve() {};

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

		// void actions_before_newton_convergence_check() {
		void actions_before_newton_step() {
			// set all node values according to bottom row of nodes
			// for (unsigned int e = 0; e < nx; e++) {
			// 	FiniteElement *elem = bulk_mesh_pt->finite_element_pt(e);
			// 	// elem->node_pt(2)->set_value(0, elem->nodal_value(0, 0));
			// 	// elem->node_pt(3)->set_value(0, elem->nodal_value(1, 0));
			// 	for (unsigned int n1 = 0; n1 < elem->nnode_1d(); n1++) {
			// 		for (unsigned int n2 = 1; n2 < elem->nnode_1d(); n2++) {
			// 			elem->node_pt(n1 + elem->nnode_1d() * n2)->set_value(0, elem->nodal_value(n1, 0));
			// 		}
			// 	}
			// }
			
			Vector<double> flux(2), flux_fd(2);

			double tot_flux = 0.0;
			double tot_flux_fd = 0.0;
			unsigned long int nelems = bulk_mesh_pt->nboundary_element(4);

			Vector<double> s(2);
			s[1] = 0.0;
			double factor = 0.0;

			for (unsigned long int e = 0; e < nelems; e++) {
				int face_index = bulk_mesh_pt->face_index_at_boundary(4, e);
				EL *elem = dynamic_cast<EL *>(bulk_mesh_pt->boundary_element_pt(4, e));

				if (face_index == 1) {
					s[0] = 1.0;
					factor = -D[0];
				} else if (face_index == -1) {
					s[0] = -1.0;
					factor = D[0] * k[1] / k[0];
				}
				elem->get_flux(s, flux);
				get_flux_fd(0, elem, flux_fd, face_index);

				tot_flux += factor * flux[0];
				tot_flux_fd += factor * flux_fd[0];
			}
			
			printf("Setting total flux to %8.6f (fd: %8.6f) interface at %16.14f\n", tot_flux / ny, tot_flux_fd / ny, geometry->get_interface());
			geometry->set_flux(tot_flux / ny);
		}

		// void actions_after_newton_step() {
		void actions_before_newton_convergence_check() {
			for (unsigned s = 0; s < bulk_mesh_pt->nspine(); s++) {
				bulk_mesh_pt->spine_pt(s)->height() = geometry->x1();
			}

			bulk_mesh_pt->node_update();
		}

		void set_initial_condition() {
			time_pt()->time() = t_shift;
			unsigned long int nnode = bulk_mesh_pt->nnode();

			Vector<double> x(2);
			Vector<double> u(1);
			Vector<double> s(2);
			Vector<double> flux(2);
			Vector<double> r(2);

			unsigned int tsteps = time_stepper_pt()->nprev_values();
			double start_iface_x = sqrt(D[2] * time_pt()->time(tsteps));
			geometry->set_interface(start_iface_x);
			bulk_mesh_pt->node_update();

			for (unsigned long int n = 0; n < nnode; n++)
				time_stepper_pt()->assign_initial_positions_impulsive(bulk_mesh_pt->node_pt(n));

			for (unsigned long int n = 0; n < nnode; n++) {
				bulk_mesh_pt->node_pt(n)->position(tsteps, x);
				get_exact_u(time_pt()->time(tsteps), x, u);
				bulk_mesh_pt->node_pt(n)->set_value(tsteps, 0, u[0]);
			}

			unsigned int step = 0;
			doc_step(step, tsteps);
			step++;

			s[0] = 1.0;
			s[1] = 0.0;
			unsigned long int nelems = bulk_mesh_pt->nboundary_element(4);

			for (int t = tsteps -1; t >= 0; t--) {
				double time = time_pt()->time((unsigned int) t);

				//
				// update_interface
				//
				double new_h = 0.0;
				if (D[2] == 0.0) {
					double tot_flux = 0.0;
					for (unsigned long int e = 0; e < nelems; e++) {
						int face_index = bulk_mesh_pt->face_index_at_boundary(4, e);
						if (face_index == 1) {
							EL *elem = dynamic_cast<EL *>(bulk_mesh_pt->boundary_element_pt(4, e));
							get_flux_ic(t+1, elem, s, flux);
							tot_flux += k[0] * flux[0];
						} else if (face_index == -1) {
							EL *elem = dynamic_cast<EL *>(bulk_mesh_pt->boundary_element_pt(4, e));
							get_flux_ic(t+1, elem, s, flux);
							tot_flux -= k[1] * flux[0];
						}
					}
					double v = tot_flux / (rho[0] * L * ny);
					new_h = geometry->get_interface() + v*dt;
					// printf("IC %u: v=%8.6f x_old=%8.6f x_new=%8.6f\n", t, v, geometry->get_interface(), new_h);
				} else {
					new_h = sqrt(D[2] * time);
				}

				printf("IC %u: x_old=%8.6f x_new=%8.6f\n", t, geometry->get_interface(), new_h);
				geometry->set_interface(new_h);
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

			double tot_error = 0.0;

			char fname[256];
			sprintf(fname, "%s/steps/step%u.dat", info.directory().c_str(), info.number());
			FILE *file = fopen(fname, "w");
			for (unsigned long int n = 0; n < nnode; n++) {
				bulk_mesh_pt->node_pt(n)->position(t, x);
				bulk_mesh_pt->node_pt(n)->value(t, numerical_u);
				get_exact_u(time, x, exact_u);

				double error = numerical_u[0] - exact_u[1];
				tot_error += error * error;

				double redim_numerical_u = T_l - (T_l - T_s) * numerical_u[0];
				double redim_exact_u = T_l - (T_l - T_s)  * exact_u[0];

				fprintf(file, "%16.14f %16.14f %16.14f %16.14f %16.14f\n", x[0], x[1], redim_exact_u, redim_numerical_u, error);
			}
			fclose(file);

			tot_error = sqrt(tot_error) / nnode;

			double max_elem_size, min_elem_size;
			bulk_mesh_pt->max_and_min_element_size(max_elem_size, min_elem_size);

			// printf("Max element size: %16.14f min element size: %16.14f\n", max_elem_size, min_elem_size);
			sprintf(fname, "%s/results.dat", info.directory().c_str());
			file = fopen(fname, "a");
			fprintf(file, "%16.14f %16.14f %16.14f %16.14f %16.14f\n", time, tot_error, geometry->get_interface(), sqrt(D[2]*time), max_elem_size/ys[1]);
			fclose(file);

			printf("[%4u] time=%8.6f error=%e iface_err=%e interface=%8.6f expected=%8.6f\n", timestep, time, tot_error, fabs(geometry->get_interface() - sqrt(D[2] * time)), geometry->get_interface(), sqrt(D[2]*time));
			info.number()++;
		}

		// modified from UnsteadyHeatEquations::get_flux from unsteady_heat_elements.h
		void get_flux_ic(const unsigned int &t, EL * elem, const Vector<double> &s, Vector<double> &flux) const {
			// Find out how many nodes there are in the element
			unsigned n_node = elem->nnode();

			// Find the index at which the variable is stored
			unsigned u_nodal_index = elem->u_index_ust_heat();

			// Set up memory for the shape and test functions
			Shape psi(n_node);
			DShape dpsidx(n_node, 2);

			// Call the derivatives of the shape and test functions
			elem->dshape_eulerian(s, psi, dpsidx);

			// Initialise to zero
			for (unsigned j = 0; j < 2; j++) {
				flux[j] = 0.0;
			}

			// Loop over nodes
			for (unsigned l = 0; l < n_node; l++) {
				// Loop over derivative directions
				for (unsigned j = 0; j < 2; j++) {
					flux[j] += elem->nodal_value(t, l, u_nodal_index) * dpsidx(l, j);
				}
			}
		}

		void get_flux_fd(const unsigned &t, EL *elem, Vector<double> &flux, int direction) {
			unsigned n_node = elem->nnode_1d();
			Vector<Vector<double>> pos(n_node);
			Vector<double> fluxes(n_node);
			Vector<Node *> nodes(n_node);

			for (unsigned int i = 0; i < n_node; i++) pos[i].reserve(2);

			// n1 is to loop in y-direction
			// n2 is main loop in x-direciton
			for (unsigned int n1 = 0; n1 < n_node; n1++) {
				for (unsigned int n2 = 0; n2 < n_node; n2++) {
					nodes[n2] = elem->node_pt(n1 * n_node + n2);
					nodes[n2]->position(t, pos[n2]);
				}

				double dx = pos[n_node-1][0] - pos[0][0];
				double dT = 0.0;
				
				if (direction > 0) {
					// solid phase 
					switch(n_node) {
						case 2: dT = - nodes[0]->value(0, t) +     nodes[1]->value(0, t);															break;
						case 3: dT =   nodes[0]->value(0, t) - 4.0*nodes[1]->value(0, t) + 3.0*nodes[2]->value(0, t);								break;
						case 4: dT = - nodes[0]->value(0, t) + 4.5*nodes[1]->value(0, t) - 9.0*nodes[2]->value(0, t) + 5.5*nodes[3]->value(0, t);	break;
					}
				} else {
					// liquid phase
					switch (n_node) {
						case 2: dT =  nodes[1]->value(0, t) -     nodes[0]->value(0, t);															break;
						case 3: dT = -nodes[2]->value(0, t) + 4.0*nodes[1]->value(0, t) - 3.0*nodes[0]->value(0, t);								break;
						case 4: dT =  nodes[3]->value(0, t) - 4.5*nodes[2]->value(0, t) + 9.0*nodes[1]->value(0, t) - 5.5*nodes[0]->value(0, t);	break;
					}
				}
				fluxes[n1] = dT / dx;
			}
			flux[0] = fluxes[0];
		}
};

int main(int argc, char **argv) {
	CommandLineArgs::setup(argc, argv);

	uint Nx1 = 0;
	uint Nx2 = 0;
	uint Ny = 10;
	uint t_steps = 100;
	double dt = 0.0;
	double t_shift = 0.0;
	bool var_dt = false;
	uint write_freq = 1;
	std::string dname;

	CommandLineArgs::specify_command_line_flag("--nx", &Nx1);
	CommandLineArgs::specify_command_line_flag("--ny", &Ny);
	CommandLineArgs::specify_command_line_flag("--tsteps", &t_steps);
	CommandLineArgs::specify_command_line_flag("--dt", &dt);
	CommandLineArgs::specify_command_line_flag("--tshift", &t_shift);
	CommandLineArgs::specify_command_line_flag("--vardt");
	CommandLineArgs::specify_command_line_flag("--write-freq", &write_freq);
	CommandLineArgs::specify_command_line_flag("--k1", &k[0]);
	CommandLineArgs::specify_command_line_flag("--k2", &k[1]);
	CommandLineArgs::specify_command_line_flag("--rho1", &rho[0]);
	CommandLineArgs::specify_command_line_flag("--rho2", &rho[1]);
	CommandLineArgs::specify_command_line_flag("--cp1", &Cp[0]);
	CommandLineArgs::specify_command_line_flag("--cp2", &Cp[1]);
	CommandLineArgs::specify_command_line_flag("--L", &L);
	CommandLineArgs::specify_command_line_flag("--De", &D[2]);
	CommandLineArgs::specify_command_line_flag("--dname", &dname, "doc");

	CommandLineArgs::parse_and_assign();

	var_dt = CommandLineArgs::command_line_flag_has_been_set("--vardt");

	if (Nx1 == 0) {
		cout << "Error: Nx not specified" << endl;
		exit(1);
	}

	if (dt == 0.0 && !var_dt) {
		cout << "Error: dt not specified" << endl;
		exit(1);
	}

	if (var_dt) {
		if (dt != 0.0) {
			cout << "Error: dt has been specified with variable dt flag, check cmdline arguments" << endl;
			exit(1);
		}

		dt = pow(((double)1/(double)Nx1), 2);
	}

	if (!CommandLineArgs::command_line_flag_has_been_set("--dname")) {
		char temp[256];
		sprintf(temp, "RESLT/%dn%u_%dt%.2e", X_ORDER, Nx1, T_ORDER, dt);
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

	// set approximate starting point for interface based on previous runs
	double _interface_time = t_shift - T_ORDER*dt;
	xs[1] = (D[2] == 0.0) ? sqrt(1.26 * _interface_time) : sqrt(D[2] * _interface_time);
	
	for (int i = 0; i < 2; i++) D[i] = k[i]/(Cp[i]*rho[i]);

	St = L / (Cp[0] * (T_l - T_s));

	printf("Problem Def:\n");
	printf("\tk1=%8.6f k2=%8.6f\n", k[0], k[1]);
	printf("\trho1=%8.6f rho2=%8.6f\n", rho[0], rho[1]);
	printf("\tCp1=%8.6f Cp2=%8.6f\n", Cp[0], Cp[1]);
	printf("\tD1=%8.6f D2=%8.6f De=%8.6f\n", D[0], D[1], D[2]);
	printf("\tx0=%8.6f x1=%8.6f x2=%8.6f\n", xs[0], xs[1], xs[2]);
	printf("\tL=%8.6f St=%8.6f\n", L, St);

	Nx2 = Nx1;
	auto problem = Erf2DProblem<QUnsteadyHeatElement<2,X_ORDER>>(Nx1, Nx2, Ny, t_steps, dt, t_shift, info);

	problem.initialise_dt(dt);
	problem.set_initial_condition();

	// save current parameters to file for future reference if needed
	char config_fname[256];
	sprintf(config_fname, "%s/config", dname.c_str());

	FILE *file = fopen(config_fname, "w");
	fprintf(file, "nx1=%d\n", Nx1);
	fprintf(file, "nx2=%d\n", Nx2);
	fprintf(file, "nx=%d\n", Nx1+Nx2);
	fprintf(file, "x_order=%d\n", X_ORDER);
	fprintf(file, "dt=%e\n", dt);
	fprintf(file, "t_shift=%e\n", t_shift);
	fprintf(file, "t_order=%d\n", T_ORDER);
	fprintf(file, "t_steps=%d\n", t_steps);
	fprintf(file, "write_freq=%u\n", write_freq);
	fprintf(file, "k1=%10.8f\n", k[0]);
	fprintf(file, "k2=%10.8f\n", k[1]);
	fprintf(file, "rho1=%10.8f\n", rho[0]);
	fprintf(file, "rho2=%10.8f\n", rho[1]);
	fprintf(file, "cp1=%10.8f\n", Cp[0]);
	fprintf(file, "cp2=%10.8f\n", Cp[1]);
	fprintf(file, "L=%10.8f\n", L);
	fprintf(file, "De=%10.8f\n", D[2]);
	fclose(file);

	int prev_steps = problem.time_stepper_pt()->nprev_values()+1;

	for (uint t = 0; t < t_steps; t++) {
		problem.unsteady_newton_solve(dt);

		if (t % write_freq == 0 || t == t_steps-1)
			problem.doc_step(t+prev_steps);

		double x_int = geometry->get_interface();

		if (x_int >= xs[2] || isnan(x_int)) {
			printf("Interface reached right boundary, exiting...\n");
			break;
		}
	}
}