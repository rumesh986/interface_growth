
// moving interface
// calculate as part of solution in newton iterations (hopefully)

#include <filesystem>

#include "includes.h"
#include "two_phase_free_boundary_mesh.h"

#ifndef RUN_SCRIPT
#define X_ORDER 2
#define T_ORDER 1
#endif

double k[2] = {1.0, 0.5};
double rho[2] = {1.0, 0.5};
double Cp[2] = {1.0, 0.5};

double D[3] {k[0]/(Cp[0]*rho[0]), k[1]/(Cp[1]*rho[1]), 0.0};
double L = 1.0;

bool ic_set = false;

static const double T_s = -1.0;
static const double T_m = 0.0;
static const double T_l = 1.0;

double xs[3] = {-1.0, 0.025, 1.0};
static const double ys[2] = {0.0, 1.0};

FreeBoundaryElement *geometry;

void no_flux_fct(const double &t, const Vector<double> &x, double &flux) {
	flux = 0.0;
}

void get_source(const double &t, const Vector<double> &x, double &source) {
	source = 0.0;
}

// analytical solution
void get_exact_u(const double &t, const Vector<double> &x, Vector<double> &u) {
	// double h = geometry->get_interface();
	// if (x[0] < h) {
	// 	u[0] = T_s + (T_m-T_s)/(1+erf(h/(2*sqrt(D[0]*t)))) * (1 + erf(x[0]/(2*sqrt(D[0]*t))));
	// } else {
	// 	u[0] = T_l - (T_l-T_m)/(1-erf(h/(2*sqrt(D[1]*t)))) * (1 - erf(x[0]/(2*sqrt(D[1]*t))));
	// }

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
		double erf_iface = erf(h/(2*sqrt(D[0]*t)));
		u[0] = ((T_m - T_s)*erf(x[0]/(2*sqrt(D[0]*t))) + T_m + T_s*erf_iface)/(1+erf_iface);

		erf_iface = erf(h_ana/(2*sqrt(D[0]*t)));
		u[1] = ((T_m - T_s)*erf(x[0]/(2*sqrt(D[0]*t))) + T_m + T_s*erf_iface)/(1+erf_iface);
	} else {
		double erf_iface = erf(h/(2*sqrt(D[1]*t)));
		u[0] = ((T_l - T_m)*erf(x[0]/(2*sqrt(D[1]*t))) + T_m - T_l*erf_iface)/(1-erf_iface);

		erf_iface = erf(h_ana/(2*sqrt(D[1]*t)));
		u[1] = ((T_l - T_m)*erf(x[0]/(2*sqrt(D[1]*t))) + T_m - T_l*erf_iface)/(1-erf_iface);
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
	{0, no_flux_fct},
	{2, no_flux_fct},
};


template<class EL>
class Erf2DProblem : public Problem {
	private:
		unsigned int nx1, nx2, nx, ny;
		unsigned int t_steps;
		double dt, t_shift;

		DocInfo info;
	
		TwoPhaseFreeBoundaryMesh<MacroElementNodeUpdateElement<EL>> *bulk_mesh_pt;
		Mesh *geometry_mesh_pt;

		FreeBoundaryDomain *domain;
	public:
		Erf2DProblem(
			unsigned int nx_,
			unsigned int ny_,
			unsigned int tsteps,
			double dt_,
			double tshift,
			DocInfo info_
		) : nx1(nx_), nx2(nx_), nx(nx_+nx_), ny(ny_), t_steps(tsteps), dt(dt_), t_shift(tshift), info(info_) {

			add_time_stepper_pt(new BDF<T_ORDER>);

			// TimeStepper *domain_ts_pt = new BDF<T_ORDER>;

			// FreeBoundaryGeometry *geometry_object = new FreeBoundaryGeometry(xs[0], xs[1], xs[2], ys[0], ys[1], time_stepper_pt());

			geometry = new FreeBoundaryElement(xs[0], xs[1], xs[2], ys[0], ys[1], rho[0] * L, time_stepper_pt());
			domain = new FreeBoundaryDomain(geometry, nx1, nx2, ny);

			printf("x0=%8.6f x1=%8.6f x2=%8.6f\n", geometry->x0(), geometry->x1(), geometry->x2());

			bulk_mesh_pt = new TwoPhaseFreeBoundaryMesh<MacroElementNodeUpdateElement<EL>>(nx1, nx2, ny, geometry, domain, time_stepper_pt());
			// bulk_mesh_pt = new TwoPhaseFreeBoundaryMesh<MacroElementNodeUpdateElement<EL>>(nx1, nx2, ny, geometry_object, time_stepper_pt());
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

			assign_eqn_numbers();

			linear_solver_pt()->disable_doc_time();
			disable_info_in_newton_solve();
			max_newton_iterations() = 1e5;
			max_residuals() = 1e7;
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

		void actions_before_newton_solve() {};
		void actions_after_newton_solve() {};

		void actions_before_implicit_timestep() {};
		void actions_after_implicit_timestep() {};

		void actions_before_newton_convergence_check() {
			// Vector<double> s(2);
			Vector<double> flux(2);

			double tot_flux = 0.0;
			
			unsigned long int nelems = bulk_mesh_pt->nboundary_element(4);
			for (unsigned long int e = 0; e < nelems; e++) {
				int face_index = bulk_mesh_pt->face_index_at_boundary(4, e);
				if (face_index == 1) {
					Vector<double> s = {1.0, 0.0};

					EL *elem = dynamic_cast<EL *>(bulk_mesh_pt->boundary_element_pt(4, e));
					elem->get_flux(s, flux);

					tot_flux += k[0] * flux[0];
				} else if (face_index == -1) {
					Vector<double> s = {-1.0, 0.0};

					EL *elem = dynamic_cast<EL *>(bulk_mesh_pt->boundary_element_pt(4, e));
					elem->get_flux(s, flux);

					tot_flux -= k[1] * flux[0];
				}
			}

			printf("Setting total flux to %8.6f\n", tot_flux / ny);
			geometry->set_flux(tot_flux / ny);
			bulk_mesh_pt->node_update();
			printf("Moved interface from %8.6f to %8.6f\n", old_int, domain->get_interface());
		}

		void set_initial_condition() {
			time_pt()->time() = t_shift;
			unsigned long int nnode = bulk_mesh_pt->nnode();

			Vector<double> x(2);
			Vector<double> u(1);
			Vector<double> s(2);
			Vector<double> flux(2);
			Vector<double> r(2);

			for (unsigned long int n = 0; n < nnode; n++)
				time_stepper_pt()->assign_initial_positions_impulsive(bulk_mesh_pt->node_pt(n));
			
			unsigned int tsteps = time_stepper_pt()->nprev_values();

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
							tot_flux += k[1] * flux[0];
						}
					}
					double v = tot_flux / (rho[0] * L * ny);
					new_h = geometry->get_interface() + v*dt;
					// printf("IC %u: v=%8.6f x_old=%8.6f x_new=%8.6f\n", t, v, geometry->get_interface(), new_h);
				} else {
					new_h = sqrt(D[2] * time);
				}

				geometry->set_interface(new_h);
				printf("IC %u: x_old=%8.6f x_new=%8.6f\n", t, geometry->get_interface(), new_h);
				bulk_mesh_pt->node_update();

				// for (unsigned int n = 0; n < nnode; n++) {
				// 	MacroElementNodeUpdateNode *node_pt = dynamic_cast<MacroElementNodeUpdateNode *>(bulk_mesh_pt->node_pt(n));
				// 	node_update_ic(t, node_pt);
				// }

				// mesh_pt()->node_update() only moves the nodes for current timestep
				// when setting initial condition, we need to move the nodes ourselves
				// the below code has been adapted from mesh_pt()->node_update()
				// std::map<Node *, bool> node_handled;

				// nelems = nx * ny;
				
				// // only work with each node once
				// for (unsigned long int n = 0; n < nnode; n++) {
				// 	Node *node_pt = bulk_mesh_pt->node_pt(n);
				// 	node_handled[node_pt] = false;
				// }

				// for (unsigned long int e = 0; e < nelems; e++) {
				// 	FiniteElement *elem_pt = dynamic_cast<FiniteElement *>(bulk_mesh_pt->element_pt(e));
				// 	unsigned long int elem_nnode = elem_pt->nnode();
				// 	for (unsigned long int n = 0; n < elem_nnode; n++) {
				// 		Node *node_pt = elem_pt->node_pt(n);

				// 		if (!node_handled[node_pt]) {
				// 			elem_pt->local_coordinate_of_node(n, s);
				// 			elem_pt->get_x(t, s, r);
				// 			if (elem_pt->macro_elem_pt() == 0) {
				// 				printf("We have big problems here !!!!\n\n");
				// 			}

				// 			for (int i = 0; i < 2; i++)
				// 				node_pt->x(t, i) = r[i];
							
				// 			node_handled[node_pt] = true;
				// 		}
				// 	}
				// }


				for (unsigned long int n = 0; n < nnode; n++) {
					bulk_mesh_pt->node_pt(n)->position(t, x);
					get_exact_u(time, x, u);
					bulk_mesh_pt->node_pt(n)->set_value(t, 0, u[0]);
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

				fprintf(file, "%16.14f %16.14f %16.14f %16.14f %16.14f\n", x[0], x[1], exact_u[1], numerical_u[0], error);
			}
			fclose(file);

			tot_error = sqrt(tot_error) / nnode;

			double max_elem_size, min_elem_size;
			bulk_mesh_pt->max_and_min_element_size(max_elem_size, min_elem_size);

			sprintf(fname, "%s/results.dat", info.directory().c_str());
			file = fopen(fname, "a");
			fprintf(file, "%16.14f %16.14f %16.14f %16.14f %16.4f\n", time, tot_error, geometry->get_interface(), sqrt(D[2]*time), max_elem_size);
			fclose(file);

			printf("[%4u] time=%8.6f error = %e interface = %8.6f expected = %8.6f\n", timestep, time, tot_error, geometry->get_interface(), sqrt(D[2]*time));
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

		void node_update_ic(const unsigned int &t, MacroElementNodeUpdateNode *node_pt) {
			static unsigned long int error_count = 0;
			if (node_pt->node_update_element_pt() == 0)  {
				printf("We might still have problems!! %lu / %lu\n", error_count, t*bulk_mesh_pt->nnode());
				error_count++;
			} else {
				Vector<double> x_new(2), s_local(2);
				s_local = node_pt->s_in_node_update_element();
				node_pt->node_update_element_pt()->get_x(t, s_local, x_new);
				for (unsigned int i = 0; i < 2; i++) {
					node_pt->x(t, i) = x_new[i];
				}
				
			}
		}
};

int main(int argc, char **argv) {
	CommandLineArgs::setup(argc, argv);

	uint Nx = 0;
	uint Ny = 10;
	uint t_steps = 100;
	double dt = 0.0;
	double t_shift = 0.0;
	bool var_dt = false;
	uint write_freq = 1;
	std::string dname;

	CommandLineArgs::specify_command_line_flag("--nx", &Nx);
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

	if (Nx == 0) {
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

		dt = pow(((double)1/(double)Nx), 2);
	}

	if (!CommandLineArgs::command_line_flag_has_been_set("--dname")) {
		char temp[256];
		sprintf(temp, "RESLT/%dn%u_%dt%.2e", X_ORDER, Nx, T_ORDER, dt);
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

	printf("Problem Def:\n");
	printf("\tk1=%8.6f k2=%8.6f\n", k[0], k[1]);
	printf("\trho1=%8.6f rho2=%8.6f\n", rho[0], rho[1]);
	printf("\tCp1=%8.6f Cp2=%8.6f\n", Cp[0], Cp[1]);
	printf("\tD1=%8.6f D2=%8.6f De=%8.6f\n", D[0], D[1], D[2]);
	printf("\tL=%8.6f\n", L);

	auto problem = Erf2DProblem<QUnsteadyHeatElement<2,X_ORDER>>(Nx, Ny, t_steps, dt, t_shift, info);

	problem.initialise_dt(dt);
	problem.set_initial_condition();

	// save current parameters to file for future reference if needed
	char config_fname[256];
	sprintf(config_fname, "%s/config", dname.c_str());

	FILE *file = fopen(config_fname, "w");
	fprintf(file, "nx1=%d\n", Nx);
	fprintf(file, "nx2=%d\n", Nx);
	fprintf(file, "nx=%d\n", Nx+Nx);
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
	fclose(file);

	int prev_steps = problem.time_stepper_pt()->nprev_values()+1;

	for (uint t = 0; t < t_steps; t++) {
		problem.unsteady_newton_solve(dt);

		if (t % write_freq == 0)
			problem.doc_step(t+prev_steps);

		double x_int = geometry->get_interface();

		if (x_int >= xs[2] || isnan(x_int)) {
			printf("Interface reached right boundary, exiting...\n");
			break;
		}
	}
}