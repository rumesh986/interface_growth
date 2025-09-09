
// moving interface 
// interface moves by calculating instantaneous mesh velocity from flux
// right boundary set to Tfr and allowed to evolve
// use dh/dt estimates for setting history values in initial condition

// domain and its boundaries are defined as follows
// 							2
// 		-------------------------------------
// 		|					|				|
// 		|					|				|
// 		|					|				|
// 		|					|				|
// 	  3 |				  4 |				| 1
// 		|					|				|
// 		|					|				|
// 		|					|				|
// 		|					|				|
// 		-------------------------------------
// 							0


#include <filesystem>

#include "includes.h"
#include "two_layer_mesh.h"

#ifndef RUN_SCRIPT
#define X_ORDER 2
#define T_ORDER 1
#endif

double k1 = 1.0;
double k2 = 0.5;

double rho1 = 1.0;
double rho2 = 0.5;

double Cp1 = 1.0;
double Cp2 = 0.5;

double D1 = k1 / (Cp1 * rho1);
double D2 = k2 / (Cp2 * rho2);

double L = 1.0;

static const double Ts = -1.0;
static const double Tm = 0.0;
static const double Tfr = 1.0;

static const double x0 = -1.0;
double x1 = 0.025; // h0, initial interface position at start time
static const double x2 = 1.0;

static const double y_0 = 0.0;
static const double y_1 = 1.0;

static const double vel = 0.0;

TwoPhaseDomain *domain;

// list of boundaries that need to match analytical values
Vector<unsigned int> analytical_boundaries = {
	1, // right side
	3, // left side
};

// list of boundaries that are pinned
Vector<unsigned int> pinned_boundaries = {
	4 // interface
};

// function to implement no flux neumann conditions
void no_flux_fct(const double &t, const Vector<double> &x, double &flux) {
	flux = 0.0;
}

// mapping of boundaries to neumann conditions
unordered_map<unsigned int, FluxFctPt> neumann_boundaries = {
	{0, no_flux_fct}, // bottom side
	{2, no_flux_fct}  // top side
};

// analytical solution
void get_exact_u(const double &t, const Vector<double> &x, Vector<double> &u) {
	double h = domain->get_interface();
	if (x[0] < h) {
		u[0] = Ts + (Tm-Ts)/(1+erf(h/(2*sqrt(D1*t)))) * (1 + erf(x[0]/(2*sqrt(D1*t))));
	} else {
		u[0] = Tfr - (Tfr-Tm)/(1-erf(h/(2*sqrt(D2*t)))) * (1 - erf(x[0]/(2*sqrt(D2*t))));
	}
}

template<class EL>
class Erf2DProblem : public Problem {
	private:
		uint Nx1, Nx2, Nx, Ny;
		uint t_steps;
		double dt, t_shift;
		
		DocInfo info;

	public:
		Erf2DProblem(uint nx, uint ny,  uint t_steps, double dt, double t_shift, DocInfo info)
		: Nx1(nx), Nx2(nx), Nx(nx+nx), Ny(ny), t_steps(t_steps), dt(dt), t_shift(t_shift), info(info) {
			// calls super constructor as Problem() has no arguments

			// prepare problem options
			add_time_stepper_pt(new BDF<T_ORDER>);
			problem_is_nonlinear(false);

			// prepare mesh
			domain = new TwoPhaseDomain(x0, x1, x2, vel, Nx1, Nx2, Ny, time_pt());
			mesh_pt() = new RefineableTwoLayer2DMesh<EL>(Nx1, Nx2, Ny, x0, x1, x2, y_0, y_1, domain, time_stepper_pt());
			mesh_pt()->setup_boundary_element_info();

			for (uint e = 0; e < mesh_pt()->nelement(); e++)
				dynamic_cast<EL *>(mesh_pt()->element_pt(e))->source_fct_pt() = get_source;
			
			// pin nodes of pinned boundaries
			for (unsigned int b : pinned_boundaries) {
				unsigned long int nnode = mesh_pt()->nboundary_node(b);
				for (unsigned long int n = 0; n < nnode; n++)
					mesh_pt()->boundary_node_pt(b, n)->pin_all();
			}

			// pin analytical boundaries
			for (unsigned int b : analytical_boundaries) {
				unsigned long int nnode = mesh_pt()->nboundary_node(b);
				for (unsigned int n = 0; n < nnode; n++)
					mesh_pt()->boundary_node_pt(b, n)->pin_all();
			}

			// set neumann boundaries
			for (auto iter : neumann_boundaries) {
				create_flux_elements(iter.first, iter.second);
			}

			// Set physical parameters for every element
			for (uint yi = 0; yi < Ny; yi++) {
				for (uint e = 0; e < Nx1; e++)
					dynamic_cast<EL *>(mesh_pt()->element_pt(yi*Nx + e))->beta_pt() = &D1;
				
				for (uint e = Nx1; e < Nx; e++)
					dynamic_cast<EL *>(mesh_pt()->element_pt(yi*Nx + e))->beta_pt() = &D2;
			}

			assign_eqn_numbers();

			// reduce verbosity in logs
			linear_solver_pt()->disable_doc_time();
			disable_info_in_newton_solve();
		}

		~Erf2DProblem() {
			delete mesh_pt();
		}

		static void get_source(const double &t, const Vector<double> &x, double &source) {
			source = 0.0;
		}

		void actions_before_newton_solve() {};
		void actions_after_newton_solve() {};

		void actions_before_implicit_timestep() {
			update_interface(); // move interface

			// reset analytical boundaries
			Vector<double> x(2);
			Vector<double> u(1);

			double cur_t = time_pt()->time() + dt;

			for (unsigned long int b : analytical_boundaries) {
				unsigned long int nnode = mesh_pt()->nboundary_node(b);
				for (unsigned int n = 0; n < nnode; n++) {
					mesh_pt()->boundary_node_pt(b, n)->position(x);
					get_exact_u(cur_t, x, u);
					mesh_pt()->boundary_node_pt(b, n)->set_value(0, u[0]);
				}
			}
		}

		void actions_after_implicit_timestep() {}

		void set_initial_conditions() {
			time_pt()->time() = t_shift;
			unsigned long int nnode = mesh_pt()->nnode();

			Vector<double> x(2);
			Vector<double> u(1);

			// set positions of nodes in past history values
			for (unsigned long int n = 0; n < nnode; n++)
				time_stepper_pt()->assign_initial_positions_impulsive(mesh_pt()->node_pt(n));

			printf("Setting Initial conditions\n");

			unsigned int tsteps = time_stepper_pt()->nprev_values();

			// set initial values
			for (unsigned long int n = 0; n < nnode; n++) {
				mesh_pt()->node_pt(n)->position(tsteps, x);
				get_exact_u(time_pt()->time(tsteps), x, u);
				mesh_pt()->node_pt(n)->set_value(tsteps, 0, u[0]);
			}

			int step = 0;
			// save current state to file
			doc_step(step, tsteps);
			step++;

			// set past values for initial conditions and interface positions
			for (int t = tsteps-1; t >= 0; t--) {
				double time = time_pt()->time((uint) t);

				update_interface(t, true);

				// set initial values
				for (unsigned long int n = 0; n < nnode; n++) {
					mesh_pt()->node_pt(n)->position(t, x);
					get_exact_u(time, x, u);
					mesh_pt()->node_pt(n)->set_value(t, 0, u[0]);
				}

				printf("[% 4d] Setting initial condition at t=%8.6f\n", step, time);
				doc_step(step, t);
				step++;
			}

			time_pt()->time() = t_shift;
		}

		// method that calculates next interface position and actually moves the interface between timesteps
		void update_interface(const unsigned int &t = 0, bool ic = false) {
			Vector<double> s(2);
			Vector<double> flux(2);

			s[0] = 1.0;
			s[1] = 0.0;

			double tot_flux = 0.0;

			unsigned long int nelems = mesh_pt()->nboundary_element(4);
			for (unsigned long int e = 0; e < nelems; e++) {
				// only choose fluxes in the horizontal direction
				int face_index = mesh_pt()->face_index_at_boundary(4, e);
				if (face_index == 1) {
					EL *elem = dynamic_cast<EL *>(mesh_pt()->boundary_element_pt(4, e));
					if (ic) get_flux_ic(t+1, elem, s, flux);
					else 	elem->get_flux(s, flux);

					tot_flux += k1 * flux[0];
				} else if (face_index == -1) {
					EL *elem = dynamic_cast<EL *>(mesh_pt()->boundary_element_pt(4, e));
					if (ic) get_flux_ic(t+1, elem, s, flux);
					else 	elem->get_flux(s, flux);

					tot_flux += k2 * flux[0];
				}
			}

			// calculate dh/dt and new interface position
			double v = tot_flux / (rho1 * L * Ny);
			double new_h = domain->get_interface() + v * dt;

			domain->set_interface(new_h);
			mesh_pt()->node_update();

			// mesh_pt()->node_update() only moves the nodes for current timestep
			// when setting initial condition, we need to move the nodes ourselves
			// the below code has been adapted from mesh_pt()->node_update()
			if (ic) {
				std::map<Node *, bool> node_handled;
				Vector<double> r(2);

				unsigned long int nelems = Nx*Ny;
				unsigned long int nnode_total = mesh_pt()->nnode();
				
				// only work with each node once
				for (unsigned long int n = 0; n < nnode_total; n++) {
					Node *node_pt = mesh_pt()->node_pt(n);
					node_handled[node_pt] = false;
				}

				for (unsigned long int e = 0; e < nelems; e++) {
					FiniteElement *elem_pt = dynamic_cast<FiniteElement *>(mesh_pt()->element_pt(e));
					unsigned long int nnode = elem_pt->nnode();
					for (unsigned long int n = 0; n < nnode; n++) {
						Node *node_pt = elem_pt->node_pt(n);

						if (!node_handled[node_pt]) {
							elem_pt->local_coordinate_of_node(n, s);
							elem_pt->get_x(t, s, r);
							if (elem_pt->macro_elem_pt() == 0) {
								printf("We have big problems here !!!!\n\n");
							}

							for (int i = 0; i < 2; i++)
								node_pt->x(t, i) = r[i];
							
							node_handled[node_pt] = true;
						}
					}
				}
			}

			printf("moved interface to x=%8.6f (v=%8.6f, flux=%8.6f)\n", new_h, v, tot_flux);
		}

		// Neumann boundaries need special flux elements associated at the boundaries
		// this method will create the required elements and attach the correct flux functions to them
		void create_flux_elements(unsigned int b, FluxFctPt &flux_pt) {
			unsigned int nelems = mesh_pt()->nboundary_element(b);
			for (unsigned int e = 0; e < nelems; e++) {
				EL *elem = dynamic_cast<EL *>(mesh_pt()->boundary_element_pt(b, e));
				int face_i = mesh_pt()->face_index_at_boundary(b, e);
				UnsteadyHeatFluxElement<EL> *flux_elem = new UnsteadyHeatFluxElement<EL>(elem, face_i);
				flux_elem->flux_fct_pt() = flux_pt;
				mesh_pt()->add_element_pt(flux_elem);
			}
		}

		// save simulation data to file
		void doc_step(const unsigned int &timestep, const unsigned int &t = 0) {
			double time = time_pt()->time(t);
			unsigned long int nnode = mesh_pt()->nnode();

			Vector<double> x(2);
			Vector<double> exact_u(1);
			Vector<double> numerical_u(1);

			double tot_error = 0.0;

			// save current step information
			char fname[512];
			sprintf(fname, "%s/steps/step%u.dat", info.directory().c_str(), info.number());
			FILE *file = fopen(fname, "w");
			for (unsigned long int n = 0; n < nnode; n++) {
				mesh_pt()->node_pt(n)->position(t, x);
				mesh_pt()->node_pt(n)->value(t, numerical_u);
				get_exact_u(time, x, exact_u);

				double error = numerical_u[0] - exact_u[0];
				tot_error += error * error;

				fprintf(file, "%16.14f %16.14f %16.14f %16.14f %16.14f\n", x[0], x[1], exact_u[0], numerical_u[0], error);
			}
			fclose(file);
			
			tot_error = sqrt(tot_error) / nnode;

			// save errors and interface positions in seperate file (has information from all timesteps)
			sprintf(fname, "%s/results.dat", info.directory().c_str());
			file = fopen(fname, "a");
			fprintf(file, "%16.14f %16.14f %16.14f\n", time, tot_error, domain->get_interface());
			fclose(file);

			printf("[%4u] time=%8.6f error = %e interface = %8.6f\n", timestep, time, tot_error, domain->get_interface());
			info.number()++;
		}

		// modified from UnsteadyHeatEquations::get_flux from unsteady_heat_elements.h
		void get_flux_ic(const unsigned int &t, EL * elem, const Vector<double>& s, Vector<double>& flux) const {
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
	CommandLineArgs::specify_command_line_flag("--k1", &k1);
	CommandLineArgs::specify_command_line_flag("--k2", &k2);
	CommandLineArgs::specify_command_line_flag("--rho1", &rho1);
	CommandLineArgs::specify_command_line_flag("--rho2", &rho2);
	CommandLineArgs::specify_command_line_flag("--cp1", &Cp1);
	CommandLineArgs::specify_command_line_flag("--cp2", &Cp2);
	CommandLineArgs::specify_command_line_flag("--L", &L);
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
	x1 = sqrt(1.26 * t_shift);

	printf("Problem Def:\n");
	printf("\tk1=%8.6f k2=%8.6f\n", k1, k2);
	printf("\trho1=%8.6f rho2=%8.6f\n", rho1, rho2);
	printf("\tCp1=%8.6f Cp2=%8.6f\n", Cp1, Cp2);
	printf("\tD1=%8.6f D2=%8.6f\n", D1, D2);

	auto problem = Erf2DProblem<RefineableQUnsteadyHeatElement<2,X_ORDER>>(Nx, Ny, t_steps, dt, t_shift, info);

	problem.initialise_dt(dt);
	problem.set_initial_conditions();

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
	fprintf(file, "k1=%10.8f\n", k1);
	fprintf(file, "k2=%10.8f\n", k2);
	fprintf(file, "rho1=%10.8f\n", rho1);
	fprintf(file, "rho2=%10.8f\n", rho2);
	fprintf(file, "cp1=%10.8f\n", Cp1);
	fprintf(file, "cp2=%10.8f\n", Cp2);
	fprintf(file, "L=%10.8f\n", L);
	fclose(file);

	int prev_steps = problem.time_stepper_pt()->nprev_values()+1;

	for (uint t = 0; t < t_steps; t++) {
		problem.unsteady_newton_solve(dt);

		if (t % write_freq == 0)
			problem.doc_step(t+prev_steps);

		double x_int = domain->get_interface();

		if (x_int >= x2 || isnan(x_int)) {
			printf("Interface reached right boundary, exiting...\n");
			break;
		}
	}
}