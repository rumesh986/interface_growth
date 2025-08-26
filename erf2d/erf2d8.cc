
// moving interface 
// interface moves by calculating instantaneous mesh velocity from flux
// right boundary set to Tfr and allowed to evolve
// use dh/dt estimates for setting history values in initial condition

#include <filesystem>

#include "includes.h"
#include "two_layer_mesh.h"

#ifndef RUN_SCRIPT
#define X_ORDER 2
#define T_ORDER 1
#endif

static const double D1 = 1.0;
static const double D2 = 0.5;

static const double Ts = -1.0;
static const double Tm = 0.0;
static const double Tfr = 1.0;

static const double alpha = 1.0;

static const double x0 = -1.0;
static const double x1 = 0.025; // h0
static const double x2 = 1.0;

static const double y_0 = 0.0;
static const double y_1 = 1.0;

static const double vel = 0.0;

TwoPhaseDomain *domain;

Vector<unsigned int> analytical_boundaries = {
	1, // right side
	3, // left side
	4  // interface
};

Vector<unsigned int> pinned_boundaries = {
	// 4 // interface
};

void no_flux_fct(const double &t, const Vector<double> &x, double &flux) {
	flux = 0.0;
}

unordered_map<unsigned int, FluxFctPt> neumann_boundaries = {
	{0, no_flux_fct},
	{2, no_flux_fct}
};

void get_exact_u(const double &t, const Vector<double> &x, Vector<double> &u) {
	double h = domain->get_interface();
	if (x[0] < h) {
		u[0] = Ts + (Tm-Ts)/(1+erf(h/(2*sqrt(D1*t)))) * (1 + erf(x[0]/(2*sqrt(D1*t))));
	} else {
		u[0] = Tfr - (Tfr-Tm)/(1-erf(h/(2*sqrt(D2*t)))) * (1 - erf(x[0]/(2*sqrt(D2*t))));
	}
}

void get_initial_u(const Vector<double> &x, Vector<double> &u) {
	// u[0] = (x[0] < x1) ? Ts : Tfr;
	// if (x[0] < x1) u[0] = Ts;
	// else if (x[0] == x1) u[0] = Tm;
	// else u[0] = Tfr;
	
	const unsigned int nx_base = 10;
	double dx = (x1 - x0) / nx_base;

	int xi = x[0] / dx;
	// printf("x0 =% 8.6f xi =%3d\n", x[0], xi);
	if (xi < 0) {
		u[0] = Ts;
	} else if (xi == 0) {
		u[0] = x[0] * (Tfr - Ts)/(2*dx);
	} else {
		u[0] = Tfr;
	}


}

void dhdt(const double &t, const Vector<double> &x, Vector<double> &v) {
	double h = domain->get_interface();
	v[0] = 2*(Tm-Ts)/(alpha * sqrt(Pi*D1*t) * (1 + erf(h/(2*sqrt(D1*t))))) * exp(-1 * h*h/(4*D1*t));
}

template<class EL>
class Erf2D6Problem : public Problem {
	private:
		uint Nx1, Nx2, Nx, Ny;
		uint t_steps;
		double dt, t_shift;
		
		DocInfo info;

	public:
		Erf2D6Problem(uint nx, uint ny,  uint t_steps, double dt, double t_shift, DocInfo info)
		: Nx1(nx), Nx2(nx), Nx(nx+nx), Ny(ny), t_steps(t_steps), dt(dt), t_shift(t_shift), info(info) {
			add_time_stepper_pt(new BDF<T_ORDER>);

			domain = new TwoPhaseDomain(x0, x1, x2, vel, Nx1, Nx2, Ny, time_pt());
			mesh_pt() = new RefineableTwoLayer2DMesh<EL>(Nx1, Nx2, Ny, x0, x1, x2, y_0, y_1, domain, time_stepper_pt());
			mesh_pt()->setup_boundary_element_info();

			printf("nnode 1d: %u\n", dynamic_cast<EL *>(mesh_pt()->element_pt(0))->nnode_1d());

			for (uint e = 0; e < mesh_pt()->nelement(); e++)
				dynamic_cast<EL *>(mesh_pt()->element_pt(e))->source_fct_pt() = get_source;
				
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

			for (uint yi = 0; yi < Ny; yi++) {
				for (uint e = 0; e < Nx1; e++)
					dynamic_cast<EL *>(mesh_pt()->element_pt(yi*Nx + e))->beta_pt() = (double *) &D1;
				
				for (uint e = Nx1; e < Nx; e++)
					dynamic_cast<EL *>(mesh_pt()->element_pt(yi*Nx + e))->beta_pt() = (double *) &D2;
			}

			assign_eqn_numbers();

			linear_solver_pt()->disable_doc_time();
			disable_info_in_newton_solve();
		}

		~Erf2D6Problem() {
			delete mesh_pt();
		}

		static void get_source(const double &t, const Vector<double> &x, double &source) {
			source = 0.0;
		}

		void actions_before_newton_solve() {};
		void actions_after_newton_solve() {};

		void actions_before_implicit_timestep() {
			update_interface();

			// reset analytical boundaries
			Vector<double> x(2);
			Vector<double> u(1);

			double cur_t = time_pt()->time();

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
			Vector<double> v(1);

			Vector<double> s(2);
			Vector<double> flux(2);

			// set positions of nodes in past history values
			for (unsigned long int n = 0; n < nnode; n++)
				time_stepper_pt()->assign_initial_positions_impulsive(mesh_pt()->node_pt(n));

			printf("Setting Initial conditions\n");

			unsigned int tsteps = time_stepper_pt()->nprev_values();

			// set initial values
			for (unsigned long int n = 0; n < nnode; n++) {
				mesh_pt()->node_pt(n)->position(tsteps, x);
				// get_initial_u(x, u);
				get_exact_u(time_pt()->time(tsteps), x, u);
				mesh_pt()->node_pt(n)->set_value(tsteps, 0, u[0]);
			}

			int step = 0;
			doc_solution(step, time_pt()->time(tsteps), tsteps);
			step++;

			for (int t = tsteps-1; t >= 0; t--) {
				double time = time_pt()->time((uint) t);

				update_interface(t, true);

				// set initial values
				for (unsigned long int n = 0; n < nnode; n++) {
					mesh_pt()->node_pt(n)->position(t, x);
					get_exact_u(time, x, u);
					mesh_pt()->node_pt(n)->set_value(t, 0, u[0]);
					// printf("[%6.4f] x=% 6.4f u=% 6.4f actual=%6.4f\n", time, x[0], u[0], mesh_pt()->node_pt(n)->value(t, 0));
				}

				printf("[% 4d] Setting initial condition at t=%8.6f\n", step, time);
				doc_solution(step, time, t);
				step++;
			}

			time_pt()->time() = t_shift;
		}

		void update_interface(const unsigned int &t = 0, bool ic = false) {
			Vector<double> s(2);
			Vector<double> flux(2);

			s[0] = 1.0;
			s[1] = 0.0;

			double tot_flux = 0.0;

			unsigned long int nelems = mesh_pt()->nboundary_element(4);
			for (unsigned long int e = 0; e < nelems; e++) {
				if (mesh_pt()->face_index_at_boundary(4, e) == 1) {
					EL *elem = dynamic_cast<EL *>(mesh_pt()->boundary_element_pt(4, e));
					if (ic) get_flux_ic(t+1, elem, s, flux);
					else 	elem->get_flux(s, flux);

					tot_flux += flux[0];
				}
			}

			double v = tot_flux / (alpha * Ny);
			double new_h = domain->get_interface() + v * dt;

			domain->set_interface(new_h);
			mesh_pt()->node_update();

			char filename[512];
			sprintf(filename, "%s/interface_velocity.dat", info.directory().c_str());
			FILE *file = fopen(filename, "a");
			fprintf(file, "%8.6f %8.6f %8.6f\n", time_pt()->time(t), new_h, v);
			fclose(file);

			printf("moved interface to x=%8.6f (v=%8.6f, flux=%8.6f)\n", new_h, v, tot_flux);
		}

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

		void doc_solution(uint timestep) {
			double time = time_pt()->time();
			doc_solution(timestep, time, 0);
		}
 
		void doc_solution(uint timestep, double time, const unsigned int &history_t) {
			cuint npts = 5;

			char filename[256];
			ofstream outfile;

			sprintf(filename, "%s/solns/soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			if (history_t == 0) {
				for (unsigned int e = 0; e < Nx*Ny; e++)
					dynamic_cast<EL *>(mesh_pt()->element_pt(e))->output(outfile, npts);
			} else {
				doc_primary(history_t, outfile, npts);
			}
			outfile.close();

			sprintf(filename, "%s/exact_solns/exact_soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			for (unsigned int e = 0; e < Nx*Ny; e++)
				dynamic_cast<EL *>(mesh_pt()->element_pt(e))->output_fct(outfile, npts, time, get_exact_u);
			outfile.close();

			double norm = 0.0;
			double error = 0.0;
			double err_tmp = 0.0;
			sprintf(filename, "%s/errors/error%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			if (history_t == 0) {
				for (unsigned int e = 0; e < Nx*Ny; e++) {
					EL *elem = dynamic_cast<EL *>(mesh_pt()->element_pt(e));
					elem->compute_error(outfile, get_exact_u, time, err_tmp, norm);
					error += err_tmp;
				}
			} else {
				for (unsigned int e = 0; e < Nx*Ny; e++) {
					EL *elem = dynamic_cast<EL *>(mesh_pt()->element_pt(e));
					compute_error_ic(history_t, elem, outfile, get_exact_u, time, err_tmp, norm);
					error += err_tmp;
				}
			}
			outfile.close();
			error = sqrt(error) / (Nx*Ny);

			double gamma = domain->get_interface() / (2 * sqrt(D2 * time));
			double D_eff = 4 * gamma * gamma * D1;

			printf("[%4u] time = %10.8f | error = %e | interface = %8.6f | gamma = %8.6f | D_eff = %8.6f\n", timestep, time, error, domain->get_interface(), gamma, D_eff);

			sprintf(filename, "%s/times.dat", info.directory().c_str());
			outfile.open(filename, ios::app);
			outfile << time << endl;
			outfile.close();

			sprintf(filename, "%s/interface.dat", info.directory().c_str());
			outfile.open(filename, ios::app);
			outfile << domain->get_interface() << endl;
			outfile.close();

			sprintf(filename, "%s/gamma.dat", info.directory().c_str());
			outfile.open(filename, ios::app);
			outfile << gamma << endl;
			outfile.close();

			info.number()++;
		}

	// modified from UnsteadyHeatEquations::output in unsteady_heat_elements.cc
	void doc_primary(const unsigned int &t, std::ostream& outfile, const unsigned int &nplot)
	{
		// Vector of local coordinates
		Vector<double> s(2);

		for (unsigned long int e = 0; e < Nx*Ny; e++) {
			// printf("Writing elem %lu\n", e);
			EL *elem = dynamic_cast<EL *>(mesh_pt()->element_pt(e));
			// Tecplot header info
			outfile << elem->tecplot_zone_string(nplot);
	
			// Loop over plot points
			unsigned num_plot_points = elem->nplot_points(nplot);
			for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
			{
				// Get local coordinates of plot point
				elem->get_s_plot(iplot, nplot, s);
	
				for (unsigned i = 0; i < 2; i++)
				{
				outfile << elem->interpolated_x(t, s, i) << " ";
				}
				outfile << elem->interpolated_u_ust_heat(t, s) << std::endl;
			}
	
			// Write tecplot footer (e.g. FE connectivity lists)
			elem->write_tecplot_zone_footer(outfile, nplot);
		}
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
      for (unsigned j = 0; j < 2; j++)
      {
        flux[j] = 0.0;
      }

      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over derivative directions
        for (unsigned j = 0; j < 2; j++)
        {
          flux[j] += elem->nodal_value(t, l, u_nodal_index) * dpsidx(l, j);
        }
      }
    }

	// modified from UnsteadyHeatEquations<2>::compute_error from unsteady_heat_elements.cc
	void compute_error_ic(
		const unsigned int &t,
		EL *elem,
		std::ostream& outfile,
		FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
		const double& time,
		double& error,
		double& norm) 
	{
		// Initialise
		error = 0.0;
		norm = 0.0;
		unsigned int DIM = 2;

		// Vector of local coordinates
		Vector<double> s(DIM);

		// Vector for coordintes
		Vector<double> x(DIM);

		// Find out how many nodes there are in the element
		unsigned n_node = elem->nnode();

		Shape psi(n_node);

		// Set the value of n_intpt
		unsigned n_intpt = elem->integral_pt()->nweight();

		// Tecplot
		outfile << "ZONE" << std::endl;

		// Exact solution Vector (here a scalar)
		Vector<double> exact_soln(1);

		// Loop over the integration points
		for (unsigned ipt = 0; ipt < n_intpt; ipt++)
		{
			// Assign values of s
			for (unsigned i = 0; i < DIM; i++)
			{
				s[i] = elem->integral_pt()->knot(ipt, i);
			}

			// Get the integral weight
			double w = elem->integral_pt()->weight(ipt);

			// Get jacobian of mapping
			double J = elem->J_eulerian(s);

			// Premultiply the weights and the Jacobian
			double W = w * J;

			// Get x position as Vector
			elem->interpolated_x(s, x);

			// Get FE function value
			double u_fe = elem->interpolated_u_ust_heat(t, s);

			// Get exact solution at this point
			(*exact_soln_pt)(time, x, exact_soln);

			// Output x,y,...,error
			for (unsigned i = 0; i < DIM; i++)
			{
				outfile << x[i] << " ";
			}
			outfile << exact_soln[0] << " " << exact_soln[0] - u_fe << std::endl;

			// Add to error and norm
			norm += exact_soln[0] * exact_soln[0] * W;
			error += (exact_soln[0] - u_fe) * (exact_soln[0] - u_fe) * W;
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

	CommandLineArgs::specify_command_line_flag("--nx", &Nx);
	CommandLineArgs::specify_command_line_flag("--ny", &Ny);
	CommandLineArgs::specify_command_line_flag("--tsteps", &t_steps);
	CommandLineArgs::specify_command_line_flag("--dt", &dt);
	CommandLineArgs::specify_command_line_flag("--tshift", &t_shift);
	CommandLineArgs::specify_command_line_flag("--vardt");
	CommandLineArgs::specify_command_line_flag("--write-freq", &write_freq);

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

	char dname[256];
	sprintf(dname, "RESLT/%dn%u_%dt%.2e", X_ORDER, Nx, T_ORDER, dt);
	printf("Saving results to %s\n", dname);
	
	if (std::filesystem::exists(dname)) {
		cout << dname << " exists" << endl;
	} else {
		std::filesystem::create_directories(dname);
	}

	char sub_dname[512];
	sprintf(sub_dname, "%s/errors", dname);
	if (std::filesystem::exists(sub_dname)) {
		cout << sub_dname << " exists" << endl;
	} else {
		std::filesystem::create_directories(sub_dname);
	}

	sprintf(sub_dname, "%s/solns", dname);
	if (std::filesystem::exists(sub_dname)) {
		cout << sub_dname << " exists" << endl;
	} else {
		std::filesystem::create_directories(sub_dname);
	}

	sprintf(sub_dname, "%s/exact_solns", dname);
	if (std::filesystem::exists(sub_dname)) {
		cout << sub_dname << " exists" << endl;
	} else {
		std::filesystem::create_directories(sub_dname);
	}

	DocInfo info;
	info.set_directory(dname);
	info.number() = 0;

	cout << "Output directory: " << info.directory() << endl;

	if (t_shift == 0.0) {
		t_shift = T_ORDER * dt;
	}

	// Ny = Nx;

	auto problem = Erf2D6Problem<RefineableQUnsteadyHeatElement<2,X_ORDER>>(Nx, Ny, t_steps, dt, t_shift, info);

	problem.initialise_dt(dt);
	problem.set_initial_conditions();

	// problem.doc_solution(0);

	char config_fname[256];
	sprintf(config_fname, "%s/config", dname);

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
	fclose(file);

	int prev_steps = problem.time_stepper_pt()->nprev_values()+1;

	for (uint t = 0; t < t_steps; t++) {
		problem.unsteady_newton_solve(dt);

		if (t % write_freq == 0)
			problem.doc_solution(t + prev_steps);

		double x_int = domain->get_interface();

		if (x_int >= x2 || isnan(x_int)) {
			printf("Interface reached right boundary, exiting...\n");
			break;
		}
	}
}