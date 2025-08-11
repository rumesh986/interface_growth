
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
static const double x1 = 0.0;
static const double x2 = 1.0;

static const double y_0 = 0.0;
static const double y_1 = 1.0;

// static const double Lx = 1.0;
// static const double Ly = 1.0;
static const double vel = 0.0;

unsigned int num_y = 0;

TwoPhaseDomain *domain;

void get_exact_u(const double &t, const Vector<double> &x, Vector<double> &u) {
	double h = domain->get_interface();
	if (x[0] < h) {
		u[0] = Ts + (Tm-Ts)/(1+erf(h/(2*sqrt(D1*t)))) * (1 + erf(x[0]/(2*sqrt(D1*t))));
	} else {
		u[0] = Tfr - (Tfr-Tm)/(1-erf(h/(2*sqrt(D2*t)))) * (1 - erf(x[0]/(2*sqrt(D2*t))));
	}
}

void dhdt(const double &t, const Vector<double> &x, Vector<double> &v) {
	double h = domain->get_interface();
	v[0] = 2*(Tm-Ts)/(alpha * sqrt(Pi*D1*t) * (1 + erf(h/(2*sqrt(D1*t))))) * exp(-1 * h*h/(4*D1*t));
	printf("t=%8.6f h=%8.6f v=%8.6f\n", t, h, v[0]);
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

			num_y = ny;
			domain = new TwoPhaseDomain(x0, x1, x2, vel, Nx1, Nx2, Ny, time_pt());

			mesh_pt() = new RefineableTwoLayer2DMesh<EL>(Nx1, Nx2, Ny, x0, x1, x2, y_0, y_1, domain, time_stepper_pt());
			mesh_pt()->setup_boundary_element_info();

			for (uint e = 0; e < mesh_pt()->nelement(); e++)
				dynamic_cast<EL *>(mesh_pt()->element_pt(e))->source_fct_pt() = get_source;

			// pin boundary 1 (right)
			for (uint n = 0; n < mesh_pt()->nboundary_node(1); n++) {
				mesh_pt()->boundary_node_pt(1, n)->pin(0);
			}

			// pin boundary 3 (left)
			for (uint n = 0; n < mesh_pt()->nboundary_node(3); n++) {
				mesh_pt()->boundary_node_pt(3, n)->pin(0);
			}

			// pin boundary 4 (interface)
			for (uint n = 0; n < mesh_pt()->nboundary_node(4); n++) {
				mesh_pt()->boundary_node_pt(4, n)->pin(0);
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
			Vector<double> flux(2);
			Vector<double> s(2);

			// local coordinates in element
			// set to middle of right hand side
			s[0] = 1.0;
			s[1] = 0.0;

			double tot_flux = 0.0;

			uint nelems = mesh_pt()->nboundary_element(4);
			for (uint e = 0; e < nelems; e++) {
				if (mesh_pt()->face_index_at_boundary(4, e) < 0)
					continue;

				dynamic_cast<EL *>(mesh_pt()->boundary_element_pt(4, e))->get_flux(s, flux);
				tot_flux += flux[0];
			}

			// take average of total flux across the y direction
			// to spread the flux across all the elements
			double new_x1 = domain->get_interface() +  tot_flux/(alpha*Ny) * dt;

			domain->set_interface(new_x1);
			mesh_pt()->node_update();

			// reset analytical boundaries
			Vector<double> x(2);
			Vector<double> u(1);

			double cur_t = time_pt()->time();

			// set boundary 1 (right)
			for (uint n = 0; n < mesh_pt()->nboundary_node(1); n++) {
				x[0] = mesh_pt()->boundary_node_pt(1, n)->x(0);
				get_exact_u(cur_t, x, u);
				mesh_pt()->boundary_node_pt(1, n)->set_value(0, u[0]);
			}

			// set boundary 3 (left)
			for (uint n = 0; n < mesh_pt()->nboundary_node(3); n++) {
				x[0] = mesh_pt()->boundary_node_pt(3, n)->x(0);
				get_exact_u(cur_t, x, u);
				mesh_pt()->boundary_node_pt(3, n)->set_value(0, u[0]);
			}

			// set boundary 4 (interface)
			for (uint n = 0; n < mesh_pt()->nboundary_node(4); n++) {
				x[0] = mesh_pt()->boundary_node_pt(4, n)->x(0);
				get_exact_u(cur_t, x, u);
				mesh_pt()->boundary_node_pt(4, n)->set_value(0, u[0]);
			}
		}

		void actions_after_implicit_timestep() {}

		void set_initial_conditions() {
			time_pt()->time() = t_shift;
			unsigned long int nnode = mesh_pt()->nnode();

			Vector<double> x(2);
			Vector<double> u(1);
			Vector<double> v(1);

			for (unsigned long int n = 0; n < nnode; n++) {
				mesh_pt()->node_pt(n)
					->position_time_stepper_pt()
					->assign_initial_positions_impulsive(mesh_pt()->node_pt(n));
			}

			for (int t = time_stepper_pt()->nprev_values(); t >= 0; t--) {
				double cur_t = time_pt()->time((uint) t);
				dhdt(cur_t, x, v);
				double new_x1 = domain->get_interface() + v[0] * dt;
				domain->set_interface(new_x1);
				mesh_pt()->node_update();
				printf("Setting interface at t=%8.6f to %8.6f (v=%8.6f)\n", cur_t, new_x1, v[0]);
			}
			
			for (int t = time_stepper_pt()->nprev_values(); t >= 0; t--) {
			// for (int t = 0; t < 1; t++) {
				double cur_t = time_pt()->time((uint) t);
				printf("[% 2d] Setting initial condition at t = %10.8f\n", t, cur_t);

				for (unsigned long int n = 0; n < nnode; n++) {
					x[0] = mesh_pt()->node_pt(n)->x(t, 0);
					x[1] = mesh_pt()->node_pt(n)->x(t, 1);

					get_exact_u(cur_t, x, u);
					mesh_pt()->node_pt(n)->set_value(t, 0, u[0]);
				}
			}

			time_pt()->time() = t_shift;
		};

		void doc_solution(uint timestep) {
			cuint npts = 5;

			char filename[100];
			ofstream outfile;

			double time = time_pt()->time();// + problem.t_shift;

			sprintf(filename, "%s/soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			for (uint e = 0; e < Nx; e++) {
				EL * el_pt = dynamic_cast<EL *>(mesh_pt()->element_pt(e));
				el_pt->output(outfile, npts);
			}
			outfile.close();

			sprintf(filename, "%s/exact_soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			for (uint e = 0; e < Nx; e++) {
				EL * el_pt = dynamic_cast<EL *>(mesh_pt()->element_pt(e));
				el_pt->output_fct(outfile, npts, time, get_exact_u);
			}
			outfile.close();

			double error, norm;
			sprintf(filename, "%s/error%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			for (uint e = 0; e < Nx; e++) {
				EL * el_pt = dynamic_cast<EL *>(mesh_pt()->element_pt(e));
				el_pt->compute_error(outfile, get_exact_u, time, error, norm);
			}
			outfile.close();

			double gamma = domain->get_interface() / (2 * sqrt(D2 * time));

			printf("[%4u] time = %10.8f | error = %e | interface = %8.6f | gamma = %8.6f\n", timestep, time_pt()->time(), error, domain->get_interface(), gamma);

			sprintf(filename, "%s/times.dat", info.directory().c_str());
			outfile.open(filename, ios::app);
			outfile << time << endl;
			outfile.close();

			sprintf(filename, "%s/interface.dat", info.directory().c_str());
			outfile.open(filename, ios::app);
			outfile << domain->get_interface() << endl;
			outfile.close();

			info.number()++;
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

	DocInfo info;
	info.set_directory(dname);
	info.number() = 0;

	cout << "Output directory: " << info.directory() << endl;

	if (t_shift == 0.0) {
		t_shift = (T_ORDER + 1) * dt;
	}

	// Ny = Nx;

	auto problem = Erf2D6Problem<RefineableQUnsteadyHeatElement<2,X_ORDER>>(Nx, Ny, t_steps, dt, t_shift, info);

	problem.initialise_dt(dt);
	problem.set_initial_conditions();

	problem.doc_solution(0);

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

	for (uint t = 0; t < t_steps; t++) {
		problem.unsteady_newton_solve(dt);

		if (t % write_freq == 0)
			problem.doc_solution(t);

		double x_int = domain->get_interface();

		if (x_int >= x2 || isnan(x_int)) {
			printf("Interface reached right boundary, exiting...\n");
			break;
		}
	}
}