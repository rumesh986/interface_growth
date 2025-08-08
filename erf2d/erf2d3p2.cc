// trial code for three phase 2d simulation
// one interface moves with prescribed velocity

#include <filesystem>

#include "includes.h"
#include "three_layer_mesh.h"

#ifndef RUN_SCRIPT
#define X_ORDER 2
#define T_ORDER 1
#endif

double D1 = 1.0;
double D2 = 0.5;
double D3 = 2.0;

static const double Ts = -1.0;
static const double Tm = 0.0;
static const double Tfr = 1.0;

static const double alpha = 10.0;

static const double x0 = -1.0;
static const double x1 = 0.0;
static double x2 = 0.1;
static const double x3 = 1.0;

static const double y_0 = 0.0;
static const double y_1 = 1.0;

static const double vel = 1.0;

TwoPhaseDomain *domain;

const Vector<unsigned int> analytical_boundaries = {
	1,
	3,
	// 4
};

void get_exact_u(const double &t, const Vector<double> &x, Vector<double> &u) {
	double h = domain->get_interface();
	if (x[0] < h) {
		u[0] = Ts + (Tm-Ts)/(1+erf(h/(2*sqrt(D1*t)))) * (1 + erf(x[0]/(2*sqrt(D1*t))));
	} else {
		u[0] = Tfr - (Tfr-Tm)/(1-erf(h/(2*sqrt(D2*t)))) * (1 - erf(x[0]/(2*sqrt(D2*t))));
	}
}

void get_source(const double &t, const Vector<double> &x, double &source) {
	source = 0.0;
}

template<class EL>
class Erf2D3P2Problem : public Problem {
	private:
		uint nx1, nx2, nx3, nx, ny;
		uint t_steps;
		double dt, t_shift;

		DocInfo info;
	public:
		Erf2D3P2Problem(
			uint Nx,
			uint Ny,
			uint t_steps,
			double dt,
			double t_shift,
			DocInfo info
		) : nx1(Nx), nx2(Nx), nx3(Nx), nx(nx1+nx2+nx3), ny(Ny), t_steps(t_steps), dt(dt), t_shift(t_shift), info(info) {
			add_time_stepper_pt(new BDF<T_ORDER>);

			x2 = x2 + vel * t_shift;

			mesh_pt() = new RefineableThreePhase2DMesh<EL>(nx1, nx2, nx3, ny, x0, x1, x2, x3, y_0, y_1, vel, time_pt(), time_stepper_pt());
			mesh_pt()->setup_boundary_element_info();
			domain = dynamic_cast<RefineableThreePhase2DMesh<EL> *>(mesh_pt())->domain;

			if (domain == nullptr) {
				printf("Domain got problems\n");
			}

			unsigned long int nelems = mesh_pt()->nelement();
			for (unsigned int e = 0; e < nelems; e++)
				dynamic_cast<EL *>(mesh_pt()->element_pt(e))->source_fct_pt() = get_source;


			for (unsigned int b : analytical_boundaries) {
				unsigned long int nnode = mesh_pt()->nboundary_node(b);
				for (unsigned long int n = 0; n < nnode; n++)
					mesh_pt()->boundary_node_pt(b, n)->pin_all();
			}

			for (unsigned int yi = 0; yi < ny; yi++) {
				for (unsigned int xi = 0; xi < nx1; xi++) {
					dynamic_cast<EL *>(mesh_pt()->element_pt(yi*nx + xi))->beta_pt() = &D1;
					printf("yi=%2u xi=%2u beta=%4.2f D=%4.3f\n", yi, xi, *(dynamic_cast<EL *>(mesh_pt()->element_pt(yi*nx + xi))->beta_pt()), D1);
				}
				
				for (unsigned int xi = nx1; xi < nx1+nx2; xi++) {
					dynamic_cast<EL *>(mesh_pt()->element_pt(yi*nx + xi))->beta_pt() = &D2;
					printf("yi=%2u xi=%2u beta=%4.2f D=%4.3f\n", yi, xi, *(dynamic_cast<EL *>(mesh_pt()->element_pt(yi*nx + xi))->beta_pt()), D2);
				}
				
				for (unsigned int xi = nx1+nx2; xi < nx1+nx2+nx3; xi++) {
					dynamic_cast<EL *>(mesh_pt()->element_pt(yi*nx + xi))->beta_pt() = &D3;
					printf("yi=%2u xi=%2u beta=%4.2f D=%4.3f\n", yi, xi, *(dynamic_cast<EL *>(mesh_pt()->element_pt(yi*nx + xi))->beta_pt()), D3);
				}
			}

			assign_eqn_numbers();
			linear_solver_pt()->disable_doc_time();
			disable_info_in_newton_solve();

			initialise_dt(dt);
			set_initial_conditions();
			doc_solution(0);
		}

		~Erf2D3P2Problem() {
			delete mesh_pt();
		}

		void actions_before_newton_solve() {}
		void actions_after_newton_solve() {}

		void actions_before_implicit_timestep() {
			mesh_pt()->node_update();
			
			Vector<double> x(2);
			Vector<double> u(1);
			double cur_t = time_pt()->time();

			for (unsigned int b : analytical_boundaries) {
				unsigned long int nnode = mesh_pt()->nboundary_node(b);
				for (unsigned long int n = 0; n < nnode; n++) {
					Node *node_pt = mesh_pt()->boundary_node_pt(b, n);
					node_pt->position(x);
					get_exact_u(cur_t, x, u);
					node_pt->set_value(0, u[0]);
				}
			}
		}

		void actions_after_implicit_timestep() {}

		void set_initial_conditions() {
			time_pt()->time() = t_shift;
			unsigned long int nnode = mesh_pt()->nnode();

			// set node positions in previous timesteps
			// oomph lib defaults to only setting node positions (x[]) for current timestep (t = 0)
			// the other timesteps (t > 0) have positions defaulted to 0 (x[] = 0)
			// This is Technically wrong for moving interface problems -> need to set past positions properly.
			for (unsigned long int n = 0; n < nnode; n++) {
				mesh_pt()->node_pt(n)
				->position_time_stepper_pt()
				->assign_initial_positions_impulsive(mesh_pt()->node_pt(n));
			}
			
			// actually set initial conditions
			Vector<double> x(2);
			Vector<double> u(1);
			// use int instead of unsigned to avoid segfault when t goes below zero
			// happens after last run through loop
			for (int t = time_stepper_pt()->nprev_values(); t>= 0; t--) {
				double cur_t = time_pt()->time((unsigned int) t);
				printf("[% 2d] Setting initial condition at t=%8.6f\n", t, cur_t);

				for (unsigned long int n = 0; n < nnode; n++) {
					Node *node_pt = mesh_pt()->node_pt(n);
					node_pt->position(t, x);
					get_exact_u(cur_t, x, u);
					node_pt->set_value(t, 0, u[0]);
				}
			}

			time_pt()->time() = t_shift;

		}

		void doc_solution(uint timestep) {
			cuint npts = 5;

			char filename[100];
			ofstream outfile;

			double time = time_pt()->time();// + problem.t_shift;

			sprintf(filename, "%s/soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			for (uint e = 0; e < nx; e++) {
				EL * el_pt = dynamic_cast<EL *>(mesh_pt()->element_pt(e));
				el_pt->output(outfile, npts);
			}
			outfile.close();

			sprintf(filename, "%s/exact_soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			for (uint e = 0; e < nx; e++) {
				EL * el_pt = dynamic_cast<EL *>(mesh_pt()->element_pt(e));
				el_pt->output_fct(outfile, npts, time, get_exact_u);
			}
			outfile.close();

			double error, norm;
			sprintf(filename, "%s/error%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			for (uint e = 0; e < nx; e++) {
				EL * el_pt = dynamic_cast<EL *>(mesh_pt()->element_pt(e));
				el_pt->compute_error(outfile, get_exact_u, time, error, norm);
			}
			outfile.close();

			printf("[%4u] time = %10.8f | error = %e | interface = %8.6f\n", timestep, time_pt()->time(), error, domain->get_interface());

			sprintf(filename, "%s/times.dat", info.directory().c_str());
			outfile.open(filename, ios::app);
			outfile << time << endl;
			outfile.close();

			sprintf(filename, "%s/interface.dat", info.directory().c_str());
			outfile.open(filename, ios::app);
			outfile << domain->get_interface() << endl;
			outfile.close();

			info.number()++;

			if (timestep == 0) {
				for (unsigned int yi = 0; yi < ny; yi++) {
					for (unsigned int xi = 0; xi < nx1; xi++)
						printf("yi=% 2u xi=% 2u beta=%5.3f exp=%5.3f\n", yi, xi, *(dynamic_cast<EL *>(mesh_pt()->element_pt(yi*nx + xi))->beta_pt()), D1);
					
					for (unsigned int xi = nx1; xi < nx1+nx2; xi++)
						printf("yi=% 2u xi=% 2u beta=%5.3f exp=%5.3f\n", yi, xi, *(dynamic_cast<EL *>(mesh_pt()->element_pt(yi*nx + xi))->beta_pt()), D2);
					
					for (unsigned int xi = nx1+nx2; xi < nx1+nx2+nx3; xi++)
						printf("yi=% 2u xi=% 2u beta=%5.3f exp=%5.3f\n", yi, xi, *(dynamic_cast<EL *>(mesh_pt()->element_pt(yi*nx + xi))->beta_pt()), D3);
				}
			}
		}
};

int main(int argc, char *argv[]) {
	uint Nx = 0;
	uint Ny = 3;
	uint t_steps = 100;
	double dt = 0.0;
	double t_shift = 0.0;
	bool var_dt = false;
	uint write_freq = 1;

	CommandLineArgs::setup(argc, argv);
	CommandLineArgs::specify_command_line_flag("--nx", &Nx);
	CommandLineArgs::specify_command_line_flag("--ny", &Nx);
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

	Ny = Nx;

	auto problem = Erf2D3P2Problem<RefineableQUnsteadyHeatElement<2, X_ORDER>>(Nx, Ny, t_steps, dt, t_shift, info);

	char config_fname[256];
	sprintf(config_fname, "%s/config", dname);

	FILE *file = fopen(config_fname, "w");
	fprintf(file, "nx1=%d\n", Nx);
	fprintf(file, "nx2=%d\n", Nx);
	fprintf(file, "nx3=%d\n", Nx);
	fprintf(file, "nx=%d\n", Nx+Nx+Nx);
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
	}
}