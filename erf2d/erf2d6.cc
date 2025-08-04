
// moving interface 
// interface moves by calculating instantaneous mesh velocity from flux
// right phase pinned to melting temperature

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
static const double alpha = 100.0;

TwoPhaseDomain *domain;

// change this 
void get_exact_u(const double &t, const Vector<double> &x, Vector<double> &u) {
	if (x[0] <= domain->get_interface()) {
		u[0] = Ts + (Tm - Ts) * erf(0.5 * x[0]/sqrt(D1 * t)) / erf(0.5 * domain->get_interface() / sqrt(D1 * t));
	} else {
		u[0] = Tm;
	}
}

template<class EL>
class Erf2D6Problem : public Problem {
	private:
		uint Nx1, Nx2, Nx, Ny;
		uint t_steps;
		double dt, t_shift;
		static constexpr double Lx = 1.0;
		static constexpr double Ly = 1.0;
		static constexpr double vel = 1.0;

		// static constexpr double D1 = 1.0;
		// static constexpr double D2 = 0.5;
		// static constexpr double Ts = -1.0;
		static constexpr double Tfr = 1.0;
		// static constexpr double Tm = 0.0;
		// static constexpr double alpha = 1.0;

		// TwoPhaseDomain *domain;
		DocInfo info;

	public:
		Erf2D6Problem(uint nx, uint ny,  uint t_steps, double dt, double t_shift, DocInfo info)
		: Nx1(nx), Nx2(nx), Nx(nx+nx), Ny(ny), t_steps(t_steps), dt(dt), t_shift(t_shift), info(info) {
			add_time_stepper_pt(new BDF<T_ORDER>);

			domain = new TwoPhaseDomain(0.0, 0.1, Lx, 0.0, Nx1, Nx2, Ny, time_pt());

			mesh_pt() = new RefineableTwoLayer2DMesh<EL>(Nx1, Nx2, Ny, 0.0, 0.1, Lx, 0.0, Ly, domain, time_stepper_pt());
			mesh_pt()->setup_boundary_element_info();

			for (uint e = 0; e < mesh_pt()->nelement(); e++)
				dynamic_cast<EL *>(mesh_pt()->element_pt(e))->source_fct_pt() = get_source;
			
			// set boundary 1 (right)
			for (uint n = 0; n < mesh_pt()->nboundary_node(1); n++) {
				mesh_pt()->boundary_node_pt(1, n)->set_value(0, Tm);
				mesh_pt()->boundary_node_pt(1, n)->pin(0);
			}

			// set boundary 3 (left)
			for (uint n = 0; n < mesh_pt()->nboundary_node(3); n++) {
				mesh_pt()->boundary_node_pt(3, n)->set_value(0, Ts);
				mesh_pt()->boundary_node_pt(3, n)->pin(0);
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
			mesh_pt()->node_update();
		}

		void actions_after_implicit_timestep() {
			// domain->set_interface(0.1); // change this
			Vector<double> flux(2);
			Vector<double> s(2);

			// local coordinates in element
			// set to middle of right hand side
			s[0] = 1.0;
			s[1] = 0.0;

			double tot_flux = 0.0;

			uint nelems = mesh_pt()->nboundary_element(4);
			for (uint e = 0; e < nelems; e++) {
				EL *elem = dynamic_cast<EL *>(mesh_pt()->boundary_element_pt(4, e));

				elem->get_flux(s, flux);

				tot_flux += flux[0];

			}

			printf("Total flux at bound4 to the right = %8.6f\n", tot_flux);
			double cur_x1 = domain->get_interface();
			printf("moving interface from %8.6f to ", cur_x1);
			double new_x1 = cur_x1 +  tot_flux/alpha * dt;
			printf("%8.6f\n", new_x1);

			domain->set_interface(new_x1);

			Vector<double> x(2), u(1);

			for (uint yi = 0; yi < Ny; yi++) {
				for (uint e = Nx1; e < Nx; e++) {
					EL *elem = dynamic_cast<EL *>(mesh_pt()->element_pt(yi*Nx + e));
					for (uint n = 0; n < elem->nnode(); n++)
						elem->node_pt(n)->set_value(0, Tm);
				}
			}
		}

		void set_initial_conditions() {
			time_pt()->time() = t_shift;
			int nnode = mesh_pt()->nnode();

			Vector<double> x(2);
			Vector<double> u(1);

			for (int n = 0; n < nnode; n++) {
				mesh_pt()->node_pt(n)
					->position_time_stepper_pt()
					->assign_initial_positions_impulsive(mesh_pt()->node_pt(n));
			}

			for (int t = time_stepper_pt()->nprev_values(); t >= 0; t--) {
				double cur_t = time_pt()->time((uint) t);
				printf("[% 2d] Setting initial condition at t = %10.8f\n", t, cur_t);

				for (int n = 0; n < nnode; n++) {
					x[0] = mesh_pt()->node_pt(n)->x(0);
					x[1] = mesh_pt()->node_pt(n)->x(1);

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

			printf("[%4u] time = %10.8f | error = %e\n", timestep, time_pt()->time(), error);

			sprintf(filename, "%s/times.dat", info.directory().c_str());
			outfile.open(filename, ios::app);
			outfile << time << endl;
			outfile.close();

			sprintf(filename, "%s/interface.dat", info.directory().c_str());
			outfile.open(filename, ios::app);
			outfile << domain->get_interface() << endl;
			outfile.close();

			info.number()++;

			EL *elem = dynamic_cast<EL *>(mesh_pt()->element_pt(Nx1-2));
			for (uint n = 0; n < elem->nnode(); n++) {
				printf("node %u x=%8.6f y=%8.6f beta=%8.6f\n", n, elem->node_pt(n)->x(0), elem->node_pt(n)->value(0), *(elem->beta_pt()));
			}

			elem = dynamic_cast<EL *>(mesh_pt()->element_pt(Nx1-1));
			for (uint n = 0; n < elem->nnode(); n++) {
				printf("node %u x=%8.6f y=%8.6f beta=%8.6f\n", n, elem->node_pt(n)->x(0), elem->node_pt(n)->value(0), *(elem->beta_pt()));
			}
		}

};

int main(int argc, char **argv) {
	CommandLineArgs::setup(argc, argv);

	uint Nx = 0;
	uint Ny = 3;
	uint t_steps = 100;
	double dt = 0.0;
	double t_shift = 0.0;
	bool var_dt = false;
	uint write_freq = 1;

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
	}
}