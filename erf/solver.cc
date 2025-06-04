#include "includes.h"

#include <filesystem>

#include "erf1.h"

template<class EL, class P> class Solver : public Problem {
	private:
		P problem;
		DocInfo info;
	
	public:
		Solver(P problem, DocInfo info) : problem(problem), info(info) {
			// add time stepper
			add_time_stepper_pt(new BDF<T_ORDER>);

			// create basic mesh and set up boundary information
			mesh_pt() = new OneDMesh<EL>(problem.Nx, problem.Lx, time_stepper_pt());
			mesh_pt()->setup_boundary_element_info();

			for (uint e = 0; e < mesh_pt()->nelement(); e++)
				dynamic_cast<EL *>(mesh_pt()->element_pt(e))->source_fct_pt() = problem.get_source;

			for (auto iter = problem.flux_boundaries.begin(); iter != problem.flux_boundaries.end(); iter++) {
				cout << "Creating flux elements at boundary " << iter->first << endl;
				create_flux_elements(iter->first, iter->second);
			}

			for (auto iter = problem.pinned_boundaries.begin(); iter != problem.pinned_boundaries.end(); iter++) {
				cout << "Assigning pinned boundary condition at " << iter->first << endl;
				int nnode = mesh_pt()->nboundary_node(iter->first);
				for (int n = 0; n < nnode; n++) {
					mesh_pt()->boundary_node_pt(iter->first, n)->set_value(0, iter->second);
					mesh_pt()->boundary_node_pt(iter->first, n)->pin(0);
				}
			}

			for (int i = 0; i < problem.analytical_boundaries.size(); i++) {
				uint b = problem.analytical_boundaries[i];
				for (uint n = 0; n < mesh_pt()->nboundary_node(b); n++) {
					mesh_pt()->boundary_node_pt(b, n)->pin(0);
				}
			}

			cout << "Number of equations " << assign_eqn_numbers() << endl;
		}

		~Solver() {
			delete mesh_pt();
		}

		void actions_before_newton_solve() {}
		void actions_after_newton_solve() {}
		
		void actions_before_implicit_timestep() {
			double cur_t = time_pt()->time();

			Vector<double> x(1);
			Vector<double> u(1);

			for (uint i = 0; i < problem.analytical_boundaries.size(); i++) {
				uint b = problem.analytical_boundaries[i];
				for (uint n = 0; n < mesh_pt()->nboundary_node(b); n++) {
					x[0] = mesh_pt()->boundary_node_pt(b, n)->x(0);

					problem.get_exact_u(cur_t, x, u);
					mesh_pt()->boundary_node_pt(b, n)->set_value(0, u[0]);
				}
			}
		}

		void actions_after_implicit_timestep() {}

		void create_flux_elements(uint b, FluxFctPt& flux_ptr) {
			uint n_elems = mesh_pt()->nboundary_element(b);

			for (uint e = 0; e < n_elems; e++) {
				EL *elem = dynamic_cast<EL *>(mesh_pt()->boundary_element_pt(b, e));
				int face_index = mesh_pt()->face_index_at_boundary(b, e);

				UnsteadyHeatFluxElement<EL> *flux_elem = new UnsteadyHeatFluxElement<EL>(elem, face_index);
				flux_elem->flux_fct_pt() = flux_ptr;
				mesh_pt()->add_element_pt(flux_elem);
			}
		}

		void set_initial_conditions() {
			time_pt()->time() = problem.t_shift;
			int nnode = mesh_pt()->nnode();

			Vector<double> x(1);
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

					problem.get_exact_u(cur_t, x, u);
					mesh_pt()->node_pt(n)->set_value(t, 0, u[0]);
					printf("[%10.8f] x=%10.8f, u=%10.8f\n", cur_t, x[0], u[0]);
				}
			}

			time_pt()->time() = problem.t_shift;
		}

		void doc_solution(uint timestep) {
			cuint npts = 5;

			char filename[100];
			ofstream outfile;

			double time = time_pt()->time();// + problem.t_shift;

			sprintf(filename, "%s/soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);

			for (uint e = 0; e < problem.Nx; e++) {
				EL * el_pt = dynamic_cast<EL *>(mesh_pt()->element_pt(e));
				el_pt->output(outfile, npts);
			}

			outfile.close();

			sprintf(filename, "%s/exact_soln%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			for (uint e = 0; e < problem.Nx; e++) {
				EL * el_pt = dynamic_cast<EL *>(mesh_pt()->element_pt(e));
				el_pt->output_fct(outfile, npts, time, problem.get_exact_u);
			}
			outfile.close();

			double error, norm;
			sprintf(filename, "%s/error%i.dat", info.directory().c_str(), info.number());
			outfile.open(filename);
			for (uint e = 0; e < problem.Nx; e++) {
				EL * el_pt = dynamic_cast<EL *>(mesh_pt()->element_pt(e));
				el_pt->compute_error(outfile, &problem.get_exact_u, time, error, norm);
			}
			outfile.close();

			printf("[% 4u] time = %10.8f | error = %e\n", timestep, time_pt()->time(), error);

			sprintf(filename, "%s/times.dat", info.directory().c_str());
			outfile.open(filename, ios::app);
			outfile << time << endl;
			outfile.close();

			info.number()++;
		}

		void quiet_solve() {
			linear_solver_pt()->disable_doc_time();
			disable_info_in_newton_solve();
		}

};

int main(int argc, char **argv) {
	CommandLineArgs::setup(argc, argv);

	uint Nx = 0;
	uint t_steps = 100;
	double dt = 0.0;
	double t_shift = 0.0;
	bool var_dt = false;

	CommandLineArgs::specify_command_line_flag("--nx", &Nx);
	CommandLineArgs::specify_command_line_flag("--tsteps", &t_steps);
	CommandLineArgs::specify_command_line_flag("--dt", &dt);
	CommandLineArgs::specify_command_line_flag("--tshift", &t_shift);
	CommandLineArgs::specify_command_line_flag("--vardt");

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

	auto problem_def = Erf1Problem(Nx, t_steps, dt, t_shift);

	cout << problem_def.Nx << " " << problem_def.t_steps << " " << problem_def.dt << " " << problem_def.t_shift << endl;
	
	Solver problem = Solver<QUnsteadyHeatElement<1, X_ORDER>, Erf1Problem>(problem_def, info);

	problem.initialise_dt(problem_def.dt);
	problem.set_initial_conditions();
	problem.quiet_solve();

	problem.doc_solution(0);

	char config_fname[256];
	sprintf(config_fname, "%s/config", dname);

	FILE *file = fopen(config_fname, "w");
	fprintf(file, "nx=%d\n", problem_def.Nx);
	fprintf(file, "x_order=%d\n", X_ORDER);
	fprintf(file, "dt=%e\n", problem_def.dt);
	fprintf(file, "t_shift=%e\n", problem_def.t_shift);
	fprintf(file, "t_order=%d\n", T_ORDER);
	fprintf(file, "t_steps=%d\n", problem_def.t_steps);
	fclose(file);

	for (uint t = 0; t < problem_def.t_steps; t++) {
		problem.unsteady_newton_solve(problem_def.dt);
		problem.doc_solution(t);
	}

}