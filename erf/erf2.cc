#include <filesystem>

#include "includes.h"
#include "solver.h"
#include "erf2.h"

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

	
	auto problem_def = Erf2Problem(Nx, Nx, t_steps, dt, t_shift);
	
	cout << problem_def.Nx1 << " " << problem_def.Nx2 << " " << problem_def.t_steps << " " << problem_def.dt << " " << problem_def.t_shift << endl;
	
	Solver problem = Solver<QUnsteadyHeatElement<1, X_ORDER>, Erf2Problem>(problem_def, info);

	problem.initialise_dt(problem_def.dt);
	problem.set_initial_conditions();
	problem.quiet_solve();

	problem.doc_solution(0);

	char config_fname[256];
	sprintf(config_fname, "%s/config", dname);

	FILE *file = fopen(config_fname, "w");
	fprintf(file, "nx1=%d\n", problem_def.Nx1);
	fprintf(file, "nx2=%d\n", problem_def.Nx2);
	fprintf(file, "nx=%d\n", problem_def.Nx1+problem_def.Nx2);
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