#include <filesystem>

#include "includes.h"
#include "solver.h"
#include "erf2d1.h"

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

	
	auto problem_def = Erf2D1Problem(Nx, Nx, Ny, t_steps, dt, t_shift);
	
	cout << problem_def.Nx << " " << problem_def.t_steps << " " << problem_def.dt << " " << problem_def.t_shift << endl;
	
	Solver problem = Solver<QUnsteadyHeatElement<2, X_ORDER>, Erf2D1Problem>(problem_def, info);

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
	fprintf(file, "write_freq=%u\n", write_freq);
	fclose(file);

	for (uint t = 0; t < problem_def.t_steps; t++) {
		problem.unsteady_newton_solve(problem_def.dt);

		if (t % write_freq == 0)
			problem.doc_solution(t);
	}
}