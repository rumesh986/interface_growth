
// moving interface
// calculate as part of solution in newton iterations (hopefully)
// move with spine meshes
// using proper non-dimensionalized form with non-dimensionalized constants
// and corrected temperature scaling
// still accepts dimensional inputs and outputs non-dimensional

#include <filesystem>

#include "includes.h"
#include "two_phase_free_boundary_mesh.h"

// dimensional parameters
double _k[2] = {1.0, 0.25}; // thermal conductivity
double _rho[2] = {1.0, 2.0}; // density
double _Cp[2] = {1.0, 0.25}; // specific heat
double L = 1.0; // latent heat of fusion

// the above parameters are used to determine the thermal diffusivity below
// the third value here is the analyical value of the effective diffusivity
double _D[3] = {0.0, 0.0, 0.0};

// dimensionless parameters
double D = _D[1] / _D[0];
double D2 = 0.0;
double k = _k[1] / _k[0];
double St = 0.0;

bool ic_set = false;

// boundary temperature values (dimensional)
static const double T_s = -1.0; // temp at x=0 (solid)
static const double T_m = 0.0; // melting point x=h(t)
static const double T_l = 1.0; // temp far away in liquid x->\infty

// simulation domain
double xs[3] = {0.0, 0.025, 2.0};
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
	double h = geometry->get_interface();
	double h_ana = sqrt(_D[2] * t);

	// dimensionless temperature far from interface
	double trans_Tl = (T_l - T_m) / (T_m - T_s);

	if (x[0] < h) {
		double denom = 0.5 / sqrt(t);
		u[0] = (erf(x[0] * denom) / erf(h * denom)) - 1.0;
	} else {
		double denom = 0.5 / sqrt(D * t);
		u[0] = trans_Tl * (erf(x[0] * denom) - erf(h * denom)) / (1.0 - erf(h * denom));
	}

	if (x[0] < h_ana) {
		double denom = 0.5 / sqrt(t);
		u[1] = (erf(x[0] * denom) / erf(h_ana * denom)) - 1.0;
	} else {
		double denom = 0.5 / sqrt(D * t);
		u[1] = trans_Tl * (erf(x[0] * denom) - erf(h_ana * denom)) / (1.0 - erf(h_ana * denom));
	}
}

Vector<unsigned int> analytical_boundaries = {
	1, // right
	// 3, // left
	// 4 //interface
};

Vector<unsigned int> pinned_boundaries = {
	3, // left
	4 // interface
};

map<unsigned int, FluxFctPt> flux_boundaries = {
	// {0, no_flux_fct},
	// {2, no_flux_fct},
	// {1, no_flux_fct},
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

			add_time_stepper_pt(new BDF<T_ORDER>(true));

			geometry = new FreeBoundaryElement(xs[0], xs[1], xs[2], ys[0], ys[1], St, k, time_stepper_pt());

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
			
			// set boundary conditions
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
				
				for (unsigned int e = nx1; e < nx; e++)
					dynamic_cast<EL *>(bulk_mesh_pt->element_pt(base + e))->beta_pt() = &D;
			}

			geometry->add_phase1_element(dynamic_cast<EL *>(bulk_mesh_pt->element_pt(nx1-1)));
			geometry->add_phase2_element(dynamic_cast<EL *>(bulk_mesh_pt->element_pt(nx1)));

			printf("Total number of equations: %lu\n", assign_eqn_numbers());
			printf("NDOF: %lu\n", ndof());

			linear_solver_pt()->disable_doc_time();
			disable_info_in_newton_solve();
			// newton_solver_tolerance() = 1e-7 * dt;
			newton_solver_tolerance() = 5e-10;
			max_newton_iterations() = 1e7;
			max_residuals() = 1e3;

			printf("tolerance for newton: %e\n", newton_solver_tolerance());
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

		void actions_before_newton_solve() {
			// double De_sqrt_estimate = geometry->get_interface() / sqrt(time_pt()->time()-dt);
			time_stepper_pt()->set_predictor_weights();
			time_stepper_pt()->calculate_predicted_values(geometry->geom_data_pt(0));
			
			// double h_est = De_sqrt_estimate * sqrt(time_pt()->time() - 0.5*dt);
			double h_pred = geometry->x1(time_stepper_pt()->predictor_storage_index());
			// double h_ana = sqrt(_D[2] * time_pt()->time());

			// printf("[BSolve] setting interface estimate to %16.14f at time=%8.6f (analytical: %16.14f, predicted: %16.14f)\n", h_est, time_pt()->time(), h_ana, h_pred);
			// printf("\terror in estimate: %e\n\terror in prediction: %e\n", fabs(h_ana - h_est), fabs(h_ana - h_pred));

			geometry->set_interface(h_pred);

			for (unsigned s = 0; s < bulk_mesh_pt->nspine(); s++) {
				bulk_mesh_pt->spine_pt(s)->height() = geometry->x1();
			}

			bulk_mesh_pt->node_update();
		};
		void actions_after_newton_solve() {};

		void actions_before_implicit_timestep() {
			Vector<double> x(2), u(2);

			double time = time_pt()->time();

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

		void actions_before_newton_step() {};
		
		// calculate flux and store flux
		void actions_before_newton_convergence_check() {
			// Vector<double> flux(2);

			// double tot_flux = 0.0;
			// unsigned long int nelems = bulk_mesh_pt->nboundary_element(4);

			// Vector<double> s(2);
			// s[1] = 0.0;
			
			// for (unsigned long int e = 0; e < nelems; e++) {
			// 	int face_index = bulk_mesh_pt->face_index_at_boundary(4, e);
			// 	EL *elem = dynamic_cast<EL *>(bulk_mesh_pt->boundary_element_pt(4, e));
				
			// 	double factor = 0.0;
			// 	if (face_index == 1) {
			// 		s[0] = 1.0;
			// 		factor = 1.0;
			// 	} else if (face_index == -1) {
			// 		s[0] = -1.0;
			// 		factor = -k;
			// 	}

			// 	elem->get_flux(s, flux);
			// 	tot_flux += factor * flux[0];
			// }
			
			// // double latent_est = St * time_stepper_pt()->time_derivative(1, geometry->internal_data_pt(0), 1);
			// // printf("[BConv1] flux=%16.14f h=%16.14f Stdhdt=%16.14f diff=%e\n", tot_flux / ny, geometry->get_interface(), latent_est, fabs(tot_flux - latent_est));
			
			// geometry->set_flux(tot_flux / ny);
		}

		// update node positions after each newton step
		void actions_after_newton_step() {
			for (unsigned s = 0; s < bulk_mesh_pt->nspine(); s++) {
				bulk_mesh_pt->spine_pt(s)->height() = geometry->x1();
			}

			bulk_mesh_pt->node_update();
		}

		void set_initial_condition() {
			Vector<double> x(2);
			Vector<double> u(2);
			
			unsigned int tsteps = time_stepper_pt()->nprev_values();
			unsigned long int nnode = bulk_mesh_pt->nnode();
			unsigned int step = 0;

			double time = t_shift - tsteps*dt;
			time_pt()->time() = time;
			time_pt()->dt() = dt;
			double h = sqrt(_D[2] * time);

			geometry->set_interface(h);
			bulk_mesh_pt->node_update();

			for (unsigned int n = 0; n < nnode; n++) {
				bulk_mesh_pt->node_pt(n)->position(x);
				get_exact_u(time, x, u);
				bulk_mesh_pt->node_pt(n)->set_value(0, u[1]);
			}

			printf("[%2u] Setting initial confition at t=%8.6f\n", step, time);

			doc_step(step);
			step++;

			for (unsigned int t = 0; t < tsteps; t++) {
				shift_time_values();
				time += dt;
				time_pt()->time() = time;
				time_pt()->dt() = dt;
				h = sqrt(_D[2] * time);
				
				geometry->set_interface(h);
				bulk_mesh_pt->node_update();

				for (unsigned int n = 0; n < nnode; n++) {
					bulk_mesh_pt->node_pt(n)->position(x);
					get_exact_u(time, x, u);
					bulk_mesh_pt->node_pt(n)->set_value(0, u[1]);
				}

				printf("[%2u] Setting initial confition at t=%8.6f\n", step, time);
				doc_step(step);
				step++;
			}

			time_pt()->time() = time;
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
			fprintf(file, "%16.14f %16.14f %16.14f %16.14f %16.14f\n", time, tot_error, geometry->get_interface(), sqrt(_D[2]*time), max_elem_size/ys[1]);
			fclose(file);

			printf("[%4u] time=%8.6f error=%10.8e iface_err=%10.8e interface=%16.14f expected=%8.6f\n", timestep, time, tot_error, fabs(geometry->get_interface() - sqrt(_D[2] * time)), geometry->get_interface(), sqrt(_D[2]*time));
			info.number()++;
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
	CommandLineArgs::specify_command_line_flag("--nx2", &Nx2);
	CommandLineArgs::specify_command_line_flag("--ny", &Ny);
	CommandLineArgs::specify_command_line_flag("--tsteps", &t_steps);
	CommandLineArgs::specify_command_line_flag("--dt", &dt);
	CommandLineArgs::specify_command_line_flag("--tshift", &t_shift);
	CommandLineArgs::specify_command_line_flag("--vardt");
	CommandLineArgs::specify_command_line_flag("--write-freq", &write_freq);
	CommandLineArgs::specify_command_line_flag("--k1", &_k[0]);
	CommandLineArgs::specify_command_line_flag("--k2", &_k[1]);
	CommandLineArgs::specify_command_line_flag("--rho1", &_rho[0]);
	CommandLineArgs::specify_command_line_flag("--rho2", &_rho[1]);
	CommandLineArgs::specify_command_line_flag("--cp1", &_Cp[0]);
	CommandLineArgs::specify_command_line_flag("--cp2", &_Cp[1]);
	CommandLineArgs::specify_command_line_flag("--L", &L);
	CommandLineArgs::specify_command_line_flag("--De", &_D[2]);
	CommandLineArgs::specify_command_line_flag("--dname", &dname, "doc");

	CommandLineArgs::parse_and_assign();

	var_dt = CommandLineArgs::command_line_flag_has_been_set("--vardt");

	if (Nx1 == 0) {
		cout << "Error: Nx not specified" << endl;
		exit(1);
	}

	if (Nx2 == 0) {
		printf("nx2 not defined, assuming nx2 = nx1\n");
		Nx2 = Nx1;
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
	xs[1] = (_D[2] == 0.0) ? sqrt(1.26 * _interface_time) : sqrt(_D[2] * _interface_time);
	
	for (int i = 0; i < 2; i++) _D[i] = _k[i]/(_Cp[i]*_rho[i]);

	D = _D[1] / _D[0];
	k = _k[1] / _k[0];
	St = L / (_Cp[0] * (T_m - T_s));

	printf("Problem Def:\n");
	printf("\tk1=%8.6f k2=%8.6f\n", _k[0], _k[1]);
	printf("\trho1=%8.6f rho2=%8.6f\n", _rho[0], _rho[1]);
	printf("\tCp1=%8.6f Cp2=%8.6f\n", _Cp[0], _Cp[1]);
	printf("\tD1=%e D2=%e De=%e\n", _D[0], _D[1], _D[2]);
	printf("\tx0=%8.6f x1=%8.6f x2=%8.6f\n", xs[0], xs[1], xs[2]);
	printf("\tL=%8.6f D=%8.6f k=%8.6f St=%8.6f\n", L, D, k, St);

	// Nx2 = Nx1;
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
	fprintf(file, "k1=%10.8f\n", _k[0]);
	fprintf(file, "k2=%10.8f\n", _k[1]);
	fprintf(file, "rho1=%10.8f\n", _rho[0]);
	fprintf(file, "rho2=%10.8f\n", _rho[1]);
	fprintf(file, "cp1=%10.8f\n", _Cp[0]);
	fprintf(file, "cp2=%10.8f\n", _Cp[1]);
	fprintf(file, "L=%10.8f\n", L);
	fprintf(file, "De=%e\n", _D[2]);
	fclose(file);

	int prev_steps = problem.time_stepper_pt()->nprev_values()+1;

	for (uint t = 1; t < t_steps+1; t++) {
		problem.unsteady_newton_solve(dt);

		// if (t % write_freq == 0 || t == t_steps - 1)
		if (t % write_freq == 0)
			problem.doc_step(t+prev_steps-1);

		double x_int = geometry->get_interface();
		// printf("[%u] h = %16.14f\n", t, x_int);

		if (x_int >= xs[2] || isnan(x_int)) {
			printf("Interface reached right boundary, exiting...\n");
			break;
		}
	}
}