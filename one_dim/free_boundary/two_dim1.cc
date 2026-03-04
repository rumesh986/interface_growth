// moving interface
// calcualte interface position as part of solution in newton iterations
// move nodes with spine mesh
// use proper non-dimensionalised and corrected temperature scaling
// 		uses the form with separate alpha/beta constants
// interface temperature is not pinned, adds flux contributions for Stefan condition
// still accepts dimensinoal inputs and outputs non-dimensional values

#include <cmath>
#include <filesystem>

#include "generic.h"
#include "unsteady_heat.h"
// #include "unsteady_heat_flux_elements.h"

#include "two_phase_free_boundary_mesh.h"

using namespace oomph;

namespace params {
	struct material {
		double k = NAN;
		double rho = NAN;
		double cp = NAN;
		double D = NAN;
	};

	struct material solid;
	struct material liquid;

	double L = NAN;
	double De = NAN;

	double alpha = 0.0;
	double beta = 0.0;
	double D = 0.0;
	double St = 0.0;

	double T_s = NAN;
	double T_m = NAN;
	double T_l = NAN;

	double xs[3] = {0.0, 0.5, 1.0};
	double ys[2] = {0.0, 1.0};

	unsigned int nxs[2] = {10, 10};
	unsigned int nx = 0;
	unsigned int ny = 1;
	double dt = 0.01;
	double tstart = 0.5;
	double tend = 1.0;
	unsigned int write_freq = 1;
	std::string dname;

	unsigned int tsteps;

	void get_source(const double &t, const Vector<double> &x, double &source) {
		source = 0.0;
	}

	void no_flux_fct(const double &t, const Vector<double> &x, double &flux) {
		flux = 0.0;
	}

	void get_initial_interface_profile(const double &zeta, double &r) {
		r = 0.2 * sin(8.0 * MathematicalConstants::Pi * zeta) + 0.2;
		// r = 0.0;
	}

	void get_exact_u(const double &t, const Vector<double> &x, double &u) {
		double h_ana = sqrt(De * t);
		double eps = 0.0;
		get_initial_interface_profile(x[1], eps);

		double h = h_ana + eps;

		if (x[0] < h) {
			double denom = 0.5 / sqrt(t);
			u = (erf(x[0] * denom) / erf(h * denom)) - 1.0;
		} else {
			double denom = 0.5 / sqrt(D * t);
			double trans_tl = (T_l - T_m) / (T_m - T_s);
			u = trans_tl * (erf(x[0] * denom) - erf(h * denom)) / (1.0 - erf(h * denom));
		}
	}

	Vector<unsigned int> pinned_boundaries = {
		1,
		3
	};

	Vector<unsigned int> analytical_boundaries = {
		// 1
	};

	std::map<unsigned int, UnsteadyHeatEquations<X_ORDER>::UnsteadyHeatSourceFctPt> flux_boundaries = {

	};

}

template<class EL>
class TwoDimStefanProblem : public Problem {
	private:
		FreeBoundaryElement *geometry;
		TwoPhaseFreeBoundarySpineMesh<SpineElement<EL>> *bulk_mesh_pt;
		Mesh *surf_mesh_pt;
		Mesh *geom_mesh_pt;
		DocInfo info;

		unsigned int interface_boundary_index;
	
	public:
		TwoDimStefanProblem() {
			add_time_stepper_pt(new BDF<T_ORDER>(true));

			bulk_mesh_pt = new TwoPhaseFreeBoundarySpineMesh<SpineElement<EL>>(params::nxs, params::ny, params::xs, params::ys, time_stepper_pt());
			bulk_mesh_pt->setup_boundary_element_info();
			add_sub_mesh(bulk_mesh_pt);
			geometry = bulk_mesh_pt->geometry();

			geom_mesh_pt = new Mesh;
			geom_mesh_pt->add_element_pt(geometry);
			add_sub_mesh(geom_mesh_pt);

			interface_boundary_index = bulk_mesh_pt->free_boundary_index();

			surf_mesh_pt = new Mesh;
			for (unsigned int e = 0; e < bulk_mesh_pt->nboundary_element(interface_boundary_index); e++) {
				EL *elem = dynamic_cast<EL *>(bulk_mesh_pt->boundary_element_pt(interface_boundary_index, e));
				int face_index = bulk_mesh_pt->face_index_at_boundary(interface_boundary_index, e);
				if (face_index == 1) {
					auto flux_elem = new FreeBoundaryFluxElement<EL>(elem, face_index, params::St);

					// each flux elem now has mapping for all boundary nodes
					// need to remove the additional data if this works as expected
					for (unsigned int n = 0; n < bulk_mesh_pt->nboundary_node(interface_boundary_index); n++) {
						Node *node = bulk_mesh_pt->boundary_node_pt(interface_boundary_index, n);
						Data *data = geometry->get_data_for_node(node);
						flux_elem->add_node_data(node, data);
					}

					surf_mesh_pt->add_element_pt(flux_elem);
				}
			}
			add_sub_mesh(surf_mesh_pt);
			
			build_global_mesh();

			for (unsigned int e = 0; e < bulk_mesh_pt->nelement(); e++) {
				dynamic_cast<EL *>(bulk_mesh_pt->element_pt(e))->source_fct_pt() = params::get_source;
			}

			for (unsigned int b : params::pinned_boundaries) {
				for (unsigned int n = 0; n < bulk_mesh_pt->nboundary_node(b); n++) {
					bulk_mesh_pt->boundary_node_pt(b, n)->pin_all();
				}
			}

			for (unsigned int b : params::analytical_boundaries) {
				for (unsigned int n = 0; n < bulk_mesh_pt->nboundary_node(b); n++) {
					bulk_mesh_pt->boundary_node_pt(b, n)->pin_all();
				}
			}

			for (auto iter : params::flux_boundaries) {
				create_flux_elements(iter.first, iter.second);
			}

			for (unsigned int i = 0; i < params::ny; i++) {
				unsigned int base = i * params::nx;

				for (unsigned int e = params::nxs[0]; e < params::nx; e++) {
					EL *elem = dynamic_cast<EL *>(bulk_mesh_pt->element_pt(base + e));
					elem->alpha_pt() = &params::alpha;
					elem->beta_pt() = &params::beta;
				}
			}

			printf("Number of equations: %lu\n", assign_eqn_numbers());
			printf("NDOF: %lu\n", ndof());

			linear_solver_pt()->disable_doc_time();
			disable_info_in_newton_solve();
			newton_solver_tolerance() = 5e-10;
			max_newton_iterations() = 1e7;
			max_residuals() = 1e3;

			printf("Tolerance for newton solver: %e\n", newton_solver_tolerance());
		}

		~TwoDimStefanProblem() {
			delete bulk_mesh_pt;
			delete surf_mesh_pt;
			delete geom_mesh_pt;
		}

		void create_flux_elements(unsigned int b, UnsteadyHeatEquations<X_ORDER>::UnsteadyHeatSourceFctPt &flux_pt) {
			for (unsigned int e = 0; e < bulk_mesh_pt->nboundary_element(b); e++) {
				EL *elem = dynamic_cast<EL *>(bulk_mesh_pt->boundary_element_pt(b, e));
				int face_index = bulk_mesh_pt->face_index_at_boundary(b, e);
				UnsteadyHeatFluxElement<EL> *flux_elem = new UnsteadyHeatFluxElement<EL>(elem, face_index);
				flux_elem->flux_fct_pt() = flux_pt;
				bulk_mesh_pt->add_element_pt(flux_elem);
			}
		}

		void actions_before_implicit_timestep() {
			Vector<double> x(2);
			double u;
			double time = time_pt()->time();

			for (unsigned int b : params::analytical_boundaries) {
				for (unsigned int n = 0; n < bulk_mesh_pt->nboundary_node(b); n++) {
					bulk_mesh_pt->boundary_node_pt(b, n)->position(x);
					params::get_exact_u(time, x, u);
					bulk_mesh_pt->boundary_node_pt(b, n)->set_value(0, u);
				}
			}
		}

		void actions_before_newton_solve() {
			time_stepper_pt()->set_predictor_weights();
			geometry->calculate_predicted_values();

			for (unsigned int s = 0; s < bulk_mesh_pt->nspine(); s++) {
				double h_pred = geometry->geom_data_pt(s)->value(time_stepper_pt()->predictor_storage_index(), 0);
				bulk_mesh_pt->spine_pt(s)->height() = h_pred;
				// *geometry->geom_data_pt(s)->value_pt(0) = h_pred;
				geometry->geom_data_pt(s)->set_value(0, h_pred);
			}

			bulk_mesh_pt->node_update();
		}

		void actions_before_newton_step() {};
		void actions_before_newton_convergence_check() {}
		void actions_after_newton_step() {};
		void actions_after_newton_solve() {};
		void actions_after_implicit_timestep() {};

		void set_initial_condition() {
			Vector<double> x(2);
			double u;

			unsigned int tsteps = time_stepper_pt()->nprev_values();
			unsigned int nnode = bulk_mesh_pt->nnode();
			unsigned int step = 0;

			double time = params::tstart - tsteps * params::dt;
			time_pt()->time() = time;
			time_pt()->dt() = params::dt;

			double h = sqrt(params::De * time);
			double eps;
			for (unsigned int s = 0; s < bulk_mesh_pt->nspine(); s++) {
				Spine *spine = bulk_mesh_pt->spine_pt(s);
				params::get_initial_interface_profile(spine->geom_parameter(0), eps);
				spine->height() = h + eps;
				geometry->geom_data_pt(s)->set_value(0, h + eps);
			}
			bulk_mesh_pt->node_update();

			for (unsigned int n = 0; n < nnode; n++) {
				bulk_mesh_pt->node_pt(n)->position(x);
				params::get_exact_u(time, x, u);
				bulk_mesh_pt->node_pt(n)->set_value(0, u);
			}

			printf("[%2u] Setting initial condition at t=%8.6f\n", step, time);
			doc_step(step);
			step++;

			for (unsigned int t = 0; t < tsteps; t++) {
				shift_time_values();
				time += params::dt;

				time_pt()->time() = time;
				time_pt()->dt() = params::dt;

				double h = sqrt(params::De * time);
				for (unsigned int s = 0; s < bulk_mesh_pt->nspine(); s++) {
					Spine *spine = bulk_mesh_pt->spine_pt(s);
					params::get_initial_interface_profile(spine->geom_parameter(0), eps);
					spine->height() = h + eps;
					geometry->geom_data_pt(s)->set_value(0, h + eps);
				}
				bulk_mesh_pt->node_update();

				for (unsigned int n = 0; n < nnode; n++) {
					bulk_mesh_pt->node_pt(n)->position(x);
					params::get_exact_u(time, x, u);
					bulk_mesh_pt->node_pt(n)->set_value(0, u);
				}

				printf("[%2u] Setting initial condition at t=%8.6f\n", step, time);
				doc_step(step);
				step++;
			}

			for (unsigned int t = 0; t < time_stepper_pt()->nprev_values(); t++) {
				printf("h0 values: %u=%e\n", t, geometry->geom_data_pt(0)->value(t, 0));

			}

			time_pt()->time() = time;
		}

		double get_interface() {
			return 0.0;
		}

		void prepare_docs(DocInfo input) {
			info = input;


			char fname[256];
			sprintf(fname, "%s/results.dat", info.directory().c_str());
			FILE *file = fopen(fname, "w");
			fprintf(file, "time error h exact_h area\n");
			fclose(file);

			sprintf(fname, "%s/dim", info.directory().c_str());
			file = fopen(fname, "w");
			fprintf(file, "2");
			fclose(file);

			sprintf(fname, "%s/interface.dat", info.directory().c_str());
			file = fopen(fname, "w");
			fprintf(file, "time");
			for (unsigned int n = 0; n < bulk_mesh_pt->nboundary_node(interface_boundary_index); n++) {
				double y = bulk_mesh_pt->boundary_node_pt(interface_boundary_index, n)->position(1);
				fprintf(file, " %16.14f", y);
			}
			fprintf(file, "\n");
			fclose(file);
		}

		void doc_step(const unsigned int &timestep, const unsigned int &t = 0) {
			double time = time_pt()->time();
			unsigned long int nnode = bulk_mesh_pt->nnode();

			Vector<double> x(2), num_u(1);
			double ana_u;

			double tot_error = 0.0;

			char fname[256];
			sprintf(fname, "%s/steps/step%u.dat", info.directory().c_str(), info.number());
			FILE *file = fopen(fname, "w");
			fprintf(file, "x y exact_u u error\n");
			for (unsigned long int n = 0; n < nnode; n++) {
				bulk_mesh_pt->node_pt(n)->position(t, x);
				bulk_mesh_pt->node_pt(n)->value(t, num_u);
				params::get_exact_u(time, x, ana_u);

				double error = num_u[0] - ana_u;
				tot_error += error * error;

				fprintf(file, "%16.14f %16.14f %16.14f %16.14f %16.14f\n", x[0], x[1], ana_u, num_u[0], error);
			}
			fclose(file);

			tot_error = sqrt(tot_error) / nnode;

			double max_elem_size, min_elem_size;
			bulk_mesh_pt->max_and_min_element_size(max_elem_size, min_elem_size);

			double num_h = geometry->geom_data_pt(0)->value(0);

			sprintf(fname, "%s/results.dat", info.directory().c_str());
			file = fopen(fname, "a");
			fprintf(file, "%16.14f %16.14f %16.14f %16.14f %16.14f\n", time, tot_error, num_h, sqrt(params::De*time), max_elem_size/params::ys[1]);
			fclose(file);

			sprintf(fname, "%s/interface.dat", info.directory().c_str());
			file = fopen(fname, "a");
			fprintf(file, "%16.14f", time);
			for (unsigned int n = 0; n < bulk_mesh_pt->nboundary_node(interface_boundary_index); n++) {
				double x = bulk_mesh_pt->boundary_node_pt(interface_boundary_index, n)->position(0);
				fprintf(file, " %16.14f", x);
			}
			fprintf(file, "\n");
			fclose(file);

			printf("[%4u] time=%8.6f error=%10.8e iface_err=%10.8e interface=%16.14f expected=%8.6f\n", timestep, time, tot_error, fabs(num_h - sqrt(params::De * time)), num_h, sqrt(params::De * time));
			printf("interface positions: ");
			for (unsigned int i = 0; i < geometry->ngeom_data(); i++) printf("%16.14f ", geometry->geom_data_pt(i)->value(0));
			printf("\n");
			info.number()++;
		}
};

int main(int argc, char **argv) {
	CommandLineArgs::setup(argc, argv);

	CommandLineArgs::specify_command_line_flag("--k1", &params::solid.k);
	CommandLineArgs::specify_command_line_flag("--rho1", &params::solid.rho);
	CommandLineArgs::specify_command_line_flag("--cp1", &params::solid.cp);
	CommandLineArgs::specify_command_line_flag("--k2", &params::liquid.k);
	CommandLineArgs::specify_command_line_flag("--rho2", &params::liquid.rho);
	CommandLineArgs::specify_command_line_flag("--cp2", &params::liquid.cp);
	CommandLineArgs::specify_command_line_flag("--L", &params::L);
	CommandLineArgs::specify_command_line_flag("--De", &params::De);

	CommandLineArgs::specify_command_line_flag("--Ts", &params::T_s);
	CommandLineArgs::specify_command_line_flag("--Tm", &params::T_m);
	CommandLineArgs::specify_command_line_flag("--Tl", &params::T_l);

	CommandLineArgs::specify_command_line_flag("--nx1", &params::nxs[0]);
	CommandLineArgs::specify_command_line_flag("--nx2", &params::nxs[1]);
	CommandLineArgs::specify_command_line_flag("--ny", &params::ny);
	CommandLineArgs::specify_command_line_flag("--dt", &params::dt);
	CommandLineArgs::specify_command_line_flag("--tstart", &params::tstart);
	CommandLineArgs::specify_command_line_flag("--tend", &params::tend);
	CommandLineArgs::specify_command_line_flag("--write-freq", &params::write_freq);
	CommandLineArgs::specify_command_line_flag("--dname", &params::dname, "doc");

	CommandLineArgs::parse_and_assign();
	CommandLineArgs::output();

	params::tsteps = (unsigned int) ((params::tend - params::tstart) / params::dt);
	printf("tsteps set to %u\n", params::tsteps);

	if (isnan(params::solid.k) || isnan(params::solid.rho) || isnan(params::solid.cp) ||
		isnan(params::liquid.k) || isnan(params::liquid.rho) || isnan(params::liquid.cp) ||
		isnan(params::L) || isnan(params::De) || 
		isnan(params::T_s) || isnan(params::T_m) || isnan(params::T_l)
	) {
		
		printf("Properties not set correctly, exiting...\n\n");
		exit(1);
	}

	if (!CommandLineArgs::command_line_flag_has_been_set("--dname")) {
		char temp[256];
		sprintf(temp, "RESLT/%dn%u+%u_%dt%.2e", X_ORDER, params::nxs[0], params::nxs[1], T_ORDER, params::dt);
		params::dname.assign(temp);
	}

	printf("Saving results to %s\n", params::dname.c_str());
	if (std::filesystem::exists(params::dname.c_str())) {
		printf("Directory %s exists\n", params::dname.c_str());
	} else {
		std::filesystem::create_directories(params::dname.c_str());
	}

	char sub_dname[512];
	sprintf(sub_dname, "%s/steps", params::dname.c_str());
	if (std::filesystem::exists(sub_dname)) {
		printf("Directory %s exists\n", sub_dname);
	} else {
		std::filesystem::create_directories(sub_dname);
	}

	DocInfo info;
	info.set_directory(params::dname.c_str());
	info.number() = 0;

	printf("Output directory: %s\n", info.directory().c_str());

	params::alpha = (params::liquid.cp * params::liquid.rho) / (params::solid.cp * params::solid.rho);
	params::beta = params::liquid.k / params::solid.k;
	params::D = params::beta / params::alpha;
	params::St = params::L / (params::solid.cp * (params::T_m - params::T_s));
	params::nx = params::nxs[0] + params::nxs[1];
	params::xs[1] = sqrt(params::De * params::tstart);

	printf("Problem Def:\n");
	printf("\tsolid: k=%8.6f rho=%8.6f cp=%8.6f\n", params::solid.k, params::solid.rho, params::solid.cp);
	printf("\tliquid: k=%8.6f rho=%8.6f cp=%8.6f\n", params::liquid.k, params::liquid.rho, params::liquid.cp);
	printf("\tL=%8.6f De=%8.6f\n", params::L, params::De);
	printf("\tTs=%8.6f Tm=%8.6f Tl=%8.6f\n", params::T_s, params::T_m, params::T_l);
	printf("\talpha=%8.6f beta=%8.6f St=%8.6f\n", params::alpha, params::beta, params::St);
	printf("\tx0=%8.6f x1=%8.6f x2=%8.6f\n", params::xs[0], params::xs[1], params::xs[2]);

	auto problem = TwoDimStefanProblem<QUnsteadyHeatElement<2, X_ORDER>>();

	problem.initialise_dt(params::dt);
	problem.prepare_docs(info);
	problem.set_initial_condition();

	int prev_steps = problem.time_stepper_pt()->nprev_values()+1;
	for (unsigned int t = 0;t < params::tsteps; t++) {
		problem.unsteady_newton_solve(params::dt);

		if (t % params::write_freq == 0) problem.doc_step(t+prev_steps);

		double h = problem.get_interface();
		if (h >= params::xs[2] || isnan(h)) {
			printf("Interface passed out of domain, exiting...\n");
			break;
		}
	}

}