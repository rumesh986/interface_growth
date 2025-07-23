#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "includes.h"
#include "erf2d1.h"
#include "erf2d2.h"
#include "erf2d3.h"
#include "two_layer_mesh.h"

#ifndef RUN_SCRIPT
#define X_ORDER 2
#define T_ORDER 1
#endif

template<class EL, class P> class Solver : public Problem {
	private:
		P problem;
		DocInfo info;

	public:
		Solver() {
			cout << "Error: constructor not specified, exiting ..." << endl;
			exit(1);
		}

		Solver(Erf2D1Problem problem, DocInfo info) : problem(problem), info(info) {
			add_time_stepper_pt(new BDF<T_ORDER>);

			mesh_pt() = new RectangularQuadMesh<EL>(problem.Nx, problem.Ny, problem.Lx, problem.Ly, time_stepper_pt());
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

			for (int b = 0; b < mesh_pt()->nboundary(); b++) {
				for  (int n = 0; n < mesh_pt()->nboundary_node(b); n++) {
					Node *node = mesh_pt()->boundary_node_pt(b, n);
					printf("b %d n %d x %10.8f y %10.8f\n", b, n, node->x(0), node->x(1));
				}
			}

			cout << "Number of equations " << assign_eqn_numbers() << endl;
		}

		Solver(Erf2D2Problem problem, DocInfo info) : problem(problem), info(info) {
			add_time_stepper_pt(new BDF<T_ORDER>);

			mesh_pt() = new TwoLayer2DMesh<EL>(problem.Nx1, problem.Nx2, problem.Ny, -problem.Lx, 0.0, problem.Lx, 0.0, problem.Ly, time_stepper_pt());
			mesh_pt()->setup_boundary_element_info();

			printf("Initial count of elements: %u\n", mesh_pt()->nelement());

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

			for (int b = 0; b < mesh_pt()->nboundary(); b++) {
				for  (int n = 0; n < mesh_pt()->nboundary_node(b); n++) {
					printf("before\n");
					Node *node = mesh_pt()->boundary_node_pt(b, n);
					printf("b %d n %d x %10.8f y %10.8f\n", b, n, node->x(0), node->x(1));
				}
			}

			// beta_pt() expects a double*, but we define it as static const
			for (uint e = 0; e < problem.Nx1; e++) {
				dynamic_cast<EL *>(mesh_pt()->element_pt(e))->beta_pt() = (double *) &problem.kappa1;
			}

			for (uint e = problem.Nx1; e < problem.Nx; e++) {
				dynamic_cast<EL *>(mesh_pt()->element_pt(e))->beta_pt() = (double *) &problem.kappa2;
			}

			for (uint e = 0; e < problem.Nx; e++)
				printf("e = %u; kappa = %f\n", e, * dynamic_cast<EL *>(mesh_pt()->element_pt(e))->beta_pt());

			cout << "Number of equations " << assign_eqn_numbers() << endl;
		}

		Solver(Erf2D3Problem problem, DocInfo info) : problem(problem), info(info) {
			add_time_stepper_pt(new BDF<T_ORDER>);

			TwoPhaseDomain *domain = new TwoPhaseDomain(-problem.Lx, 0.0, problem.Lx, problem.vel, problem.Nx1, problem.Nx2, problem.Ny, time_pt());

			mesh_pt() = new RefineableTwoLayer2DMesh<EL>(problem.Nx1, problem.Nx2, problem.Ny, -problem.Lx, 0.0, problem.Lx, 0.0, problem.Ly, domain, time_stepper_pt());
			mesh_pt()->setup_boundary_element_info();

			printf("Initial count of elements: %u\n", mesh_pt()->nelement());

			for (uint e = 0; e < mesh_pt()->nelement(); e++)
				dynamic_cast<EL *>(mesh_pt()->element_pt(e))->source_fct_pt() = problem.get_source;

			for (auto iter = problem.flux_boundaries.begin(); iter != problem.flux_boundaries.end(); iter++) {
				cout << "Creating flux elements at boundary " << iter->first << endl;
				create_flux_elements(iter->first, iter->second);
			}

			for (auto iter = problem.pinned_boundaries.begin(); iter != problem.pinned_boundaries.end(); iter++) {
				cout << "Assigning pinned boundary condition at " << iter->first << endl;
				uint nnode = mesh_pt()->nboundary_node(iter->first);
				for (uint n = 0; n < nnode; n++) {
					mesh_pt()->boundary_node_pt(iter->first, n)->set_value(0, iter->second);
					mesh_pt()->boundary_node_pt(iter->first, n)->pin(0);
				}
			}

			for (uint i = 0; i < problem.analytical_boundaries.size(); i++) {
				uint b = problem.analytical_boundaries[i];
				for (uint n = 0; n < mesh_pt()->nboundary_node(b); n++) {
					mesh_pt()->boundary_node_pt(b, n)->pin(0);
				}
			}

			// for (uint b = 0; b < mesh_pt()->nboundary(); b++) {
			// 	for  (uint n = 0; n < mesh_pt()->nboundary_node(b); n++) {
			// 		Node *node = mesh_pt()->boundary_node_pt(b, n);
			// 		// printf("b %d n %d x %10.8f y %10.8f\n", b, n, node->x(0), node->x(1));
			// 	}
			// }

			// beta_pt() expects a double*, but we define it as static const

			for (uint yi = 0; yi < problem.Ny; yi++) {
				for (uint e = 0; e < problem.Nx1; e++) {
					EL *elem = dynamic_cast<EL *>(mesh_pt()->element_pt(yi * problem.Nx + e));
					elem->beta_pt() = (double *) &problem.kappa1;
					// elem->set_macro_elem_pt(domain->macro_element_pt(0));
				}
	
				for (uint e = problem.Nx1; e < problem.Nx; e++) {
					EL *elem = dynamic_cast<EL *>(mesh_pt()->element_pt(yi * problem.Nx + e));
					elem->beta_pt() = (double *) &problem.kappa2;
					// elem->set_macro_elem_pt(domain->macro_element_pt(1));
				}
			}

			for (uint e = 0; e < problem.Nx * problem.Ny; e++)
				printf("e = %u; kappa = %f\n", e, * dynamic_cast<EL *>(mesh_pt()->element_pt(e))->beta_pt());


			// mesh_pt()->node_update();

			// domain->doc_domain(info, "");

			// Vector<double> s;
			// Vector<double> r;

			// for (uint e = 0; e < mesh_pt()->nelement(); e++) {
			// 	FiniteElement *elem = dynamic_cast<FiniteElement *>(mesh_pt()->element_pt(e));
			// 	if (elem->macro_elem_pt() == 0) {
			// 		printf("macro not set for element %d\n", e);
			// 		continue;
			// 	}

			// 	s.resize(elem->dim());
			// 	// r.resize(2);

			// 	uint nnode = elem->nnode();
			// 	for (uint j = 0; j < nnode; j++) {
			// 		Node *node = elem->node_pt(j);

			// 		r.resize(node->ndim());

			// 		elem->local_coordinate_of_node(j,s);
			// 		// cout << "(s: ";
			// 		printf("(%d, %u)  s: ", e, j);
			// 		for (auto i: s)	cout << i << " ";

			// 		elem->get_x(s, r);
			// 		cout << "\tr: ";
			// 		// printf("(%d, %u) r: ", e, j);
			// 		for (auto i: r)	cout << i << " ";
			// 		cout << endl;
			// 	}
			// }

			cout << "Number of equations " << assign_eqn_numbers() << endl;
		}

		~Solver() {
			delete mesh_pt();
		}

		void actions_before_newton_solve() {};
		void actions_after_newton_solve() {};

		void actions_before_implicit_timestep() {
			if (typeid(P) == typeid(Erf2D3Problem)) {
				mesh_pt()->node_update();
				printf("b4 loc: x=%10.8f\n", mesh_pt()->boundary_node_pt(4, 0)->x(0));
			}

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

		void actions_after_implicit_timestep() {};

		void create_flux_elements(uint b, FluxFctPt &flux_ptr) {
			uint n_elems = mesh_pt()->nboundary_element(b);

			for (uint e = 0; e < n_elems; e++) {
				EL *elem = dynamic_cast<EL *>(mesh_pt()->boundary_element_pt(b, e));
				int face_index = mesh_pt()->face_index_at_boundary(b, e);

				UnsteadyHeatFluxElement<EL> *flux_elem = new UnsteadyHeatFluxElement<EL>(elem, face_index);
				flux_elem->flux_fct_pt() = flux_ptr;
				mesh_pt()->add_element_pt(flux_elem);

				printf("at b=%d, e=%d, f_i=%d\n", b, e, face_index);
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
					x[1] = mesh_pt()->node_pt(n)->x(1);

					problem.get_exact_u(cur_t, x, u);
					mesh_pt()->node_pt(n)->set_value(t, 0, u[0]);
					printf("[%10.8f] x=%10.8f, y=%10.8f, u=%10.8f\n", cur_t, x[0], x[1], u[0]);
				}
			}

			// if (typeid(P) == typeid(Erf2D3Problem)) {
			// 	printf("Updating mesh in set_initial_conditions\n");
			// 	Vector<double> s;
			// 	Vector<double> r;

			// 	for (uint e = 0; e < mesh_pt()->nelement(); e++) {
			// 		EL *elem = dynamic_cast<EL *>(mesh_pt()->element_pt(e));
			// 		if (elem->macro_elem_pt() == 0) printf("macro not set for element %d\n", e);
					
			// 		// if (elem == 0)	printf("elem %d is 0\n", e);
			// 		// else 			printf("elem %d is not 0\n", e);

			// 		// if (elem->macro_elem_pt() != 0) printf("Got Macro\n");
					
			// 		s.resize(elem->dim());
			// 		r.resize(2);

			// 		elem->local_coordinate_of_node(0, s);
			// 		printf("Got local coordinates for elem %d\n", e);
			// 		printf("s: %10.8f %10.8f\n", s[0], s[1]);

			// 		elem->get_x(s, r);
			// 		printf("Got x for elem %d\n", e);
			// 		printf("r: %10.8f %10.8f\n", r[0], r[1]);
			// 	}

				// mesh_pt()->node_update();
				// printf("b4 loc: x=%10.8f\n", mesh_pt()->boundary_node_pt(4, 0)->x(0));
			// }

			time_pt()->time() = problem.t_shift;
		};

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

#endif