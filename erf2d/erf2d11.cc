
// moving interface
// calculate as part of solution in newton iterations (hopefully)
// move with spine meshes
// liquid phase pinned

#include <filesystem>

#include "includes.h"
// #include "two_phase_free_boundary_mesh.h"

#ifndef RUN_SCRIPT
#define X_ORDER 2
#define T_ORDER 1
#endif

double k[2] = {1.0, 0.25};
double rho[2] = {1.0, 2.0};
double Cp[2] = {1.0, 0.25};

double D[3] {k[0]/(Cp[0]*rho[0]), k[1]/(Cp[1]*rho[1]), 0.0};
double L = 1.0;

double e[2] = {0.0, 0.0};

bool ic_set = false;

static const double T_s = -1.0;
static const double T_m = 0.0;
static const double T_l = 1.0;
static const double T_0 = 0.0;

double xs[3] = {-1.0, 0.0, 1.0};
static const double ys[2] = {0.0, 1.0};


void no_flux_fct(const double &t, const Vector<double> &x, double &flux) {
	flux = 0.0;
}

void get_source(const double &t, const Vector<double> &x, double &source) {
	source = 0.0;
}

class TwoPhaseTestGeometry : public GeomObject {
	protected:
		Vector<Data *> data_pt;
		const unsigned int free_boundary_index = 2;

		bool destroy_geom_data = false;
		TimeStepper *ts_pt;
	
	public:
		TwoPhaseTestGeometry(
			const double &x0,
			const double &x1,
			const double &x2,
			const double &y0,
			const double &y1,
			TimeStepper *timestepper = new Steady<0>
		) : GeomObject(1,1), ts_pt(timestepper) {

			data_pt.resize(1);
			data_pt[0] = new Data(ts_pt, 5);

			for (unsigned int t = 0; t < ts_pt->nprev_values(); t++) {
				data_pt[0]->set_value(t, 0, x0);
				data_pt[0]->set_value(t, 1, x1);
				data_pt[0]->set_value(t, free_boundary_index, x2);
				data_pt[0]->set_value(t, 3, y0);
				data_pt[0]->set_value(t, 4, y1);
			}

			data_pt[0]->pin_all();
			destroy_geom_data = true;
		}

		~TwoPhaseTestGeometry() {
			if (destroy_geom_data)
				for (unsigned int i = 0; i < data_pt.size(); i++) delete data_pt[i];
		}

		double& x0(const unsigned int &t = 0) const {return *data_pt[0]->value_pt(t, 0);}
		double& x1(const unsigned int &t = 0) const {return *data_pt[0]->value_pt(t, 1);}
		double& x2(const unsigned int &t = 0) const {return *data_pt[0]->value_pt(t, free_boundary_index);}
		double& y0(const unsigned int &t = 0) const {return *data_pt[0]->value_pt(t, 3);}
		double& y1(const unsigned int &t = 0) const {return *data_pt[0]->value_pt(t, 4);}

		unsigned int ngeom_data() const {
			return data_pt.size();
		}

		Data *geom_data_pt(const unsigned int &j) {
			return data_pt[0];
		}

		double get_interface() {
			return x2();
		}

		void set_interface(double &x) {
			x2() = x;
		}

		void position(const unsigned int &t, const Vector<double> &zeta, Vector<double> &r) const {
			r[0] = x2(t);
		}

		void position(const Vector<double> &zeta, Vector<double> &r) const {
			position(0, zeta, r);
		}
};

class TwoPhaseTestElement : public GeneralisedElement,
							public TwoPhaseTestGeometry {
	private:
		double factor;
		unsigned int geometry_index;
		unsigned int flux_index;
		Data *flux_data_pt;

		TimeStepper *ts_pt;

	public:
		TwoPhaseTestElement(
			const double &x0,
			const double &x1,
			const double &x2,
			const double &y0,
			const double &y1,
			const double &factor,
			TimeStepper *timestepper = new Steady<0>
		) : TwoPhaseTestGeometry(x0, x1, x2, y0, y1, timestepper),
			factor(factor), ts_pt(timestepper) {

			geometry_index = add_internal_data(data_pt[0]);
			unpin_free_boundary();

			flux_data_pt = new Data(ts_pt, 1);
			flux_data_pt->set_value(0, 0.0);
			flux_data_pt->pin_all();

			flux_index = add_external_data(flux_data_pt);

			destroy_geom_data = false;
		}

		void pin_free_boundary() {
			internal_data_pt(geometry_index)->pin(free_boundary_index);
		}

		void unpin_free_boundary() {
			internal_data_pt(geometry_index)->unpin(free_boundary_index);
		}

		void set_flux(const unsigned int &t, const double &flux) {
			external_data_pt(flux_index)->set_value(t, 0, flux);
		}
		
		void set_flux(const double &flux) {
			set_flux(0, flux);
		}

		void get_residuals(Vector<double>& residuals) {
			residuals.initialise(0.0);
			DenseMatrix<double> placeholder(1);

			fill_in_generic_residual_contribution(residuals, placeholder, false);
		}

		void get_jacobian(Vector<double>& residuals, DenseMatrix<double>& jacobian) {
			residuals.initialise(0.0);
			jacobian.initialise(0.0);

			fill_in_generic_residual_contribution(residuals, jacobian, true);
		}
	
	protected:
	void fill_in_generic_residual_contribution(Vector<double>& residuals, DenseMatrix<double>& jacobian, bool compute_jacobian) {
			unsigned int ndofs = ndof();
			if (ndofs == 0) return;

			int free_boundary_local_eqn_number = internal_local_eqn(geometry_index, free_boundary_index);

			Data *interface_data_pt = internal_data_pt(geometry_index);
			TimeStepper *interface_ts_pt = interface_data_pt->time_stepper_pt();
			residuals[free_boundary_local_eqn_number] = factor * interface_ts_pt->time_derivative(1, interface_data_pt, free_boundary_index) - external_data_pt(flux_index)->value(0);

			if (compute_jacobian)
				jacobian(free_boundary_local_eqn_number, free_boundary_local_eqn_number) = factor * interface_ts_pt->weight(1, 0);
		}
};

template<class EL>
class TwoPhaseTestSpineMesh : public RectangularQuadMesh<EL>,
								public SpineMesh {
	private:
		unsigned int nx1, nx2, nx, ny;
		TwoPhaseTestGeometry *geometry;
		TimeStepper *ts_pt;

	public:
		TwoPhaseTestSpineMesh(
			const unsigned int &nx1,
			const unsigned int &nx2,
			const unsigned int &ny,
			TwoPhaseTestGeometry *geometry,
			TimeStepper *timestepper = &Mesh::Default_TimeStepper
		) : RectangularQuadMesh<EL>(nx1+nx2, ny, geometry->x0(), geometry->x2(), geometry->y0(), geometry->y1(), timestepper),
			nx1(nx1), nx2(nx2), nx(nx1+nx2), ny(ny), geometry(geometry), ts_pt(timestepper) {

			this->set_nboundary(5);
			for (unsigned int e = 0; e < ny; e++) {
				EL *elem = dynamic_cast<EL *>(this->element_pt(nx1*(1+e) + nx2*e));
				unsigned int nnode = elem->nnode_1d();

				for (unsigned int n = 0; n < nnode; n++) {
					Node *node = elem->node_pt(nnode * n);
					this->convert_to_boundary_node(node);
					this->add_boundary_node(4, node);
				}
			}

			construct_spines();

			for (unsigned int n = 0; n < nnode(); n++) spine_node_update(node_pt(n));

			this->setup_boundary_element_info();
		}

		void construct_spines() {
			unsigned int np = finite_element_pt(0)->nnode_1d();
			unsigned int nspine = (np - 1) * ny + 1;
			Spine_pt.reserve(nspine);

			unsigned int yi = 0;
			for (yi = 0; yi < ny-1; yi++) {
				for (unsigned int s = 0; s < np-1; s++) {
					Spine *spine = new Spine(geometry->x2());
					spine->spine_height_pt()->pin_all();
					Spine_pt.push_back(spine);

					// irrelevant for now
					Vector<double> parameters = {((double)yi + (double)s/(double)(np-1)) / (double)ny};
					spine->set_geom_parameter(parameters);

					Vector<GeomObject *> geom_object_pt = {geometry};
					spine->set_geom_object_pt(geom_object_pt);

					unsigned int xi = 0;
					for (xi = 0; xi < nx1-1; xi++) {
						for (unsigned int n = 0; n < np-1; n++) {
							SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
							node_pt->spine_pt() = spine;
							node_pt->fraction() = ((double)xi + (double)n/(double)(np-1)) / (double)nx1;
							node_pt->spine_mesh_pt() = this;
							node_pt->node_update_fct_id() = 0;
						}
					}

					xi = nx1-1;
					for (unsigned int n = 0; n < np; n++) {
						SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
						node_pt->spine_pt() = spine;
						node_pt->fraction() = ((double)xi + (double)n/(double)(np-1)) / (double)nx1;
						node_pt->spine_mesh_pt() = this;
						node_pt->node_update_fct_id() = 0;
					}

					// phase 2
					for (xi = nx1; xi < nx-1; xi++) {
						for (unsigned int n = 0; n < np-1; n++) {
							SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
							node_pt->spine_pt() = spine;
							node_pt->fraction() = ((double)(xi-nx1) + (double)n/(double)(np-1)) / (double)nx2;
							node_pt->spine_mesh_pt() = this;
							node_pt->node_update_fct_id() = 1;
						}
					}

					xi = nx-1;
					for (unsigned int n = 0; n < np; n++) {
						SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
						node_pt->spine_pt() = spine;
						node_pt->fraction() = ((double)(xi-nx1) + (double)n/(double)(np-1)) / (double)nx2;
						node_pt->spine_mesh_pt() = this;
						node_pt->node_update_fct_id() = 1;
					}
				}
			}

			yi = ny-1;
			for (unsigned int s = 0; s < np; s++) {
				Spine *spine = new Spine(geometry->x2());
				spine->spine_height_pt()->pin_all();
				Spine_pt.push_back(spine);

				// irrelevant for now
				Vector<double> parameters = {((double)yi + (double)s/(double)(np-1)) / (double)ny};
				spine->set_geom_parameter(parameters);

				Vector<GeomObject *> geom_object_pt = {geometry};
				spine->set_geom_object_pt(geom_object_pt);

				unsigned int xi = 0;
				for (xi = 0; xi < nx1-1; xi++) {
					for (unsigned int n = 0; n < np-1; n++) {
						SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
						node_pt->spine_pt() = spine;
						node_pt->fraction() = ((double)xi + (double)n/(double)(np-1)) / (double)nx1;
						node_pt->spine_mesh_pt() = this;
						node_pt->node_update_fct_id() = 0;
					}
				}

				xi = nx1-1;
				for (unsigned int n = 0; n < np; n++) {
					SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
					node_pt->spine_pt() = spine;
					node_pt->fraction() = ((double)xi + (double)n/(double)(np-1)) / (double)nx1;
					node_pt->spine_mesh_pt() = this;
					node_pt->node_update_fct_id() = 0;
				}

				// phase 2
				for (xi = nx1; xi < nx-1; xi++) {
					for (unsigned int n = 0; n < np-1; n++) {
						SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
						node_pt->spine_pt() = spine;
						node_pt->fraction() = ((double)(xi-nx1) + (double)n/(double)(np-1)) / (double)nx2;
						node_pt->spine_mesh_pt() = this;
						node_pt->node_update_fct_id() = 1;
					}
				}

				xi = nx-1;
				for (unsigned int n = 0; n < np; n++) {
					SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
					node_pt->spine_pt() = spine;
					node_pt->fraction() = ((double)(xi-nx1) + (double)n/(double)(np-1)) / (double)nx2;
					node_pt->spine_mesh_pt() = this;
					node_pt->node_update_fct_id() = 1;
				}
			}

			printf("Created %lu/%u spines\n\n", Spine_pt.size(), nspine);
		}

		void spine_node_update(SpineNode *node_pt) {
			TwoPhaseTestGeometry *geom = dynamic_cast<TwoPhaseTestGeometry *>(node_pt->spine_pt()->geom_object_pt(0));
			double frac = node_pt->fraction();
			Vector<double> zeta(1), r(1);
			zeta[0] = 0.0;
			geom->position(zeta, r);

			switch(node_pt->node_update_fct_id()) {
				case 0: node_pt->x(0) = geom->x0() + frac * (geom->x1() - geom->x0());	break;
				case 1: node_pt->x(0) = geom->x1() + frac * (r[0] - geom->x1());		break;
				default:	printf("Invalid node update function\n");					break;
			}
		}
};

TwoPhaseTestElement *geometry;

// analytical solution
// u[0] -> solution with interface in current position (from numerical solution)
// u[1] -> solution with interface position from analytical solution
void get_exact_u(const double &t, const Vector<double> &x, Vector<double> &u) {
	// double h = 0.0;
	// if (ic_set) {
	// 	h = (D[2] == 0.0) ? geometry->get_interface() : sqrt(D[2] * t);
	// } else {
	// 	h = geometry->get_interface();
	// }

	if (x[0] < 0.0) {
		u[0] = T_0 + (T_0 - T_s) * erf(x[0]/(2*sqrt(D[0] * t)));
	} else {
		u[0] = T_0 + (T_0 - T_s) * erf(x[0]/(2*sqrt(D[1] * t))) * (e[0]/e[1]);
	}

	u[1] = u[0];
}

Vector<unsigned int> analytical_boundaries = {
	1, // right
	3, // left
	4 // interface
};

Vector<unsigned int> pinned_boundaries = {
	// 1, // right
	// 3, // left
	// 4 // interface
};

map<unsigned int, FluxFctPt> flux_boundaries = {
	// {0, no_flux_fct},
	// {2, no_flux_fct},
};

template<class EL>
class Erf2DProblem : public Problem {
	private:
		unsigned int nx1, nx2, nx, ny;
		unsigned int t_steps;
		double dt, t_shift;

		DocInfo info;
	
		TwoPhaseTestSpineMesh<SpineElement<EL>> *bulk_mesh_pt;
		Mesh *geometry_mesh_pt;

		public:
		Erf2DProblem(
			unsigned int nx_,
			unsigned int ny_,
			unsigned int tsteps,
			double dt_,
			double tshift,
			DocInfo info_
		) : nx1(nx_), nx2(nx_), nx(nx_+nx_), ny(ny_), t_steps(tsteps), dt(dt_), t_shift(tshift), info(info_) {

			add_time_stepper_pt(new BDF<T_ORDER>);

			geometry = new TwoPhaseTestElement(xs[0], xs[1], xs[2], ys[0], ys[1], rho[1] * L, time_stepper_pt());

			printf("x0=%8.6f x1=%8.6f x2=%8.6f\n", geometry->x0(), geometry->x1(), geometry->x2());

			bulk_mesh_pt = new TwoPhaseTestSpineMesh<SpineElement<EL>>(nx1, nx2, ny, geometry, time_stepper_pt());
			bulk_mesh_pt->setup_boundary_element_info();
			add_sub_mesh(bulk_mesh_pt);

			geometry_mesh_pt = new Mesh;
			geometry_mesh_pt->add_element_pt(geometry);
			add_sub_mesh(geometry_mesh_pt);

			build_global_mesh();

			for (unsigned int e = 0; e < bulk_mesh_pt->nelement(); e++)
				dynamic_cast<EL *>(bulk_mesh_pt->element_pt(e))->source_fct_pt() = get_source;
			
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
				for (unsigned int e = 0; e < nx1; e++)
					dynamic_cast<EL *>(bulk_mesh_pt->element_pt(base + e))->beta_pt() = &D[0];
				for (unsigned int e = nx1; e < nx; e++) {
					dynamic_cast<EL *>(bulk_mesh_pt->element_pt(base + e))->beta_pt() = &D[1];

					// EL *elem = dynamic_cast<EL *>(bulk_mesh_pt->element_pt(base + e));
					// unsigned int nnode = elem->nnode();
					// for (unsigned int n = 0; n < nnode; n++) {
					// 	elem->node_pt(n)->pin_all();
					// }
				}
			}

			printf("Total number of equations: %lu\n", assign_eqn_numbers());
			printf("NDOF: %lu\n", ndof());

			// for (unsigned int e = 0; e < nx*ny; e++) {
			// 	printf("[%u] ngeom=%u\n", e, dynamic_cast<SpineElement<EL>*>(bulk_mesh_pt->element_pt(e))->ngeom_data());
			// }

			linear_solver_pt()->disable_doc_time();
			disable_info_in_newton_solve();
			// newton_solver_tolerance() = 1e-9;
			max_newton_iterations() = 1e5;
			max_residuals() = 1e3;
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

		void actions_before_newton_solve() {};
		void actions_after_newton_solve() {};

		void actions_before_implicit_timestep() {}
		void actions_after_implicit_timestep() {};

		// void actions_before_newton_convergence_check() {
		void actions_before_newton_step() {
			Vector<double> flux(2);

			double tot_flux = 0.0;
			unsigned long int nelems = bulk_mesh_pt->nboundary_element(1);
			for (unsigned long int e = 0; e < nelems; e++) {
				int face_index = bulk_mesh_pt->face_index_at_boundary(1, e);
				if (face_index == 1) {
					Vector<double> s = {1.0, 0.0};

					EL *elem = dynamic_cast<EL *>(bulk_mesh_pt->boundary_element_pt(1, e));
					elem->get_flux(s, flux);

					Vector<double> x(2);
					elem->node_pt(0)->position(x);

					// printf("e=%u x0=%6.4f x1=%6.4f flux0=%8.6f\n", e, x[0], x[1], flux[0]);

					tot_flux += k[1] * flux[0];
				}
			}

			printf("Setting total flux to %8.6f interface at %8.6f\n", tot_flux / ny, geometry->x2());
			geometry->set_flux(tot_flux / (double) ny);

			
			for (unsigned s = 0; s < bulk_mesh_pt->nspine(); s++) {
				bulk_mesh_pt->spine_pt(s)->height() = geometry->x2();
			}

			bulk_mesh_pt->node_update();

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

		void set_initial_condition() {
			time_pt()->time() = t_shift;
			unsigned long int nnode = bulk_mesh_pt->nnode();

			Vector<double> x(2);
			Vector<double> u(1);

			for (unsigned long int n = 0; n < nnode; n++)
				time_stepper_pt()->assign_initial_positions_impulsive(bulk_mesh_pt->node_pt(n));
			
			unsigned int tsteps = time_stepper_pt()->nprev_values();

			for (unsigned long int n = 0; n < nnode; n++) {
				bulk_mesh_pt->node_pt(n)->position(tsteps, x);
				get_exact_u(time_pt()->time(tsteps), x, u);
				bulk_mesh_pt->node_pt(n)->set_value(tsteps, 0, u[0]);
			}

			unsigned int step = 0;
			doc_step(step, tsteps);
			step++;

			for (int t = tsteps -1; t >= 0; t--) {
				double time = time_pt()->time((unsigned int) t);

				//
				// update_interface
				//
				double new_h = sqrt(D[2] * time);

				geometry->set_interface(new_h);
				printf("IC %u: x_old=%8.6f x_new=%8.6f\n", t, geometry->get_interface(), new_h);
				bulk_mesh_pt->node_update();

				for (unsigned long int n = 0; n < nnode; n++) {
					bulk_mesh_pt->node_pt(n)->position(t, x);
					get_exact_u(time, x, u);
					bulk_mesh_pt->node_pt(n)->set_value(t, 0, u[0]);
				}

				printf("[% 4d] Setting initial condition at t=%8.6f\n", step, time);
				doc_step(step, t);
				step++;
			}

			time_pt()->time() = t_shift;
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
			fprintf(file, "%16.14f %16.14f %16.14f %16.14f %16.14f\n", time, tot_error, geometry->get_interface(), sqrt(D[2]*time), max_elem_size);
			fclose(file);

			printf("[%4u] time=%8.6f error = %e interface = %8.6f expected = %8.6f\n", timestep, time, tot_error, geometry->get_interface(), sqrt(D[2]*time));
			info.number()++;
		}

		// modified from UnsteadyHeatEquations::get_flux from unsteady_heat_elements.h
		void get_flux_ic(const unsigned int &t, EL * elem, const Vector<double> &s, Vector<double> &flux) const {
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
			for (unsigned j = 0; j < 2; j++) {
				flux[j] = 0.0;
			}

			// Loop over nodes
			for (unsigned l = 0; l < n_node; l++) {
				// Loop over derivative directions
				for (unsigned j = 0; j < 2; j++) {
					flux[j] += elem->nodal_value(t, l, u_nodal_index) * dpsidx(l, j);
				}
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
	std::string dname;

	CommandLineArgs::specify_command_line_flag("--nx", &Nx);
	CommandLineArgs::specify_command_line_flag("--ny", &Ny);
	CommandLineArgs::specify_command_line_flag("--tsteps", &t_steps);
	CommandLineArgs::specify_command_line_flag("--dt", &dt);
	CommandLineArgs::specify_command_line_flag("--tshift", &t_shift);
	CommandLineArgs::specify_command_line_flag("--vardt");
	CommandLineArgs::specify_command_line_flag("--write-freq", &write_freq);
	CommandLineArgs::specify_command_line_flag("--k1", &k[0]);
	CommandLineArgs::specify_command_line_flag("--k2", &k[1]);
	CommandLineArgs::specify_command_line_flag("--rho1", &rho[0]);
	CommandLineArgs::specify_command_line_flag("--rho2", &rho[1]);
	CommandLineArgs::specify_command_line_flag("--cp1", &Cp[0]);
	CommandLineArgs::specify_command_line_flag("--cp2", &Cp[1]);
	CommandLineArgs::specify_command_line_flag("--L", &L);
	CommandLineArgs::specify_command_line_flag("--De", &D[2]);
	CommandLineArgs::specify_command_line_flag("--dname", &dname, "doc");

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

	if (!CommandLineArgs::command_line_flag_has_been_set("--dname")) {
		char temp[256];
		sprintf(temp, "RESLT/%dn%u_%dt%.2e", X_ORDER, Nx, T_ORDER, dt);
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
	xs[2] = (D[2] == 0.0) ? sqrt(1.26 * _interface_time) : sqrt(D[2] * _interface_time);
	
	for (int i = 0; i < 2; i++) {
		D[i] = k[i]/(Cp[i]*rho[i]);
		e[i] = sqrt(k[i] * rho[i] * Cp[i]);
	}

	printf("Problem Def:\n");
	printf("\tk1=%8.6f k2=%8.6f\n", k[0], k[1]);
	printf("\trho1=%8.6f rho2=%8.6f\n", rho[0], rho[1]);
	printf("\tCp1=%8.6f Cp2=%8.6f\n", Cp[0], Cp[1]);
	printf("\tD1=%8.6f D2=%8.6f De=%8.6f\n", D[0], D[1], D[2]);
	printf("\te1=%8.6f e2=%8.6f\n", e[0], e[1]);
	printf("\tL=%8.6f\n", L);

	auto problem = Erf2DProblem<QUnsteadyHeatElement<2,X_ORDER>>(Nx, Ny, t_steps, dt, t_shift, info);

	problem.initialise_dt(dt);
	problem.set_initial_condition();

	// save current parameters to file for future reference if needed
	char config_fname[256];
	sprintf(config_fname, "%s/config", dname.c_str());

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
	fprintf(file, "k1=%10.8f\n", k[0]);
	fprintf(file, "k2=%10.8f\n", k[1]);
	fprintf(file, "rho1=%10.8f\n", rho[0]);
	fprintf(file, "rho2=%10.8f\n", rho[1]);
	fprintf(file, "cp1=%10.8f\n", Cp[0]);
	fprintf(file, "cp2=%10.8f\n", Cp[1]);
	fprintf(file, "L=%10.8f\n", L);
	fclose(file);

	int prev_steps = problem.time_stepper_pt()->nprev_values()+1;

	for (uint t = 0; t < t_steps; t++) {
		problem.unsteady_newton_solve(dt);

		if (t % write_freq == 0)
			problem.doc_step(t+prev_steps);

		double x_int = geometry->get_interface();

		if (isnan(x_int)) {
			printf("Interface reached right boundary, exiting...\n");
			break;
		}
	}
}