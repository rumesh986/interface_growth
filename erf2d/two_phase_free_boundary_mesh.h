#ifndef __TWO_PHASE_FREE_BOUNDARY_MESH_H__
#define __TWO_PHASE_FREE_BOUNDARY_MESH_H__

#include "generic.h"
#include "meshes/rectangular_quadmesh.h"

#include "includes.h"

// geometry object for (pseudo 1D) phase change problem
// Just stores boundary locations, including interface
class FreeBoundaryGeometry : public GeomObject {
	protected:
		Vector<Data *> data_pt;
		const unsigned int free_boundary_index = 1;

		bool destroy_geom_data = false;
	public:
		FreeBoundaryGeometry(
			const double x0,
			const double x1,
			const double x2,
			const double y0,
			const double y1,
			TimeStepper *timestepper = new Steady<0>
		) : GeomObject(2,2, timestepper) {
			
			data_pt.resize(1);
			data_pt[0] = new Data(time_stepper_pt(), 5);

			// assume impulsive conditions
			for (unsigned int t = 0; t < time_stepper_pt()->nprev_values(); t++) {
				data_pt[0]->set_value(t, 0, x0);
				data_pt[0]->set_value(t, free_boundary_index, x1);
				data_pt[0]->set_value(t, 2, x2);
				data_pt[0]->set_value(t, 3, y0);
				data_pt[0]->set_value(t, 4, y1);
			}

			data_pt[0]->pin_all();

			destroy_geom_data = true;
		}

		~FreeBoundaryGeometry() {
			if (destroy_geom_data) {
				for (unsigned int i = 0; i < data_pt.size(); i++) {
					delete data_pt[i];
				}
			}
		}

		double& x0(const unsigned int &t = 0) const {return *data_pt[0]->value_pt(t, 0);}
		double& x1(const unsigned int &t = 0) const {return *data_pt[0]->value_pt(t, free_boundary_index);}
		double& x2(const unsigned int &t = 0) const {return *data_pt[0]->value_pt(t, 2);}
		double& y0(const unsigned int &t = 0) const {return *data_pt[0]->value_pt(t, 3);}
		double& y1(const unsigned int &t = 0) const {return *data_pt[0]->value_pt(t, 4);}

		unsigned int ngeom_data() const {return data_pt.size();}

		Data* geom_data_pt(const unsigned &j) {return data_pt[0];}
		
		double get_interface(const unsigned int &t) {return x1(t);}
		double get_interface() {return get_interface(0);}

		void set_interface(const unsigned int &t, double &x) {x1(t) = x;}
		void set_interface(double &x) {set_interface(0, x);}

		void position(const unsigned int &t, const Vector<double> &zeta, Vector<double> &r) const {
			r[0] = x1(t);
		}

		void position(const Vector<double> &zeta, Vector<double> &r) const {
			position(0, zeta, r);
		}
};

// Combine FreeBoundaryGeometry with a GeneralisedElement
// allows us to add an equation to the main FEM solver
class FreeBoundaryElement : public GeneralisedElement, 
							public FreeBoundaryGeometry {

	friend QUnsteadyHeatElement<2, X_ORDER>;
	private:
		double St, k;
		unsigned int geometry_index;
		unsigned int flux_index;
		Data *flux_data_pt;
		Vector<QUnsteadyHeatElement<2, X_ORDER> *> phase1Elements = {};
		Vector<QUnsteadyHeatElement<2, X_ORDER> *> phase2Elements = {};
		Vector<unsigned int> phase1_indexes;
		Vector<unsigned int> phase2_indexes;

		TimeStepper *ts_pt;

	public:
		FreeBoundaryElement(
			const double x0,
			const double x1,
			const double x2,
			const double y0,
			const double y1,
			const double factor_, // This is the Stefan number
			const double k_, // ratio of thermal conductivities
			TimeStepper *timestepper = new Steady<0>
		) : FreeBoundaryGeometry(x0, x1, x2, y0, y1, timestepper), St(factor_), k(k_), ts_pt(timestepper) {

			geometry_index = add_internal_data(data_pt[0]);
			unpin_free_boundary();

			flux_data_pt = new Data(ts_pt, 1);

			for (unsigned int t = 0; t < ts_pt->nprev_values(); t++) flux_data_pt->set_value(t, 0, 0.0);
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

		void add_phase1_element(QUnsteadyHeatElement<2, X_ORDER> *elem) {
			phase1Elements.push_back(elem);
			for (unsigned int i = 0; i < elem->nnode(); i++)
				phase1_indexes.push_back(add_external_data(elem->node_pt(i)));
		}

		void add_phase2_element(QUnsteadyHeatElement<2, X_ORDER> *elem) {
			phase2Elements.push_back(elem);
			for (unsigned int i = 0; i < elem->nnode(); i++)
				phase2_indexes.push_back(add_external_data(elem->node_pt(i)));
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
		// residual equation is set as St * dh/dt - flux = 0
		// flux has to be set in before each newton step
		void fill_in_generic_residual_contribution(Vector<double>& residuals, DenseMatrix<double>& jacobian, bool compute_jacobian) {
			unsigned int ndofs = ndof();
			if (ndofs == 0) return;

			// printf("ndofs in FreeBoundaryElement: %u\n", ndofs);			
			// printf("phase1 indexes: %lu\n\t", phase1_indexes.size());
			// for (unsigned int i = 0; i < phase1_indexes.size(); i++)
			// 	printf("%u ", phase1_indexes[i]);
			// printf("\n");
			// printf("phase2 indexes: %lu\n\t", phase2_indexes.size());
			// for (unsigned int i = 0; i < phase2_indexes.size(); i++)
			// 	printf("%u ", phase2_indexes[i]);
			// printf("\n");

			Data *interface_data_pt = internal_data_pt(geometry_index);
			TimeStepper *interface_ts_pt = interface_data_pt->time_stepper_pt();

			int free_boundary_local_eqn_number = internal_local_eqn(geometry_index, free_boundary_index);

			if (compute_jacobian)
				jacobian(free_boundary_local_eqn_number, free_boundary_local_eqn_number) = St * interface_ts_pt->weight(1, 0);

			double tot_flux = 0.0;
			Vector<double> s = {1.0, 0.0};
			Vector<double> s2(2, 0.0);
			for (unsigned int e = 0; e < phase1Elements.size(); e++) {
				QUnsteadyHeatElement<2, X_ORDER> *elem = phase1Elements[e];
				unsigned int nnode = elem->nnode();
				
				// for jacobian stuff
				Shape psi(nnode), test(nnode);
				DShape dpsidx(nnode, 2), dtestdx(nnode, 2);

				// for flux calculations
				Shape psi_f(nnode);
				DShape dpsidx_f(nnode, 2);

				// elem->dshape_eulerian_at_knot(0, psi, dpsidx);
				// for (unsigned int n = 0; n < elem->nnode_1d(); n++)
				// 	test[n] = psi[n];

				elem->dshape_eulerian(s, psi_f, dpsidx_f);
				double flux = 0.0;

				for (unsigned int n = 0; n < nnode; n++) {
					flux += elem->nodal_value(n, 0) * dpsidx_f(n, 0);

					if (compute_jacobian) {
						// int local_unknown = elem->nodal_local_eqn(n, 0);
						int local_unknown = external_local_eqn(phase1_indexes[n], 0);
						if (local_unknown >= 0) {
							// printf("phase1 n=%u dpsidx=%f local_unknown=%d\n", n, dpsidx(n, 0), local_unknown);
							jacobian(free_boundary_local_eqn_number, local_unknown) = dpsidx_f(n, 0);

							double fraction = dynamic_cast<SpineNode *>(elem->node_pt(n))->fraction();
							// jacobian(local_unknown, free_boundary_local_eqn_number) = -fraction * interface_ts_pt->weight(1, 0) * psi(n) * dpsidx(n, 0);
							// printf("phase1 jac addition: %16.14f psi=%8.6f dpsidx=%8.6f frac=%5.3f\n", jacobian(local_unknown, free_boundary_local_eqn_number), psi(n), dpsidx(n, 0), fraction);

							double jac_contrib = 0.0;
							unsigned int n_ipt = elem->integral_pt()->nweight();
							for (unsigned int ipt = 0; ipt < n_ipt; ipt++) {
								for (unsigned d = 0; d < 2; d++) s2[d] = elem->integral_pt()->knot(ipt, d);
								double w = elem->integral_pt()->weight(ipt);

								double J = elem->dshape_eulerian_at_knot(ipt, psi, dpsidx);
								test = psi;
								dtestdx = dpsidx;

								double W = w * J;

								double mesh_vel_term = 0.0;
								double dudx_term = 0.0;

								for (unsigned int n2 = 0; n2 < nnode; n2++) {
									mesh_vel_term += interface_ts_pt->weight(1, 0) * psi(n2);
									dudx_term += elem->node_pt(n2)->value(0) * dpsidx(n2, 0);
									// printf("\te=%u n=%u n2=%u psi=%f dpsidx=%f T=%f tau=%f test=%f w=%e\n", e, n, n2, psi(n2), dpsidx(n2, 0), elem->node_pt(n2)->value(0), interface_ts_pt->weight(1, 0), test(n), W);
								}

								jac_contrib -= fraction * mesh_vel_term * dudx_term * test(n) * W;

							}

							jacobian(local_unknown, free_boundary_local_eqn_number) = jac_contrib;


						}
					}
				}

				tot_flux += flux;
			}

			s[0] = -1.0;
			for (unsigned int e = 0; e < phase2Elements.size(); e++) {
				QUnsteadyHeatElement<2, X_ORDER> *elem = phase2Elements[e];
				unsigned int nnode = elem->nnode();
				
				Shape psi(nnode), test(nnode), psi_f(nnode);
				DShape dpsidx(nnode, 2), dtestdx(nnode, 2), dpsidx_f(nnode, 2);

				// elem->dshape_eulerian_at_knot(0, psi, dpsidx);
				elem->dshape_eulerian(s, psi_f, dpsidx_f);
				double flux = 0.0;

				for (unsigned int n = 0; n < nnode; n++) {
					flux += elem->nodal_value(n, 0) * dpsidx_f(n, 0);
					if (compute_jacobian) {
						// int local_unknown = elem->nodal_local_eqn(n, 0);
						int local_unknown = external_local_eqn(phase2_indexes[n], 0);
						if (local_unknown >= 0) {
							// printf("phase2 n=%u dpsidx=%f local_unknown=%d\n", n, dpsidx(n, 0), local_unknown);
							jacobian(free_boundary_local_eqn_number, local_unknown) = -k * dpsidx_f(n, 0);
							
							double fraction = 1.0 - dynamic_cast<SpineNode *>(elem->node_pt(n))->fraction();
							// jacobian(local_unknown, free_boundary_local_eqn_number) = (fraction - 1.0) * interface_ts_pt->weight(1, 0) * psi(n) * dpsidx(n, 0);
							// printf("phase2 jac addition: %16.14f psi=%8.6f dpsidx=%8.6f frac=%5.3f\n", jacobian(local_unknown, free_boundary_local_eqn_number), psi(n), dpsidx(n, 0), fraction);

							double jac_contrib = 0.0;
							unsigned int n_ipt = elem->integral_pt()->nweight();
							for (unsigned int ipt = 0; ipt < n_ipt; ipt++) {
								for (unsigned d = 0; d < 2; d++) s2[d] = elem->integral_pt()->knot(ipt, d);
								double w = elem->integral_pt()->weight(ipt);

								double J = elem->dshape_eulerian_at_knot(ipt, psi, dpsidx);
								test = psi;
								dtestdx = dpsidx;

								double W = w * J;

								double mesh_vel_term = 0.0;
								double dudx_term = 0.0;

								for (unsigned int n2 = 0; n2 < nnode; n2++) {
									mesh_vel_term += interface_ts_pt->weight(1, 0) * psi(n2);
									dudx_term += elem->node_pt(n2)->value(0) * dpsidx(n2, 0);
									// printf("\te=%u n=%u n2=%u psi=%f dpsidx=%f T=%f tau=%f test=%f w=%e\n", e, n, n2, psi(n2), dpsidx(n2, 0), elem->node_pt(n2)->value(0), interface_ts_pt->weight(1, 0), test(n), W);
								}

								jac_contrib -= fraction * mesh_vel_term * dudx_term * test(n) * W;

							}

							jacobian(local_unknown, free_boundary_local_eqn_number) = jac_contrib;


						}
					}
				}

				tot_flux -= k * flux;
			}

			residuals[free_boundary_local_eqn_number] = St * interface_ts_pt->time_derivative(1, interface_data_pt, free_boundary_index) - tot_flux;
			

			// printf("residual: tot_flux: %16.14f\n", tot_flux);
			// residuals[free_boundary_local_eqn_number] = St * interface_ts_pt->time_derivative(1, interface_data_pt, free_boundary_index) - external_data_pt(flux_index)->value(0);

			// if (compute_jacobian) {
			// 	printf("FreeBoundaryElement Jacobian\n");
			// 	for (unsigned int i = 0; i < ndof(); i++) {
			// 		printf("\t");
			// 		for (unsigned int j = 0; j < ndof(); j++) {
			// 			printf("%8.6f\t", jacobian(i, j));
			// 		}
			// 		printf("\n");
			// 	}
			// }
		}
	};

class FreeBoundaryDomain : public Domain {
	private:
		FreeBoundaryElement * geometry;
		const unsigned int nelems;
		const unsigned int nx1;
		const unsigned int nx2;
		const unsigned int nx;
		const unsigned int ny;
	
	public:
		FreeBoundaryDomain(
			FreeBoundaryElement *geom,
			unsigned int nx1,
			unsigned int nx2,
			unsigned int ny
		) : geometry(geom), nelems((nx1+nx2)*ny), nx1(nx1), nx2(nx2), nx(nx1+nx2), ny(ny) {

			Macro_element_pt.resize(nelems);

			for (unsigned int i = 0; i < nelems; i++)
				Macro_element_pt[i] = new QMacroElement<2>(this, i);
		}

		void macro_element_boundary(
			const unsigned int &t, 
			const unsigned int &macro_i, 
			const unsigned int &dir, 
			const Vector<double> &zeta, 
			Vector<double> &r) 
		{
			using namespace QuadTreeNames;

			const double h = geometry->get_interface();

			const unsigned int yi = macro_i / nx;
			const unsigned int xi = macro_i % nx;

			const double dx1 = (h - geometry->x0()) / (double) nx1;
			const double dx2 = (geometry->x2() - h) / (double) nx2;
			const double dy = 1.0 / (double) ny;

			bool growing = xi < nx1;

			double trans_x = 0.5 * (zeta[0] + 1.0);

			double start = (growing) ? geometry->x0() + ((double) xi) * dx1 : h + ((double) (xi - nx1)) * dx2;
			double end = (growing) ? geometry->x0() + ((double) (xi + 1)) * dx1 : h + ((double) (xi - nx1 + 1)) * dx2;
			double x = start + trans_x * (end - start);

			switch (dir) {
				case N:
					r[0] = x;
					r[1] = ((double) (yi + 1)) * dy;
					break;
				case E:
					r[0] = end;
					r[1] = (trans_x + ((double) yi)) * dy;
					break;
				case S:
					r[0] = x;
					r[1] = ((double) yi) * dy;
					break;
				case W:
					r[0] = start;
					r[1] = (trans_x + ((double) yi)) * dy;
					break;
				default:
					printf("Invalid direction given to FreeBoundaryDomani\n");
					return;
			}
		}
};

template<class EL>
class TwoPhaseFreeBoundaryMesh : public RectangularQuadMesh<EL>,
								 public MacroElementNodeUpdateMesh {
	private:
		TimeStepper *ts_pt;
	public:
		TwoPhaseFreeBoundaryMesh(
			const unsigned int &nx1,
			const unsigned int &nx2,
			const unsigned int &ny,
			FreeBoundaryGeometry *geometry,
			FreeBoundaryDomain *domain,
			TimeStepper *timestepper = &Mesh::Default_TimeStepper
		) : RectangularQuadMesh<EL>(nx1+nx2, ny, geometry->x0(), geometry->x2(), geometry->y0(), geometry->y1(), timestepper), 
			ts_pt(timestepper) {

			this->set_nboundary(5);

			// new FreeBoundaryGeometry(x0, x1, x2, y0, y1, ts_pt);
			
			for (unsigned int e = 0; e < ny; e++) {
				EL *elem = dynamic_cast<EL *>(this->element_pt(nx1*(1+e) + nx2*e));
				unsigned int nnode = elem->nnode_1d();
				
				for (unsigned int n = 0; n < nnode; n++) {
					Node *node = elem->node_pt(nnode*n);
					
					this->convert_to_boundary_node(node);
					this->add_boundary_node(4, node);
				}
			}

			Vector<double> s(2), r(2);

			for (unsigned int e = 0; e < (nx1+nx2)*ny; e++) {
				EL *elem = dynamic_cast<EL *>(this->element_pt(e));

				elem->set_macro_elem_pt(domain->macro_element_pt(e));
				for (unsigned int n = 0; n < elem->nnode(); n++) {
					elem->local_fraction_of_node(n, s);

					s[0] = 2.0 * s[0] - 1.0;
					s[1] = 2.0 * s[1] - 1.0;

					domain->macro_element_pt(e)->macro_map(s, r);

					for (unsigned int i = 0; i < 2; i++)
						elem->node_pt(n)->x(i) = r[i];
				}

				Vector<GeomObject *> geom_object_pt(1);
				geom_object_pt[0] = geometry;

				elem->set_node_update_info(geom_object_pt);
			}
			
			this->setup_boundary_element_info();
		}
};

// Spine based mesh to be used with FreeBoundaryElement / FreeBoundaryGeometry
template<class EL>
class TwoPhaseFreeBoundarySpineMesh : public RectangularQuadMesh<EL>,
									  public SpineMesh {
	private:
		unsigned int nx1, nx2, nx, ny;
		FreeBoundaryGeometry *geometry;
		TimeStepper *ts_pt;
	
	public:
		TwoPhaseFreeBoundarySpineMesh(
			const unsigned int &nx1,
			const unsigned int &nx2,
			const unsigned int &ny,
			FreeBoundaryGeometry *geometry,
			TimeStepper *timestepper = &Mesh::Default_TimeStepper
		) : RectangularQuadMesh<EL>(nx1+nx2, ny, geometry->x0(), geometry->x2(), geometry->y0(), geometry->y1(), timestepper),
			nx1(nx1), nx2(nx2), nx(nx1+nx2), ny(ny), geometry(geometry), ts_pt(timestepper) {

			// set interface nodes
			this->set_nboundary(5);
			for (unsigned int e = 0; e < ny; e++) {
				EL *elem = dynamic_cast<EL *>(this->element_pt(nx1*(1+e) + nx2*e));
				unsigned int nnode = elem->nnode_1d();

				for (unsigned int n = 0; n < nnode; n++) {
					Node *node = elem->node_pt(nnode*n);

					this->convert_to_boundary_node(node);
					this->add_boundary_node(4, node);
				}
			}

			construct_spines();

			for (unsigned int n = 0; n < nnode(); n++) {
				spine_node_update(node_pt(n));
			}

			this->setup_boundary_element_info();
		}

		// create spines and associate each node with a spine
		// node_update_fct_id is used to differentiate between phases
		void construct_spines() {
			unsigned int np = finite_element_pt(0)->nnode_1d();
			unsigned int nspine = (np - 1) * ny + 1;
			Spine_pt.reserve(nspine);

			// loop through vertical elements
			unsigned int yi = 0;
			for (yi = 0; yi < ny-1; yi++) {
				// loop through nodes in vertical direction
				for (unsigned int s = 0; s < np-1; s++) {
					// create spine and set parameters
					Spine *spine = new Spine(geometry->x1());
					spine->spine_height_pt()->pin(0);
					Spine_pt.push_back(spine);

					Vector<double> parameters = {((double)yi + (double)s/(double)(np-1)) / (double)ny};
					spine->set_geom_parameter(parameters);

					Vector<GeomObject *> geom_object_pt = {geometry};
					spine->set_geom_object_pt(geom_object_pt);

					// loop through nodes in phase 1
					unsigned int xi = 0;
					for (xi = 0; xi < nx1-1; xi++) {
						for (unsigned n = 0; n < np-1; n++) {
							SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
							node_pt->spine_pt() = spine;
							node_pt->fraction() = ((double)xi + (double)n/(double)(np-1)) / (double)nx1;
							node_pt->spine_mesh_pt() = this;
							node_pt->node_update_fct_id() = 0;
						}
					}

					// do all nodes in last element in phase 1
					xi = nx1-1;
					for (unsigned n = 0; n < np; n++) {
						SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
						node_pt->spine_pt() = spine;
						node_pt->fraction() = ((double)xi + (double)n/(double)(np-1)) / (double)nx1;
						node_pt->spine_mesh_pt() = this;
						node_pt->node_update_fct_id() = 0;
					}

					// loop through nodes in phase 2
					for (unsigned int xi = nx1; xi < nx-1; xi++) {
						for (unsigned n = 0; n < np-1; n++) {
							SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
							node_pt->spine_pt() = spine;
							node_pt->fraction() = ((double)(xi-nx1) + (double)n/(double)(np-1)) / (double)nx2;
							node_pt->spine_mesh_pt() = this;
							node_pt->node_update_fct_id() = 1;
						}
					}

					// do all nodes in last element in phase 2
					xi = nx-1;
					for (unsigned n = 0; n < np; n++) {
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
				// create spine and set parameters
				Spine *spine = new Spine(geometry->x1());
				spine->spine_height_pt()->pin(0);
				Spine_pt.push_back(spine);

				Vector<double> parameters = {((double)yi + (double)s/(double)(np-1)) / (double)ny};
				spine->set_geom_parameter(parameters);

				Vector<GeomObject *> geom_object_pt = {geometry};
				spine->set_geom_object_pt(geom_object_pt);

				// loop through nodes in phase 1
				unsigned int xi = 0;
				for (xi = 0; xi < nx1-1; xi++) {
					for (unsigned n = 0; n < np-1; n++) {
						SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
						node_pt->spine_pt() = spine;
						node_pt->fraction() = ((double)xi + (double)n/(double)(np-1)) / (double)nx1;
						node_pt->spine_mesh_pt() = this;
						node_pt->node_update_fct_id() = 0;
					}
				}

				// do all nodes in last element in phase 1
				xi = nx1-1;
				for (unsigned n = 0; n < np; n++) {
					SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
					node_pt->spine_pt() = spine;
					node_pt->fraction() = ((double)xi + (double)n/(double)(np-1)) / (double)nx1;
					node_pt->spine_mesh_pt() = this;
					node_pt->node_update_fct_id() = 0;
				}

				// loop through nodes in phase 2
				for (unsigned int xi = nx1; xi < nx-1; xi++) {
					for (unsigned n = 0; n < np-1; n++) {
						SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
						node_pt->spine_pt() = spine;
						node_pt->fraction() = ((double)(xi-nx1) + (double)n/(double)(np-1)) / (double)nx2;
						node_pt->spine_mesh_pt() = this;
						node_pt->node_update_fct_id() = 1;
					}
				}

				// do all nodes in last element in phase 2
				xi = nx-1;
				for (unsigned n = 0; n < np; n++) {
					SpineNode *node_pt = element_node_pt(yi*nx + xi, s*np + n);
					node_pt->spine_pt() = spine;
					node_pt->fraction() = ((double)(xi-nx1) + (double)n/(double)(np-1)) / (double)nx2;
					node_pt->spine_mesh_pt() = this;
					node_pt->node_update_fct_id() = 1;
				}
			}

			printf("Created %lu/%u spines\n\n", Spine_pt.size(), nspine);
		}

		// assign position for each node based on location of interface
		void spine_node_update(SpineNode *node_pt) {
			FreeBoundaryGeometry *geom = dynamic_cast<FreeBoundaryGeometry *>(node_pt->spine_pt()->geom_object_pt(0));
			double frac = node_pt->fraction();
			Vector<double> zeta(1), r(1);

			zeta[0] = 0.0;
			geom->position(zeta, r);

			switch (node_pt->node_update_fct_id()) {
				case 0: node_pt->x(0) = (1.0 - frac) * geom->x0() + frac * r[0];	break; // solid phase
				case 1: node_pt->x(0) = (1.0 - frac) * r[0] + frac * geom->x2();	break; // liquid phase
				default:
					printf("Invalid node update function\n");
					break;
			}
		}
};

#endif