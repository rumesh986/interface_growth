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

		const unsigned int& free_index() {return free_boundary_index;}
		
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

		void unpin_free_boundary() {
			data_pt[0]->unpin(free_boundary_index);
		}

		void pin_free_boundary() {
			data_pt[0]->pin(free_boundary_index);
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

template<class EL>
class FreeBoundaryFluxElement : public UnsteadyHeatFluxElement<EL> {
	private:
		double St;
		// FreeBoundaryGeometry *geom;

		// Vector<unsigned int> data_indexes;
		// unsigned int data_index;
		EL *solid_elem, *liquid_elem;

		Data *position;
	
	public:
		FreeBoundaryFluxElement(
			EL *_solid_elem,
			EL *_liquid_elem, 
			unsigned int face_index,
			double _St,
			// FreeBoundaryGeometry *_geom,
			TimeStepper *_ts_pt
		) : UnsteadyHeatFluxElement<EL>(_solid_elem, face_index), St(_St), solid_elem(_solid_elem), liquid_elem(_liquid_elem) {//, geom(_geom) {
			this->time_stepper_pt() = _ts_pt;
			printf("in constructor\n");
			printf("time stepper pt: %p\n", this->time_stepper_pt());
			// data_index = this->add_external_data(geom->geom_data_pt(0));

			position = new Data(this->time_stepper_pt(), 1, true);
		}

		void trial() {

			printf("Flux elem positions\n");
			for (unsigned int n = 0; n < this->nnode(); n++) {
				Vector<double> x(2);
				this->node_pt(n)->position(x);
				printf("\tNode %u x0=%f x1=%f\n", n, x[0], x[1]);
			}
			// printf("In FreeBoundaryFluxElement trial function\n");

			// printf("ndofs: %u\n", this->ndof());

			// int local_eqn = this->external_local_eqn(data_index, geom->free_index());
			// printf("Got local equation %d\n", local_eqn);
			// int global_eqn = this->eqn_number(local_eqn);


			// printf("Flux element global eqn number: %d\n", global_eqn);
		}

		 /// Compute the element residual vector
		inline void fill_in_contribution_to_residuals(Vector<double>& residuals) {
			// Call the generic residuals function with flag set to 0
			// using a dummy matrix argument
			fill_in_generic_residual_contribution_ust_heat_flux(residuals, GeneralisedElement::Dummy_matrix, 0);
		}


		/// Compute the element's residual vector and its Jacobian matrix
		inline void fill_in_contribution_to_jacobian(Vector<double>& residuals, DenseMatrix<double>& jacobian) {
			// Call the generic routine with the flag set to 1
			fill_in_generic_residual_contribution_ust_heat_flux(residuals, jacobian, 1);
		}

		void fill_in_generic_residual_contribution_ust_heat_flux(Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned int compute_jacobian) {
			unsigned int ndofs = this->ndof();
			if (ndofs == 0) return;

			unsigned int n_node = this->nnode();
			unsigned int n_ipt = this->integral_pt()->nweight();
			unsigned int n_dim = this->node_pt(0)->ndim();

			Vector<double> s(n_dim-1);
			Shape shape(n_node), test(n_node);

			const unsigned int T_idx = solid_elem->u_index_ust_heat();

			for (unsigned int ipt = 0; ipt < n_ipt; ipt++) {
				for (unsigned int i = 0; i < n_dim-1; i++) 
					s[i] = this->integral_pt()->knot(ipt, i);

				double J = this->shape_and_test(s, shape, test);
				double W = J * this->integral_pt()->weight(ipt);

				Vector<double> x_int(n_dim, 0.0);

				for (unsigned int n = 0; n < n_node; n++) {
					for (unsigned int i = 0; i < n_dim; i++)
						x_int[i] += this->nodal_position(n, i) * shape[n];
				}

				Vector<double> flux(2);
				solid_elem->get_flux(x_int, flux);

				for (unsigned int n1 = 0; n1 < n_node; n1++) {
					int local_eqn = this->nodal_local_eqn(n1, T_idx);
					if (local_eqn < 0) continue;

					residuals[local_eqn] -= solid_elem->beta() * test(n1) * flux[0] * W;
					
					if (compute_jacobian) {
						for (unsigned int n2 = 0; n2 < n_node; n2++) {
							int local_unknown = this->nodal_local_eqn(n2, T_idx);
							if (local_unknown < 0) continue;
							
							jacobian(local_eqn, local_unknown) -= solid_elem->beta() * test(n1) * W; // needs dshape(n2)
							// printf("n1=%u n2=%u res=%f jac=%f beta=%f\n", n1, n2, residuals[local_eqn], jacobian(local_eqn, local_unknown), solid_elem->beta());
						}
					}
				}
			}
		}
};

template<class EL>
class TwoPhaseElasticMesh : public ElasticRectangularQuadMesh<EL> {
	private:
		unsigned int nx1, nx2, nx, ny;
		double x0, x1, x2, y0, y1;
		TimeStepper *ts_pt;
	
	public:
		TwoPhaseElasticMesh(
			const unsigned int &nx1,
			const unsigned int &nx2,
			const unsigned int &ny,
			const double &x0,
			const double &x1,
			const double &x2,
			const double &y0,
			const double &y1,
			TimeStepper *timestepper = &Mesh::Default_TimeStepper
		) : RectangularQuadMesh<EL>(nx1+nx2, ny, x0, x2, y0, y1, timestepper),
			ElasticRectangularQuadMesh<EL>(nx1+nx2, ny, x2, y1, timestepper),
			nx1(nx1), nx2(nx2), nx(nx1+nx2), ny(ny),  x0(x0), x1(x1), x2(x2), y0(y0), y1(y1), ts_pt(timestepper) {

			this->set_nboundary(5);
			for (unsigned int e = 0; e < ny; e++) {
				EL *elem = dynamic_cast<EL *>(this->element_pt(nx1 + nx*e));
				unsigned int nnode = elem->nnode_1d();

				for (unsigned int n = 0; n < nnode; n++) {
					Node *node = elem->node_pt(nnode * n);

					this->convert_to_boundary_node(node);
					this->add_boundary_node(4, node);
				}
			}

			this->setup_boundary_element_info();
		}
};

// doesn't work, don't bother
// compiler doesn't know if it should use SolidNode::node_pt(...) or SpineNode::node_pt(...)
// for element->node_pt(...)
template<class EL>
class TwoPhaseElasticSpineMesh : public TwoPhaseFreeBoundarySpineMesh<EL>,
								 public TwoPhaseElasticMesh<EL> {
	public:
		TwoPhaseElasticSpineMesh(
			const unsigned int &nx1,
			const unsigned int &nx2,
			const unsigned int &ny,
			FreeBoundaryGeometry *geometry,
			TimeStepper *timestepper = &Mesh::Default_TimeStepper
		) : RectangularQuadMesh<EL>(nx1+nx2, ny, geometry->x0(), geometry->x2(), geometry->y0(), geometry->y1(), timestepper),
			TwoPhaseFreeBoundarySpineMesh<EL>(nx1, nx2, ny, geometry, timestepper),
			TwoPhaseElasticMesh<EL>(nx1, nx2, ny, geometry->x0(), geometry->x1(), geometry->x2(), geometry->y0(), geometry->y1(), timestepper) {

			printf("IT compiles I guess");
		}
};

#endif