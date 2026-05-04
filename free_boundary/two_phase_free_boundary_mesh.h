#ifndef __TWO_PHASE_FREE_BOUNDARY_MESH_H__
#define __TWO_PHASE_FREE_BOUNDARY_MESH_H__

#include "generic.h"
#include "unsteady_heat.h"

using namespace oomph;

class FreeBoundaryGeometry : public GeomObject {
	protected:
		Vector<Data *> _data_pt;
		bool destroy_geom_data = false;
		
		double _xs[2];
		double _ys[2];

	public:
		FreeBoundaryGeometry(
			double *xs,
			double *ys,
			unsigned int nnode,
			TimeStepper *timestepper
		) : GeomObject(2, 2, timestepper) {
			
			_xs[0] = xs[0];
			_xs[1] = xs[2];
			_ys[0] = ys[0];
			_ys[1] = ys[1];

			_data_pt.resize(nnode);

			for (unsigned int i = 0; i < nnode; i++) {
				_data_pt[i] = new Data(time_stepper_pt(), 1);

				for (unsigned int t = 0; t < time_stepper_pt()->nprev_values(); t++)
					_data_pt[i]->set_value(t, 0, xs[1]);
				
				_data_pt[i]->pin_all();
			}

			destroy_geom_data = true;
		}

		Data *geom_data_pt(const unsigned int &i) {return _data_pt[i];}
		unsigned int ngeom_data() const {return _data_pt.size();}

		double& x_min() {return _xs[0];}
		double& x_max() {return _xs[1];}
		double& y_min() {return _ys[0];}
		double& y_max() {return _ys[1];}

		double& get_interface(const unsigned int &t, const unsigned int &i) {return *_data_pt[i]->value_pt(t, 0);}
		double& get_interface(const unsigned int &i) {return get_interface(0, i);}

		void position(const unsigned int &t, const Vector<double> &zeta, Vector<double> &r) const {
			r[0] = _data_pt[0]->value(0);
		}

		void position(const Vector<double> &zeta, Vector<double>& r) const {
			position(0, zeta, r);
		}

		void pin_free_boundary() {
			for (unsigned int i = 0; i < ngeom_data(); i++)
				_data_pt[i]->pin_all();
		}

		void unpin_free_boundary() {
			for (unsigned int i = 0; i < ngeom_data(); i++)
				_data_pt[i]->unpin_all();
		}
};

class FreeBoundaryElement : public GeneralisedElement,
							public FreeBoundaryGeometry {
	private:
		unsigned int boundary_index, geometry_index;
		Vector<unsigned int> node_indices;
		std::map<unsigned int, unsigned int> node_int_map;
		std::map<unsigned int, unsigned int> node_ext_map;
		std::map<Node *, Data *> node_data_map;

	public:
		FreeBoundaryElement(
			double *xs,
			double *ys,
			unsigned int nnode,
			TimeStepper *timestepper
		) : FreeBoundaryGeometry(xs, ys, nnode, timestepper) {

			unpin_free_boundary();
			destroy_geom_data = true;
		}

		void add_node(Node *node_pt, unsigned int i) {
			node_int_map[i] = add_internal_data(geom_data_pt(i));
			node_ext_map[i] = add_external_data(node_pt);
			node_data_map[node_pt] = geom_data_pt(i);
		}

		Data*& get_data_for_node(Node *node_pt) {
			return node_data_map[node_pt];
		}

		unsigned int& free_boundary_index() {
			return boundary_index;
		}

		void calculate_predicted_values() {
			for (unsigned int i = 0; i < ngeom_data(); i++) {
				time_stepper_pt()->calculate_predicted_values(geom_data_pt(i));
			}
		}
	
	protected:
		inline void get_residuals(Vector<double> &residuals) {
			residuals.initialise(0.0);
			DenseMatrix<double> placeholder(1);
			fill_in_generic_residual_contribution(residuals, placeholder, false);
		}

		inline void get_jacobian(Vector<double> &residuals, DenseMatrix<double> &jacobian) {
			residuals.initialise(0.0);
			jacobian.initialise(0.0);
			fill_in_generic_residual_contribution(residuals, jacobian, true);
		}

		inline void fill_in_generic_residual_contribution(Vector<double> &residuals, DenseMatrix<double> &jacobian, bool compute_jacobian) {
			if (ndof() == 0) return;

			for (unsigned int i = 0; i < ngeom_data(); i++) {
				int local_eqn = internal_local_eqn(node_int_map[i], 0);
				if (local_eqn < 0) continue;

				Node *node_pt = dynamic_cast<Node *>(external_data_pt(node_ext_map[i]));
				if (node_pt == NULL || !node_pt->is_on_boundary(boundary_index)) continue; 

				residuals[local_eqn] = node_pt->value(0);

				if (compute_jacobian) {
					int local_unknown = external_local_eqn(node_ext_map[i], 0);
					if (local_unknown < 0) continue;

					jacobian(local_eqn, local_unknown) = 1.0;
				}
			}
		}
};

template<class EL>
class FreeBoundaryFluxElement : public UnsteadyHeatFluxElement<EL> {
	private:
		double St;
		std::map<Node *, int> geom_indices;

	public:
		FreeBoundaryFluxElement(
			EL *bulk_elem,
			unsigned int face_index,
			double _St
		) : UnsteadyHeatFluxElement<EL>(bulk_elem, face_index), St(_St) {}

		void add_node_data(Node *node_pt, Data *data_pt) {
			geom_indices[node_pt] = this->add_external_data(data_pt);
		}
	
	protected:
		inline void fill_in_contribution_to_residuals(Vector<double> &residuals) {
			fill_in_generic_residual_contribution_ust_heat_flux(residuals, GeneralisedElement::Dummy_matrix, false);
		}

		inline void fill_in_contribution_to_jacobian(Vector<double> &residuals, DenseMatrix<double> &jacobian) {
			fill_in_generic_residual_contribution_ust_heat_flux(residuals, jacobian, true);
		}

		inline void fill_in_generic_residual_contribution_ust_heat_flux(Vector<double> &residuals, DenseMatrix<double> &jacobian, bool compute_jacobian) {
			if (this->ndof() == 0) return;

			Vector<double> s(1), normal(2);
			Shape phi(this->nnode()), psi(this->nnode());
			
			for (unsigned int ipt = 0; ipt < this->integral_pt()->nweight(); ipt++) {
				for (unsigned int i = 0; i < s.size(); i++) s[i] = this->integral_pt()->knot(ipt, i);

				double J = this->shape_and_test(s, phi, psi);
				double W = J * this->integral_pt()->weight(ipt);
				this->outer_unit_normal(s, normal);
				double n_mag = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
				// printf("normal: n0=%16.14f n1=%16.14f mag=%16.14f\n", n[0], n[1], n_mag);
				// for error checking...
				if (n_mag - 1.0 > 1e-6) {
					printf("WARNING: normal not quite a unit\n\n\n");
				}


				for (unsigned int n = 0; n < this->nnode(); n++) {
					int local_eqn = this->nodal_local_eqn(n, 0);
					if (local_eqn < 0) continue;

					Node *node = this->node_pt(n);
					Data *geom = this->external_data_pt(geom_indices[node]);

					double dhdt = geom->time_stepper_pt()->time_derivative(1, geom, 0);
					residuals[local_eqn] -= phi(n) * St * dhdt * normal[0] * W;

					if (compute_jacobian) {
						int h_eqn = this->external_local_eqn(geom_indices[node], 0);
						if (h_eqn < 0) continue;

						jacobian(local_eqn, h_eqn) -= phi(n) * St * geom->time_stepper_pt()->weight(1, 0) * normal[0] * W;
					}
				}
			}
		}
};

template<class EL>
class TwoPhaseFreeBoundarySpineMesh : public RectangularQuadMesh<EL>,
									  public SpineMesh {
	private:
		const unsigned int _free_boundary_index = 4;
		const unsigned int nx1, nx2, nx, ny;
		const bool _periodic;

		FreeBoundaryElement *_geometry;
	public:
		TwoPhaseFreeBoundarySpineMesh(
			unsigned int *nxs,
			unsigned int ny,
			double *xs,
			double *ys,
			double gamma,
			bool periodic,
			TimeStepper *timestepper
		) : RectangularQuadMesh<EL>(nxs[0]+nxs[1], ny, xs[0], xs[2], ys[0], ys[1], timestepper),
			nx1(nxs[0]), nx2(nxs[1]), nx(nxs[0]+nxs[1]), ny(ny), _periodic(periodic) {

			// periodicity in y
			if (_periodic) {
				for (unsigned int n = 0; n < nboundary_node(0); n++) {
					boundary_node_pt(2, n)->make_periodic(boundary_node_pt(0, n));
				}
			}

			this->set_nboundary(5);
			for (unsigned int e = 0; e < ny; e++) {
				EL *elem = dynamic_cast<EL *>(this->element_pt((nxs[0]+nxs[1])*e + nxs[0]));
				unsigned int nnode_1d = elem->nnode_1d();
				for (unsigned int n = 0; n < nnode_1d; n++) {
					Node *node = elem->node_pt(n*nnode_1d);
					
					this->convert_to_boundary_node(node);
					this->add_boundary_node(_free_boundary_index, node);	
				}
			}
			
			unsigned int nspine = nboundary_node(_free_boundary_index);
			if (_periodic) nspine--;
			_geometry = new FreeBoundaryElement(xs, ys, nspine, gamma, timestepper);
			_geometry->free_boundary_index() = _free_boundary_index;

			construct_spines(xs[1]);

			for (unsigned int n = 0; n < nnode(); n++) {
				spine_node_update(node_pt(n));
			}

			for (unsigned int s = 0; s < Spine_pt.size(); s++) {
				printf("Spine %u: %p\n", s, Spine_pt[s]->geom_data_pt(1));
			}
		}

		FreeBoundaryElement* geometry() const {
			return _geometry;
		}

		Spine *_create_new_spine(unsigned int s, const double h, const Vector<double> params) {
			Spine *spine = new Spine(h);
			// spine->spine_height_pt()->pin(0);
			spine->spine_height_pt()->pin_all();
			spine->add_geom_data_pt(_geometry->geom_data_pt(s));
			spine->set_geom_parameter(params);
			Spine_pt.push_back(spine);

			// printf("Creating spine %d (%p) with data %p -> %p\n", s, spine, _geometry->geom_data_pt(s), spine->geom_data_pt(1));

			return spine;
		}

		void _add_spine_to_node(
			const unsigned int e, 
			const unsigned int n, 
			const double fraction,
			const unsigned int fct_id,
			Spine *spine
		) {
			SpineNode *node_pt = element_node_pt(e, n);
			node_pt->spine_pt() = spine;
			node_pt->fraction() = fraction;
			node_pt->spine_mesh_pt() = this;
			node_pt->node_update_fct_id() = fct_id;
		}

		void construct_spines(double h) {
			unsigned int nnode_1d = finite_element_pt(0)->nnode_1d();
			unsigned int nspine = nboundary_node(_free_boundary_index);
			if (_periodic) nspine--;
			Spine_pt.reserve(nspine);

			unsigned int yi;
			for (yi = 0; yi < ny; yi++) {
				for (unsigned int s = 0; s < nnode_1d-1; s++) {
					Vector<double> parameters = {((double)yi + (double)s/(double)(nnode_1d-1)) / (double)ny};
					Spine *spine = _create_new_spine(yi * (nnode_1d-1) + s, h, parameters);

					unsigned int xi;
					for (xi = 0; xi < nx1-1; xi++) {
						for (unsigned int n = 0; n < nnode_1d-1; n++) {
							_add_spine_to_node(yi*nx + xi, s*nnode_1d+n, ((double)xi + (double)n/(double)(nnode_1d-1))/(double)nx1, 0, spine);
						}
					}

					xi = nx1-1;
					for (unsigned int n = 0; n < nnode_1d; n++) {
						_add_spine_to_node(yi*nx + xi, s*nnode_1d+n, ((double)xi + (double)n/(double)(nnode_1d-1))/(double)nx1, 0, spine);
					}

					for (xi = nx1; xi < nx-1; xi++) {
						for (unsigned int n = 0; n < nnode_1d-1; n++) {
							_add_spine_to_node(yi*nx + xi, s*nnode_1d+n, ((double)(xi-nx1) + (double)n/(double)(nnode_1d-1))/(double)nx2, 1, spine);
						}
					}

					xi = nx-1;
					for (unsigned int n = 0; n < nnode_1d; n++) {
						_add_spine_to_node(yi*nx + xi, s*nnode_1d+n, ((double)(xi-nx1) + (double)n/(double)(nnode_1d-1))/(double)nx2, 1, spine);
					}
				}
			}

			yi = ny-1;
			unsigned int s = nnode_1d-1;
			Vector<double> parameters = {((double)yi + (double)s/(double)(nnode_1d-1)) / (double)ny};
			Spine *spine;
			if (_periodic) {
				spine = Spine_pt[0];
			} else {
				spine = _create_new_spine(yi * (nnode_1d-1) + s, h, parameters);
			}

			unsigned int xi;
			for (xi = 0; xi < nx1-1; xi++) {
				for (unsigned int n = 0; n < nnode_1d-1; n++) {
					_add_spine_to_node(yi*nx + xi, s*nnode_1d+n, ((double)xi + (double)n/(double)(nnode_1d-1))/(double)nx1, 0, spine);
				}
			}

			xi = nx1-1;
			for (unsigned int n = 0; n < nnode_1d; n++) {
				_add_spine_to_node(yi*nx + xi, s*nnode_1d+n, ((double)xi + (double)n/(double)(nnode_1d-1))/(double)nx1, 0, spine);
			}

			for (xi = nx1; xi < nx-1; xi++) {
				for (unsigned int n = 0; n < nnode_1d-1; n++) {
					_add_spine_to_node(yi*nx + xi, s*nnode_1d+n, ((double)(xi-nx1) + (double)n/(double)(nnode_1d-1))/(double)nx2, 1, spine);
				}
			}

			xi = nx-1;
			for (unsigned int n = 0; n < nnode_1d; n++) {
				_add_spine_to_node(yi*nx + xi, s*nnode_1d+n, ((double)(xi-nx1) + (double)n/(double)(nnode_1d-1))/(double)nx2, 1, spine);
			}

			printf("Created %lu/%u spines\n\n", Spine_pt.size(), nspine);
		}

		const unsigned int& free_boundary_index() const {
			return _free_boundary_index;
		}

		void spine_node_update(SpineNode * node_pt) {
			Data *geom_data = node_pt->spine_pt()->geom_data_pt(1);

			double h = geom_data->value(0);
			double frac = node_pt->fraction();

			double x0 = _geometry->x_min();
			double x1 = _geometry->x_max();

			switch (node_pt->node_update_fct_id()) {
				case 0: node_pt->x(0) = x0 + frac * (h - x0); break;
				case 1: node_pt->x(0) = h  + frac * (x1 - h); break;
				default:
					printf("Invalid node updating function\n\n");
					break;
			}
		}
};

#endif