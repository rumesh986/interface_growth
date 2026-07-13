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

			for (unsigned int i = 0; i < nnode; i++)
				add_internal_data(geom_data_pt(i));
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

		inline void fill_in_generic_residual_contribution(Vector<double> &residuals, DenseMatrix<double> &jacobian, bool compute_jacobian) {}
};

template<class EL>
class FreeBoundaryFluxElement : public UnsteadyHeatFluxElement<EL> {
	private:
		double _St = 1.0;
		double _gamma = 0.0;
		Vector<int> geom_indices;
		const unsigned int T_index = 0;
		const unsigned int Kx_index = 1;
		const unsigned int Ky_index = 2;

	public:
		FreeBoundaryFluxElement(
			EL *bulk_elem,
			unsigned int face_index
		) : UnsteadyHeatFluxElement<EL>(bulk_elem, face_index) {
			
			geom_indices.reserve(this->nnode());
			for (unsigned int n = 0; n < this->nnode(); n++) {
				SpineNode *node = dynamic_cast<SpineNode *>(this->node_pt(n));
				geom_indices[n] = this->add_external_data(node->spine_pt()->geom_data_pt(1));
			}

			Vector<unsigned int> temp(this->nnode(), 2);
			this->resize_nodes(temp);
		}

		~FreeBoundaryFluxElement() {}

		double & St() {
			return _St;
		}
	
		double & gamma() {
			return _gamma;
		}

	protected:
		inline void fill_in_contribution_to_residuals(Vector<double> &residuals) {
			fill_in_generic_residual_contribution_ust_heat_flux(residuals, GeneralisedElement::Dummy_matrix, false);
		}

		inline void fill_in_contribution_to_jacobian(Vector<double> &residuals, DenseMatrix<double> &jacobian) {
			fill_in_generic_residual_contribution_ust_heat_flux(residuals, jacobian, true);
		}

		double get_angle(Vector<double> &a, Vector<double> &b) {
			double dot = VectorHelpers::dot(a, b);
			double mag = sqrt((a[0]*a[0] + a[1]*a[1])) * sqrt(b[0]*b[0] + b[1]* b[1]);
			double inp = dot / mag;
			if (fabs(inp) > 1.0) {
				printf("\nWeird input received!\n");
				printf("\tinp=%e dot=%e mag=%e\n", inp, dot, mag);
			}
			return std::acos(inp);
		}

	public:
		inline void fill_in_generic_residual_contribution_ust_heat_flux(Vector<double> &residuals, DenseMatrix<double> &jacobian, bool compute_jacobian) {
			if (this->ndof() == 0) return;

			Vector<double> s(1, 0.0);
			Shape phi(this->nnode()), psi(this->nnode());
			DShape dphi(this->nnode(), 1);

			for (unsigned int ipt = 0; ipt < this->integral_pt()->nweight(); ipt++) {
				for (unsigned int i = 0; i < s.size(); i++) s[i] = this->integral_pt()->knot(ipt, i);

				double J = this->shape_and_test(s, psi, phi);
				this->dshape_local(s, phi, dphi);
				double w = this->integral_pt()->weight(ipt);
				double JW = J * w;

				Vector<double> tangent(2, 0.0), th(2, 0.0), nh(2, 0.0);
				for (unsigned int l = 0; l < this->nnode(); l++) {
					for (unsigned int i = 0; i < 2; i++) {
						tangent[i] += this->nodal_position(l, i) * dphi(l, 0);
					}
				}
				
				double tangent_mag = VectorHelpers::magnitude(tangent);
				for (unsigned int i = 0; i < 2; i++) th[i] = tangent[i] / tangent_mag;
				nh[0] =  th[1];
				nh[1] = -th[0];

				Vector<DenseMatrix<double>> dthdX(this->nnode()), dnhdX(this->nnode());
				if (compute_jacobian) {
					for (unsigned int l = 0; l < this->nnode(); l++) {
						dthdX[l].resize(2, 2, 0.0);
						dnhdX[l].resize(2, 2, 0.0);

						dthdX[l](0, 0) = dphi(l, 0) * (1.0 - th[0] * th[0]) / tangent_mag;
						dthdX[l](0, 1) = dphi(l, 0) * (0.0 - th[0] * th[1]) / tangent_mag;
						dthdX[l](1, 0) = dphi(l, 0) * (0.0 - th[1] * th[0]) / tangent_mag;
						dthdX[l](1, 1) = dphi(l, 0) * (1.0 - th[1] * th[1]) / tangent_mag;

						for (unsigned int j = 0; j < 2; j++) {
							dnhdX[l](0, j) =  dthdX[l](1, j);
							dnhdX[l](1, j) = -dthdX[l](0, j);
						}
					}
				}

				for (unsigned int l = 0; l < this->nnode(); l++) {
					Data *geom = this->external_data_pt(geom_indices[l]);
					double dhdt = geom->time_stepper_pt()->time_derivative(1, geom, 0);
					Vector<double> Kappa(2, 0.0);
					Kappa[0] = this->nodal_value(l, Kx_index);
					Kappa[1] = this->nodal_value(l, Ky_index);
					double kappa = VectorHelpers::dot(Kappa, nh);

					int X_eqn = this->external_local_eqn(geom_indices[l], 0);
					int T_eqn = this->nodal_local_eqn(l, T_index);
					int Kx_eqn = this->nodal_local_eqn(l, Kx_index);
					int Ky_eqn = this->nodal_local_eqn(l, Ky_index);

					residuals[X_eqn] += phi(l) * this->nodal_value(l, T_index) * JW - phi(l) * _gamma * kappa * JW;
					residuals[T_eqn] -= phi(l) * _St * dhdt * nh[0] * JW;
					residuals[Kx_eqn] += phi(l) * this->nodal_value(l, Kx_index) * JW + dphi(l, 0) * th[0] * w;
					residuals[Ky_eqn] += phi(l) * this->nodal_value(l, Ky_index) * JW + dphi(l, 0) * th[1] * w;
					
					if (compute_jacobian) {
						jacobian(X_eqn, T_eqn) += phi(l) * JW;
						
						jacobian(X_eqn, Kx_eqn) -= phi(l) * _gamma * nh[0] * JW;
						jacobian(X_eqn, Ky_eqn) -= phi(l) * _gamma * nh[1] * JW;
						for (unsigned int p = 0; p < this->nnode(); p++) {
							int P_eqn = this->external_local_eqn(geom_indices[p], 0);
							jacobian(X_eqn, P_eqn) -= phi(l) * _gamma * (Kappa[0] * dnhdX[p](0, 0) + Kappa[1] * dnhdX[p](1, 0)) * JW;
						}

						jacobian(T_eqn, X_eqn) -= phi(l) * _St * geom->time_stepper_pt()->weight(1, 0) * nh[0] * JW;
						for (unsigned int p = 0; p < this->nnode(); p++) {
							int P_eqn = this->external_local_eqn(geom_indices[p], 0);
							jacobian(T_eqn, P_eqn) -= phi(l) * _St * dhdt * dnhdX[p](0, 0) * JW;
						}

						jacobian(Kx_eqn, Kx_eqn) += phi(l) * JW;
						jacobian(Ky_eqn, Ky_eqn) += phi(l) * JW;
						for (unsigned int p = 0; p < this->nnode(); p++) {
							int P_eqn = this->external_local_eqn(geom_indices[p], 0);
							jacobian(Kx_eqn, P_eqn) += dphi(l, 0) * dthdX[p](0, 0) * w;
							jacobian(Ky_eqn, P_eqn) += dphi(l, 0) * dthdX[p](1, 0) * w;
						}
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
			_geometry = new FreeBoundaryElement(xs, ys, nspine, timestepper);
			_geometry->free_boundary_index() = _free_boundary_index;

			construct_spines(xs[1]);

			for (unsigned int n = 0; n < nnode(); n++) {
				spine_node_update(node_pt(n));
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