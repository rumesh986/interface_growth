#ifndef __TWO_PHASE_FREE_BOUNDARY_MESH_H__
#define __TWO_PHASE_FREE_BOUNDARY_MESH_H__

#include "generic.h"
#include "meshes/rectangular_quadmesh.h"

using namespace oomph;

class FreeBoundaryGeometry : public GeomObject {
	protected:
		Vector<Data *> data_pt;
		unsigned int free_boundary_index = 1;

		bool destroy_geom_data = false;
		TimeStepper *ts_pt;
	public:
		FreeBoundaryGeometry(
			const double x0,
			const double x1,
			const double x2,
			const double y0,
			const double y1,
			TimeStepper *timestepper = new Steady<0>
		) : GeomObject(2,2), ts_pt(timestepper) {
			
			data_pt.resize(1);
			data_pt[0] = new Data(ts_pt, 5);

			// assume impulsive conditions
			for (unsigned int t = 0; t < ts_pt->nprev_values(); t++) {
				data_pt[0]->set_value(t, 0, x0);
				data_pt[0]->set_value(t, 1, x1);
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

		double& x0(const unsigned int &t) const {return *data_pt[0]->value_pt(t, 0);}
		double& x1(const unsigned int &t) const {return *data_pt[0]->value_pt(t, free_boundary_index);}
		double& x2(const unsigned int &t) const {return *data_pt[0]->value_pt(t, 2);}
		double& y0(const unsigned int &t) const {return *data_pt[0]->value_pt(t, 3);}
		double& y1(const unsigned int &t) const {return *data_pt[0]->value_pt(t, 4);}

		double& x0() const {return x0(0);}
		double& x1() const {return x1(0);}
		double& x2() const {return x2(0);}
		double& y0() const {return y0(0);}
		double& y1() const {return y1(0);}

		unsigned int ngeom_data() const {return data_pt.size();}

		Data* geom_data_pt(const unsigned &j) {return data_pt[0];}
		
		double get_interface() {return x1();}
		void set_interface(double &x) {x1() = x;}

		void position(const unsigned int &t, const Vector<double> &zeta, Vector<double> &r) const {
			r[0] = x0(t) + (x2(t) - x0(t)) * zeta[0];
			r[1] = y0(t) + (y1(t) - y0(t)) * zeta[1];

			printf("FreeBoundaryGeometry - position: zeta0: %8.6f zeta1: %8.6f\n", zeta[0], zeta[1]);
		}

		void position(const Vector<double> &zeta, Vector<double> &r) const {
			position(0, zeta, r);
		}
};

class FreeBoundaryElement : public GeneralisedElement, 
							public FreeBoundaryGeometry {
	private:
		double factor;
		unsigned int geometry_index;
		unsigned int flux_index;
		Data *flux_data_pt;

		TimeStepper *ts_pt;

	public:
		FreeBoundaryElement(
			const double x0,
			const double x1,
			const double x2,
			const double y0,
			const double y1,
			const double factor_,
			TimeStepper *timestepper = new Steady<0>
		) : FreeBoundaryGeometry(x0, x1, x2, y0, y1, timestepper), factor(factor_), ts_pt(timestepper) {

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

		void set_flux(const unsigned int &t, const double &flux) {external_data_pt(flux_index)->set_value(t, 0, flux);}
		void set_flux(const double &flux) {
			printf("Setting flux to %8.6f\n", flux);
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
			unsigned ndofs = ndof();

			if (ndofs == 0) return;

			int free_boundary_local_eqn_number = internal_local_eqn(geometry_index, free_boundary_index);

			Data *interface_data_pt = internal_data_pt(geometry_index);
			TimeStepper *interface_ts_pt = interface_data_pt->time_stepper_pt();

			residuals[free_boundary_local_eqn_number] = factor * interface_ts_pt->time_derivative(1, interface_data_pt, free_boundary_index) - external_data_pt(flux_index)->value(0);
			printf("Factor: %8.6f dhdt=%8.6f\n", factor, interface_ts_pt->time_derivative(1, interface_data_pt, free_boundary_index));
			if (compute_jacobian) {
				printf("Computing jacobian\n");
				jacobian(free_boundary_local_eqn_number, free_boundary_local_eqn_number) = factor;//interface_ts_pt->time_derivative(1, interface_data_pt, free_boundary_index);
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
				
				Vector<GeomObject *> geom_object_pt(1);
				geom_object_pt[0] = geometry;

				elem->set_node_update_info(geom_object_pt);
			}
			
			this->setup_boundary_element_info();
		}
};

#endif