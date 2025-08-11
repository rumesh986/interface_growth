#ifndef __THREE_LAYER_MESH_H__
#define __THREE_LAYER_MESH_H_

#include "generic.h"
#include "meshes/rectangular_quadmesh.h"
#include "two_layer_mesh.h"

using namespace oomph;

template<class EL> class ThreePhase2DMesh : public virtual RectangularQuadMesh<EL> {
	public:
		ThreePhase2DMesh(
			const unsigned int &nx1,
			const unsigned int &nx2,
			const unsigned int &nx3,
			const unsigned int &ny,
			const double x0,
			const double x1,
			const double x2,
			const double x3,
			const double y0,
			const double y1,
			TimeStepper *ts_pt = &Mesh::Default_TimeStepper
		) : RectangularQuadMesh<EL>(nx1+nx2+nx3, ny, x0, x3, y0, y1, ts_pt) {

			printf("In regular mesh constructor\n");

			this->set_nboundary(6);

			uint nx = nx1 + nx2 + nx3;

			// set boundary 5 (between phase 1 and 2)
			for (uint e = 0; e < ny; e++) {
				FiniteElement *elem = this->finite_element_pt(nx*e + nx1);
				uint nnode = elem->nnode_1d();

				for (uint n = 0; n < nnode; n++) {
					Node *node = elem->node_pt(nnode*n);

					this->convert_to_boundary_node(node);
					this->add_boundary_node(4, node);
				}
			}

			// set boundary 6 (between phase 2 and 3)
			for (uint e = 0; e < ny; e++) {
				FiniteElement *elem = this->finite_element_pt(nx*e + nx1 + nx2);
				uint nnode = elem->nnode_1d();

				for (uint n = 0; n < nnode; n++) {
					Node *node = elem->node_pt(nnode*n);

					this->convert_to_boundary_node(node);
					this->add_boundary_node(5, node);
				}
			}
		}
};

class ThreePhaseDomain : public Domain {
	private:
		const double x0;
		const double x1;
		double x2;
		const double x3;

		const unsigned int nx1;
		const unsigned int nx2;
		const unsigned int nx3;
		const unsigned int nx;
		const unsigned int ny;
		const unsigned int nmacro;

		Time *time_pt;
	
	public:
		ThreePhaseDomain(
			const double x0,
			const double x1,
			const double x2,
			const double x3,
			const unsigned int nx1,
			const unsigned int nx2,
			const unsigned int nx3,
			const unsigned int ny,
			Time *time_pt
		) : x0(x0), x1(x1), x2(x2), x3(x3), nx1(nx1), nx2(nx2), nx3(nx3), nx(nx1+nx2+nx3), ny(ny), nmacro(nx*ny), time_pt(time_pt) {
			Macro_element_pt.resize(nmacro);

			for (unsigned int i = 0; i < nmacro; i++)
				Macro_element_pt[i] = new QMacroElement<2>(this, i);
		}

		// zeta is a 1D vector with values between -1 and 1
		void macro_element_boundary(const unsigned &t, const unsigned &macro_i, const unsigned &dir_i, const Vector<double> &zeta, Vector<double> &r) {
			if (macro_i >= nmacro) {
				printf("Invalid Macro index given in TwoPhaseDomain");
				return;
			}

			using namespace QuadTreeNames;

			const unsigned int yi = macro_i / nx;
			const unsigned int xi = macro_i % nx;

			const double dx1 = ((double)(x1 - x0)) / (double) nx1;
			const double dx2 = ((double)(x2 - x1)) / (double) nx2;
			const double dx3 = ((double)(x3 - x2)) / (double) nx3;
			const double dy = 1.0 / (double) ny;

			double start, end;

			if (xi < nx1) {
				// substrate phase
				start = x0 + ((double) (xi)) * dx1;
				end = x0 + ((double) (xi + 1)) * dx1;
			} else if (xi < nx1+nx2) {
				// growing solid phase
				start = x1 + ((double) (xi - nx1)) * dx2;
				end = x1 + ((double) (xi - nx1 + 1)) * dx2;
			} else {
				// shrinking liquid phase
				start = x2 + ((double) (xi - nx1 - nx2)) * dx3;
				end = x2 + ((double) (xi - nx1 - nx2 + 1)) * dx3;
			}
			
			double trans_x = 0.5 * (zeta[0] + 1.0);
			double x = trans_x  * (end - start) + start;

			switch (dir_i) {
				case N:	
					r[0] = x;
					r[1] = (((double) yi) + 1) * dy;
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
					printf("Invalid direction given in ThreePhaseDomain");
					return;
			}
		}

		double get_interface() {
			return x2;
		}

		void set_interface(double x) {
			x2 = x;
		} 
};

template<class EL> class RefineableThreePhase2DMesh 
	: public virtual RefineableRectangularQuadMesh<EL>,
	  public virtual ThreePhase2DMesh<EL> {
		public:
			ThreePhaseDomain *domain;

			RefineableThreePhase2DMesh(
				const unsigned int &nx1,
				const unsigned int &nx2,
				const unsigned int &nx3,
				const unsigned int &ny,
				const double x0,
				const double x1,
				const double x2,
				const double x3,
				const double y0,
				const double y1,
				const double vel,
				Time *time_pt,
				TimeStepper *ts_pt = &Mesh::Default_TimeStepper
			) : 
				RectangularQuadMesh<EL>(nx1+nx2+nx3, ny, x0, x3, y0, y1, ts_pt),
				RefineableRectangularQuadMesh<EL>(nx1+nx2+nx3, ny, x0, x3, y0, y1, ts_pt),
				ThreePhase2DMesh<EL>(nx1, nx2, nx3, ny, x0, x1, x2, x3, y0, y1, ts_pt) {

					printf("In refineable mesh constructor\n");
					domain = new ThreePhaseDomain(x0, x1, x2, x3, nx1, nx2, nx3, ny, time_pt);

					// modified from two_layer_mesh.cc (refineabletwolayermesh)
					Vector<double> s_fraction(2);
					Vector<double> s(2);
					Vector<double> r(2);

					printf("CONFIG: \n");
					printf("\tnx1 = %u\n", nx1);
					printf("\tnx2 = %u\n", nx2);
					printf("\tnx3 = %u\n", nx3);
					printf("\tnx = %u\n", nx1+nx2+nx3);
					printf("\tny = %u\n", ny);

					for (uint yi = 0; yi < ny; yi++) {
						for (uint xi = 0; xi < nx1+nx2+nx3; xi++) {
							uint e_loc = xi + yi * (nx1+nx2+nx3);
							FiniteElement *elem = this->finite_element_pt(e_loc);
							elem->set_macro_elem_pt(domain->macro_element_pt(e_loc));
							for (uint n = 0; n < elem->nnode(); n++) {
								Node *node = elem->node_pt(n);
								elem->local_fraction_of_node(n, s_fraction);
								
								s[0] = -1.0 + 2.0 * s_fraction[0];
								s[1] = -1.0 + 2.0 * s_fraction[1];
								
								domain->macro_element_pt(e_loc)->macro_map(s, r);
								
								node->x(0) = r[0];
								node->x(1) = r[1];
							}
						}
					}
				}
};

#endif