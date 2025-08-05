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

			this->setup_boundary_element_info();
		}
};

// template<class EL> class RefineableThreePhase2DMesh 
// 	: public virtual RefineableQuadMesh<EL>,
// 	  public virtual ThreePhase2DMesh<EL> {
// 		public:
// 			TwoPhaseDomain *domain;

// 			RefineableThreePhase2DMesh(
// 				const unsigned int &nx1,
// 				const unsigned int &nx2,
// 				const unsigned int &nx3,
// 				const unsigned int &ny,
// 				const double x0,
// 				const double x1,
// 				const double x2,
// 				const double x3,
// 				const double y0,
// 				const double y1,
// 				const double vel,
// 				Time *time_pt,
// 				TimeStepper *ts_pt = &Mesh::Default_TimeStepper
// 			) : RectangularQuadMesh(nx1+nx2+nx3, ny, x0, x3, y0, y1, ts_pt),
// 				RefineableRectangularQuadMesh(nx1+nx2+nx3, ny, x0, x3, y0, y1, ts_pt),
// 				ThreePhase2DMesh(nx1, nx2, nx3, ny, x0, x1, x2, x3, y0, y1, ts_pt) {

// 					domain = new TwoPhaseDomain(x0, x2, x3, vel, nx1+nx2, nx3, ny, time_pt);

// 					// modified from two_layer_mesh.cc (refineabletwolayermesh)
// 					Vector<double> s_fraction(2);
// 					Vector<double> s(2);
// 					Vector<double> r(2);

// 					printf("CONFIG: \n");
// 					printf("\tnx1 = %u\n", nx1);
// 					printf("\tnx2 = %u\n", nx2);
// 					printf("\tnx3 = %u\n", nx3);
// 					printf("\tnx = %u\n", nx1+nx2+nx3);
// 					printf("\tny = %u\n", ny);

// 					for (uint yi = 0; yi < ny; yi++) {
// 						for (uint xi = 0; xi < nx1+nx2+nx3; xi++) {
// 							uint e_loc = xi + yi * (nx1+nx2+nx3);
// 							FiniteElement *elem = this->finite_element_pt(e_loc);
// 							elem->set_macro_elem_pt(domain->macro_element_pt(e_loc));
// 							for (uint n = 0; n < elem->nnode(); n++) {
// 								Node *node = elem->node_pt(n);
// 								elem->local_fraction_of_node(n, s_fraction);
								
// 								s[0] = -1.0 + 2.0 * s_fraction[0];
// 								s[1] = -1.0 + 2.0 * s_fraction[1];
								
// 								domain->macro_element_pt(e_loc)->macro_map(s, r);
								
// 								node->x(0) = r[0];
// 								node->x(1) = r[1];
// 							}
// 						}
// 					}
// 				}

// };


#endif