#ifndef __TWO_LAYER_MESH_H__
#define __TWO_LAYER_MESH_H__

#include "generic.h"
#include "meshes/rectangular_quadmesh.h"

using namespace oomph;

template<class EL> class TwoLayer2DMesh : public virtual RectangularQuadMesh<EL> {
	public:
		TwoLayer2DMesh(
			const unsigned int &nx1, 
			const unsigned int &nx2,
			const unsigned int &ny,
			const double x0,
			const double x1,
			const double x2,
			const double y0,
			const double y1,
			TimeStepper *ts_pt = &Mesh::Default_TimeStepper
		) : RectangularQuadMesh<EL>(nx1+nx2, ny, x0, x2, y0, y1, ts_pt) {
			this->set_nboundary(5);

			for (int e = 0; e < ny; e++) {
				EL *elem = dynamic_cast<EL *>(this->finite_element_pt(nx1*(1+e) + nx2*e));
				int nnode = elem->nnode_1d();
				for (int n = 0; n < nnode; n++) {
					Node *node = elem->node_pt((nnode-1)*n);

					this->convert_to_boundary_node(node);
					this->add_boundary_node(4, node);
				}
			}

			this->setup_boundary_element_info();
		}
};

#endif
