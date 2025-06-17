#ifndef __TWO_LAYER_MESH_H__
#define __TWO_LAYER_MESH_H__

#include "generic.h"
#include "meshes/one_d_mesh.h"

using namespace oomph;

template<class EL> class TwoLayerMesh : public virtual OneDMesh<EL> {
	public:
		TwoLayerMesh(
			const unsigned int &nx1, 
			const unsigned int &nx2,
			const double x0,
			const double x1,
			const double x2,
			TimeStepper *ts_pt = &Mesh::Default_TimeStepper
		) : OneDMesh<EL>(nx1+nx2, x0, x2, ts_pt) {
			this->set_nboundary(3);

			Node *node = this->finite_element_pt(nx1)->node_pt(0);

			this->convert_to_boundary_node(node);
			this->add_boundary_node(2, node);

			this->setup_boundary_element_info();
		}
};

#endif
