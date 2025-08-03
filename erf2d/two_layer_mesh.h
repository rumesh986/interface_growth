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

class TwoPhaseDomain : public Domain {
	private:
		const double x0;
		double x1;
		const double x2;
		const double vel;

		
		const unsigned int nx1;
		const unsigned int nx2;
		const unsigned int ny;
		const unsigned int nx;
		const unsigned int nmacro;
		
		Time *time_pt;
	public:
		TwoPhaseDomain(
			const double x0,
			const double x1,
			const double x2,
			const double vel,
			const unsigned int nx1,
			const unsigned int nx2,
			const unsigned int ny,
			Time *time_pt
		) : x0(x0), x1(x1), x2(x2), vel(vel), nx1(nx1), nx2(nx2), ny(ny), nx(nx1+nx2), nmacro(nx*ny), time_pt(time_pt) {
			printf("Creating %u macro elements in twophasedomain\n", nmacro);
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

			double time = time_pt->time(t);

			// interface location
			const double x_int = (vel > 0.0) ? x1 + vel * time : x1;

			const unsigned int yi = macro_i / nx;
			const unsigned int xi = macro_i % nx;

			const double dx1 = ((double)(x_int - x0)) / (double) nx1;
			const double dx2 = ((double)(x2 - x_int)) / (double) nx2;
			const double dy = 1.0 / (double) ny;
			
			bool growing = xi < nx1;

			double trans_x = 0.5 * (zeta[0] + 1.0);
			
			double start = (growing) ? x0 + (((double) xi) * dx1) : x_int + ((((double) xi)-nx1) * dx2);
			double end = (growing) ? x0 + ((((double) xi)+1) * dx1) : x_int + ((((double) xi)-nx1+1) * dx2);
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
					printf("Invalid direction given in TwoPhaseDomain");
					return;
			}
		}

		double get_interface() {
			return x1;
		}

		void get_interface(Vector<double> &x) {
			x[0] = x1;
		}

		void get_interface(const double &t, Vector<double> &x) {
			x[0] = (vel > 0.0) ? x1 + vel * t : x1;
		}

		void set_interface(double x) {
			x1 = x;
		}

		void doc_domain(DocInfo info, char *label) {
			char filename[100];
			ofstream outfile;

			sprintf(filename, "%s/domain-%s_%i.dat", info.directory().c_str(), label, info.number());
			outfile.open(filename);
			output(outfile, 5);
			outfile.close();
		}
};

template<class EL> class RefineableTwoLayer2DMesh : public virtual RefineableRectangularQuadMesh<EL> {
	public:

		TwoPhaseDomain *domain;

		RefineableTwoLayer2DMesh(
			const unsigned int &nx1, 
			const unsigned int &nx2,
			const unsigned int &ny,
			const double x0,
			const double x1,
			const double x2,
			const double y0,
			const double y1,
			TwoPhaseDomain *domain,
			TimeStepper *ts_pt = &Mesh::Default_TimeStepper
		) :	RectangularQuadMesh<EL>(nx1+nx2, ny, x0, x2, y0, y1, ts_pt),
			RefineableRectangularQuadMesh<EL>(nx1+nx2, ny, x0, x2, y0, y1, ts_pt) {

			this->domain = domain;

			// printf("Looping through current elements\n");
			Vector<double> s_fraction(2);

			// for (uint e = 0; e < this->nelement(); e++) {
			// 	FiniteElement *elem = this->finite_element_pt(e);
			// 	for (uint n = 0; n < elem->nnode(); n++) {
			// 		printf("e=%u n=%u ", e, n);

			// 		elem->local_fraction_of_node(n, s_fraction);

			// 		for (uint i = 0; i < elem->nnode_1d(); i++) {
			// 			printf("x%u=%10.8f s_frac=%10.8f ", i, elem->node_pt(n)->x(i), s_fraction[i]);
			// 		}
			// 		printf("\n");
			// 	}
			// }

			// modified from fish_mesh.template.cc
			Vector<double> s(2);
    		Vector<double> r(2);

			printf("CONFIG: \n");
			printf("\tnx1 = %u\n", nx1);
			printf("\tnx2 = %u\n", nx2);
			printf("\tnx = %u\n", nx1+nx2);
			printf("\tny = %u\n", ny);

			for (uint yi = 0; yi < ny; yi++) {
				for (uint xi = 0; xi < nx1+nx2; xi++) {
					uint e_loc = xi + yi * (nx1+nx2);
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

			this->set_nboundary(5);
			// printf("Back to regular bit for adding boundary\n");

			for (uint e = 0; e < ny; e++) {
				EL *elem = dynamic_cast<EL *>(this->finite_element_pt(nx1*(1+e) + nx2*e));
				uint nnode_1d = elem->nnode_1d();
				for (uint n = 0; n < nnode_1d; n++) {
					Node *node = elem->node_pt((nnode_1d-1)*n);

					this->convert_to_boundary_node(node);
					this->add_boundary_node(4, node);
				}
			}

			this->setup_boundary_element_info();
		}

};

#endif
