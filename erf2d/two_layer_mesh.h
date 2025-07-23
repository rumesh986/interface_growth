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


// have two macro elements, 0 for growing and 1 for shrinking
// keeps it to one domain and macromap should take care of the rest?

class TwoPhaseDomain : public Domain {
	private:
		const double x0;
		const double x1;
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

		void macro_element_boundary(const unsigned &t, const unsigned &macro_i, const unsigned &dir_i, const Vector<double> &zeta, Vector<double> &r) {
			if (macro_i >= nmacro) {
				printf("Invalid Macro index given in TwoPhaseDomain");
				return;
			}

			using namespace QuadTreeNames;
			
			static double min_zeta0 = 0.0;
			static double max_zeta0 = 0.0;
			static double min_zeta1 = 0.0;
			static double max_zeta1 = 0.0;

			if (zeta[0] < min_zeta0)	min_zeta0 = zeta[0];
			if (zeta[0] > max_zeta0)	max_zeta0 = zeta[0];
			if (zeta[1] < min_zeta1)	min_zeta1 = zeta[1];
			if (zeta[1] > max_zeta1)	max_zeta1 = zeta[1];

			double time = time_pt->time(t);

			// interface location
			const double x_int = x1 + vel * time;

			const unsigned int yi = macro_i / nx;
			const unsigned int xi = macro_i % nx;

			const double dx1 = ((double)(x_int - x0)) / (double) nx1;
			const double dx2 = ((double)(x2 - x_int)) / (double) nx2;
			
			// printf("e=%u x=%u y=%u dx1=%8.6f dx2=%8.6f ", macro_i, xi, yi, dx1, dx2);

			bool growing = xi < nx1;

			double trans_x = 0.5 * (zeta[0] + 1.0);

			// printf("MACRO_ELEM_BOUNDS minzeta0=%8.6f maxzeta0=%8.6f minzeta1=%8.6f maxzeta1=%8.6f  \n", min_zeta0, max_zeta0, min_zeta1, max_zeta1);
			// printf("MACRO_ELEM_BOUNDS zeta0=%8.6f zeta1=%8.6f\n", zeta[0], zeta[1]);

			
			// double start = (macro_i == 0) ? x0 : x1 + vel * time;
			// double end = (macro_i == 0) ? x1 + vel * time : x2;
			// double x = trans_x * (start - end) + end;
			
			double start = (growing) ? x0 + (xi * dx1) : x_int + ((xi-nx1) * dx2);
			double end = (growing) ? x0 + ((xi+1) * dx1) : x_int + ((xi-nx1+1) * dx2);
			double x = trans_x  * (end - start) + start;

			// printf("start=%8.6f end=%8.6f zeta=%8.6f x=%8.6f x_int=%8.6f\n", start, end, zeta[0], x, x_int);

			r[1] = zeta[1];
			switch (dir_i) {
				case N:	r[0] = x;		break;
				case E:	r[0] = end;		break;
				case S:	r[0] = x;		break;
				case W:	r[0] = start;	break;
				default: 
					printf("Invalid direction given in TwoPhaseDomain");
					return;
			}
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
		) : RefineableRectangularQuadMesh<EL>(nx1+nx2, ny, x0, x2, y0, y1, ts_pt),
			RectangularQuadMesh<EL>(nx1+nx2, ny, x0, x2, y0, y1, ts_pt) {

			this->domain = domain;

			printf("building mesh i guess \n");
			this->build_mesh(ts_pt);

			printf("Looping through current elements\n");
			Vector<double> s_fraction(2);

			for (uint e = 0; e < this->nelement(); e++) {
				FiniteElement *elem = this->finite_element_pt(e);
				for (uint n = 0; n < elem->nnode(); n++) {
					printf("e=%u n=%u ", e, n);

					elem->local_fraction_of_node(n, s_fraction);

					for (uint i = 0; i < elem->nnode_1d(); i++) {
						printf("x%u=%10.8f s_frac=%10.8f ", i, elem->node_pt(n)->x(i), s_fraction[i]);
					}
					printf("\n");
				}
			}

			printf("In refineabletwolayer2dmesh constructor - new bit\n");

			
			// modified from fish_mesh.template.cc
			// uint nelem = this->nelement();
			Vector<double> s(2);
			// Vector<double> s_fraction(2);
    		Vector<double> r(2);

			printf("CONFIG: \n");
			printf("\tnx1 = %u\n", nx1);
			printf("\tnx2 = %u\n", nx2);
			printf("\tnx = %u\n", nx1+nx2);
			printf("\tny = %u\n", ny);

			// printf("In refineabletwolayer2dmesh constructor - starting loop\n");
			for (uint yi = 0; yi < ny; yi++) {
				// printf("In refineabletwolayer2dmesh constructor - Working on left block\n");
				for (uint xi = 0; xi < nx1; xi++) {
					uint e_loc = xi + yi * (nx1+nx2);
					FiniteElement *elem = this->finite_element_pt(e_loc);
					elem->set_macro_elem_pt(domain->macro_element_pt(e_loc));
					// printf("In refineabletwolayer2dmesh constructor - yi=%u, xi=%u, e_loc=%u\n", yi, xi, e_loc);
					for (uint n = 0; n < elem->nnode(); n++) {
						Node *node = elem->node_pt(n);
						elem->local_coordinate_of_node(n, s);

						// s_fraction[0] = (node->x(0) - x0) / (x1 - x0);
						// s_fraction[1] = (node->x(1) - y0) / (y1 - y0);

						// s[0] = -1.0 + 2.0 * s_fraction[0];
						// s[1] = -1.0 + 2.0 * s_fraction[1];

						// s[0] = s_fraction[0];
						// s[1] = s_fraction[1];

						// printf("Node %u: sfrac0=%8.6f sfrac1=%8.6f s0=%8.6f s1=&8.6f\n", n, s_fraction[0], s_fraction[1], s[0], s[1]);

						domain->macro_element_pt(e_loc)->macro_map(s, r);

						node->x(0) = r[0];
						// node->x(1) = r[1];

						// printf("node=%u\ts_frac: ", n);
						// for (auto i: s_fraction)	cout << i << " ";
						// printf("\ts: ");
						// for (auto i: s)	cout << i << " ";
						// printf("\tr: ");
						// for (auto i: r)	cout << i << " ";
						// cout << endl;
					}
				}

				// printf("In refineabletwolayer2dmesh constructor - Working on right block\n");
				for (uint xi = nx1; xi < nx1+nx2; xi++) {
					uint e_loc = xi + yi * (nx1+nx2);
					FiniteElement *elem = this->finite_element_pt(e_loc);
					elem->set_macro_elem_pt(domain->macro_element_pt(e_loc));
					// printf("In refineabletwolayer2dmesh constructor - yi=%u, xi=%u, e_loc=%u\n", yi, xi, e_loc);
					for (uint n = 0; n < elem->nnode(); n++) {
						Node *node = elem->node_pt(n);
						elem->local_fraction_of_node(n, s);

						// s_fraction[0] = (node->x(0) - x1) / (x2 - x1);
						// s_fraction[1] = (node->x(1) - y0) / (y1 - y0);

						// s[0] = -1.0 + 2.0 * s_fraction[0];
						// s[1] = -1.0 + 2.0 * s_fraction[1];

						domain->macro_element_pt(e_loc)->macro_map(s, r);

						node->x(0) = r[0];
						// node->x(1) = r[1];

						// printf("node=%u\ts_frac: ", n);
						// for (auto i: s_fraction)	cout << i << " ";
						// printf("\ts: ");
						// for (auto i: s)	cout << i << " ";
						// printf("\tr: ");
						// for (auto i: r)	cout << i << " ";
						// cout << endl;
					}
				}
			}

			this->set_nboundary(5);
			printf("Back to regular bit for adding boundary\n");
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

		void build_mesh(TimeStepper *time_stepper_pt) {
			MeshChecker::assert_geometric_element<QElementGeometricBase, EL>(2);
		}
};

#endif
