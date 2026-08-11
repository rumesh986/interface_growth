#ifndef __UNSTEADY_HEAT_MOVING_H__
#define __UNSTEADY_HEAT_MOVING_H__

#include "generic.h"
#include "unsteady_heat.h"

using namespace oomph;

template<unsigned int DIM>
class UnsteadyHeatMovingEquations : public virtual UnsteadyHeatEquations<DIM> {
	public:
		UnsteadyHeatMovingEquations() : UnsteadyHeatEquations<DIM>() {
			V_pt = &Default_v_parameter;
			V_direction_pt = &Default_v_direction;
		};

		const double &v() const {
			return *V_pt;
		}

		double*& v_pt() {
			return V_pt;
		}

		const unsigned int &v_direction() const {
			return *V_direction_pt;
		}

		unsigned int*& v_direction_pt() {
			return V_direction_pt;
		}

		void fill_in_contribution_to_residuals(Vector<double> &residuals) {
			// fill_in_generic_residual_contribution_ust_heat(residuals, GeneralisedElement::Dummy_matrix, false);
			fill_in_generic_residual_contribution_ust_heat(residuals, GeneralisedElement::Dummy_matrix, false);

		}

		void fill_in_contribution_to_jacobian(Vector<double> &residuals, DenseMatrix<double> &jacobian) {
			// fill_in_generic_residual_contribution_ust_heat(residuals, jacobian, true);
			fill_in_generic_residual_contribution_ust_heat(residuals, jacobian, true);

		}
	
	protected:
		double *V_pt;
		unsigned int *V_direction_pt;
	
		void fill_in_generic_residual_contribution_ust_heat(Vector<double>& residuals, DenseMatrix<double>& jacobian, bool flag) {
    // Find out how many nodes there are
    unsigned n_node = this->nnode();

    // Get continuous time from timestepper of first node
    double time = this->node_pt(0)->time_stepper_pt()->time_pt()->time();

    // Find the index at which the variable is stored
    unsigned u_nodal_index = this->u_index_ust_heat();

    // Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node, DIM), dtestdx(n_node, DIM);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    // Set the Vector to hold local coordinates
    Vector<double> s(DIM);

    // Get Alpha and beta parameters number
    double alpha_local = this->alpha();
    double beta_local = this->beta();

    // Integers to hold the local equation and unknowns
    int local_eqn = 0, local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < DIM; i++) s[i] = this->integral_pt()->knot(ipt, i);

      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Call the derivatives of the shape and test functions
      double J = this->dshape_and_dtest_eulerian_at_knot_ust_heat(
        ipt, psi, dpsidx, test, dtestdx);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Allocate memory for local quantities and initialise to zero
      double interpolated_u = 0.0;
      double dudt = 0.0;
      Vector<double> interpolated_x(DIM, 0.0);
      Vector<double> interpolated_dudx(DIM, 0.0);
      Vector<double> mesh_velocity(DIM, 0.0);

      // Calculate function value and derivatives:
      // Loop over nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Calculate the value at the nodes
        double u_value = this->raw_nodal_value(l, u_nodal_index);
        interpolated_u += u_value * psi(l);
        dudt += this->du_dt_ust_heat(l) * psi(l);
        // Loop over directions
        for (unsigned j = 0; j < DIM; j++)
        {
          interpolated_x[j] += this->raw_nodal_position(l, j) * psi(l);
          interpolated_dudx[j] += u_value * dpsidx(l, j);
        }
      }

      // Mesh velocity?
      if (!this->ALE_is_disabled)
      {
        for (unsigned l = 0; l < n_node; l++)
        {
          for (unsigned j = 0; j < DIM; j++)
          {
            mesh_velocity[j] += this->raw_dnodal_position_dt(l, j) * psi(l);
          }
        }
      }

      // Get source function
      //-------------------
      double source = 0.0;
      this->get_source_ust_heat(time, ipt, interpolated_x, source);

      // Assemble residuals and Jacobian
      //--------------------------------

      // Loop over the test functions
      for (unsigned l = 0; l < n_node; l++)
      {
        local_eqn = this->nodal_local_eqn(l, u_nodal_index);
        /*IF it's not a boundary condition*/
        if (local_eqn >= 0)
        {
          // Add body force/source term and time derivative
          residuals[local_eqn] += (alpha_local * (dudt - v() * interpolated_dudx[v_direction()]) + source) * test(l) * W;

          // The mesh velocity bit
          if (!this->ALE_is_disabled)
          {
            for (unsigned k = 0; k < DIM; k++)
            {
              residuals[local_eqn] -= alpha_local * mesh_velocity[k] *
                                      interpolated_dudx[k] * test(l) * W;
            }
          }

          // Laplace operator
          for (unsigned k = 0; k < DIM; k++)
          {
            residuals[local_eqn] +=
              beta_local * interpolated_dudx[k] * dtestdx(l, k) * W;
          }


          // Calculate the jacobian
          //-----------------------
          if (flag)
          {
            // Loop over the velocity shape functions again
            for (unsigned l2 = 0; l2 < n_node; l2++)
            {
              local_unknown = this->nodal_local_eqn(l2, u_nodal_index);
              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                // Mass matrix
                jacobian(local_eqn, local_unknown) +=
                  alpha_local * test(l) * (psi(l2) *
                  this->node_pt(l2)->time_stepper_pt()->weight(1, 0) - v() * dpsidx(l2, v_direction())) * W;

                // Laplace operator & mesh velocity bit
                for (unsigned i = 0; i < DIM; i++)
                {
                  double tmp = beta_local * dtestdx(l, i);
                  if (!this->ALE_is_disabled)
                    tmp -= alpha_local * mesh_velocity[i] * test(l);
                  jacobian(local_eqn, local_unknown) += dpsidx(l2, i) * tmp * W;
                }
              }
            }
          }
        }
      }


    } // End of loop over integration points
  }
	
	private:
		static double Default_v_parameter;
		static unsigned int Default_v_direction;
};

template<unsigned int DIM>
double UnsteadyHeatMovingEquations<DIM>::Default_v_parameter = 0.0;

template<unsigned int DIM>
unsigned int UnsteadyHeatMovingEquations<DIM>::Default_v_direction = 0;

template<unsigned int DIM, unsigned int NNODE_1D>
class QUnsteadyHeatMovingElement :  
public virtual QElement<DIM, NNODE_1D>,			
public virtual QUnsteadyHeatElement<DIM, NNODE_1D>,
									public virtual UnsteadyHeatMovingEquations<DIM> {

	public:
		QUnsteadyHeatMovingElement() : 	QUnsteadyHeatElement<DIM, NNODE_1D>(),
										UnsteadyHeatMovingEquations<DIM>() {};
};


// Flux elements (Face Elements) dont have the required FiniteElement::shape(...) overrides
// This bit is required to ensure that Flux Elements inherit from QElement of a dimension lower than bulk
// Will probably need to do this again for TElement in the future
template<unsigned DIM, unsigned NNODE_1D>
class FaceGeometry<QUnsteadyHeatMovingElement<DIM, NNODE_1D>> : public virtual QElement<DIM - 1, NNODE_1D> {
	public:
		/// Constructor: Call the constructor for the
		/// appropriate lower-dimensional QElement
		FaceGeometry() : QElement<DIM - 1, NNODE_1D>() {}
};

#endif