#include "generic.h"
#include "solid.h"
#include "linear_elasticity.h"

#include "meshes/triangle_mesh.h"
#include "meshes/rectangular_quadmesh.h"

#include "unsteady_heat.h"

using namespace std;
using namespace oomph;
using namespace MathematicalConstants;

namespace ExactSolnForUnsteadyHeat {
    double K = 10;
    double phi = 1.0;

    void get_exact_u(const double &time, const Vector<double> &x, Vector<double> &u) {
        double zeta = cos(phi) * x[0] + sin(phi) * x[1];
        u[0] = exp(-K*time)*sin(zeta*sqrt(K));
    }

    void get_source(const double &time, const Vector<double> &x, double &source) {
        source = 0.0;
    }
}

template<class ELEMENT> class UnsteadyHeatProblem : public Problem {
    public:
        UnsteadyHeatProblem(UnsteadyHeatEquations<2>::UnsteadyHeatSourceFctPt source_fct_ptr);
        ~UnsteadyHeatProblem(){};

        void actions_after_newton_solve() {}
        void actions_before_newton_solce() {}
        void actions_after_implicit_timestep() {}
        void actions_before_implicit_timestep();

        void set_initial_condition();

        void doc_solution(DocInfo &doc_info, ofstream &trace_file);

    private:
        UnsteadyHeatEquations<2>::UnsteadyHeatSourceFctPt Source_fct_ptr;

        Node * Control_node_ptr;
};

template<class ELEMENT> UnsteadyHeatProblem<ELEMENT>::UnsteadyHeatProblem(UnsteadyHeatEquations<2> :: UnsteadyHeatSourceFctPt source_fct_ptr) :
    Source_fct_ptr(source_fct_ptr) {
    
    add_time_stepper_pt(new BDF<2>);

    unsigned int nx = 5;
    unsigned int ny = 5;

    double lx = 1.0;
    double ly = 1.0;

    mesh_pt() = new RectangularQuadMesh<ELEMENT>(nx, ny, lx, ly, time_stepper_pt());

    unsigned n_bound = mesh_pt()->nboundary();
    for (unsigned b = 0; b < n_bound; b++) {
        unsigned n_node = mesh_pt()->nboundary_node(b);
        for (unsigned n = 0; n < n_node; n++) {
            mesh_pt()->boundary_node_pt(b,n)->pin(0);
        }
    }

    unsigned n_element = mesh_pt()->nelement();
    for(unsigned i = 0; i < n_element; i++) {
        ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
        el_pt->source_fct_pt() = Source_fct_ptr;
    }

    unsigned control_el = unsigned(n_element/2);
    Control_node_ptr = mesh_pt()->finite_element_pt(control_el)->node_pt(0);

    cout << "Number of equations: " << assign_eqn_numbers() << endl;
}

template<class ELEMENT> void UnsteadyHeatProblem<ELEMENT>::actions_before_implicit_timestep() {
    double time = time_pt()->time();

    unsigned num_bound = mesh_pt()->nboundary();
    for (unsigned ibound = 0; ibound < num_bound; ibound++) {
        unsigned num_node = mesh_pt()->nboundary_node(ibound);
        for (unsigned inod = 0; inod < num_node; inod++) {
            Node *nod_ptr = mesh_pt()->boundary_node_pt(ibound, inod);
            Vector<double> u(1);
            Vector<double> x(2);

            x[0] = nod_ptr->x(0);
            x[1] = nod_ptr->x(1);

            ExactSolnForUnsteadyHeat::get_exact_u(time, x, u);
            nod_ptr->set_value(0, u[0]);
        }
    }
}

template<class ELEMENT> void UnsteadyHeatProblem<ELEMENT>::set_initial_condition() {
    double backed_up_time = time_pt()->time();

    Vector<double> soln(1);
    // double soln;
    Vector<double> x(2);

    unsigned num_node = mesh_pt()->nnode();

    int nprev_steps = time_stepper_pt()->nprev_values();
    Vector<double> prev_time(nprev_steps+1);

    for (int t = nprev_steps; t >=0; t--) {
        prev_time[t] = time_pt()->time(unsigned(t));
    }

    for (int t = nprev_steps; t >= 0; t--) {
        double time = prev_time[t];
        cout << "setting IC at time = " << time << endl;

        for (unsigned n = 0; n < num_node; n++) {
            x[0] = mesh_pt()->node_pt(n)->x(0);
            x[1] = mesh_pt()->node_pt(n)->x(1);

            ExactSolnForUnsteadyHeat::get_exact_u(time, x, soln);

            mesh_pt()->node_pt(n)->set_value(t, 0, soln[0]);

            for (unsigned i = 0; i < 2; i++) {
                mesh_pt()->node_pt(n)->x(t, i) = x[i];
            }
        }
    }

    time_pt()->time() = backed_up_time;
}

template<class ELEMENT> void UnsteadyHeatProblem<ELEMENT>::doc_solution(DocInfo &doc_info, ofstream &trace_file) {
    ofstream some_file;
    char filename[100];

    unsigned npts;
    npts = 5;

    cout << endl;
    cout << "===========================================" << endl;
    cout << "docing solution for t = " << time_pt()->time() << endl;
    cout << "===========================================" << endl;
    
    sprintf(filename, "%s/soln%i.dat", doc_info.directory().c_str(), doc_info.number());
    some_file.open(filename);
    mesh_pt()-> output(some_file, npts);

    some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = " << time_pt()->time() << "\"";
    some_file << "GEOMETRY X=2.5,Y=98,T=LINE,C=BLUE,LT=0.4" << endl;
    some_file << "1" << endl;
    some_file << "2" << endl;
    some_file << " 0 0" << endl;
    some_file << time_pt()->time()*20.0 << " 0" << endl;
    some_file.close();

    double error, norm;
    sprintf(filename, "%s/error%i.dat", doc_info.directory().c_str(), doc_info.number());

    some_file.open(filename);
    mesh_pt()->compute_error(some_file, ExactSolnForUnsteadyHeat::get_exact_u, time_pt()->time(), error, norm); 
    some_file.close();
    
    // Doc solution and error
    //-----------------------
    cout << "error: " << error << std::endl; 
    cout << "norm : " << norm << std::endl << std::endl;
    
    // Get exact solution at control node
    Vector<double> x_ctrl(2);
    x_ctrl[0]=Control_node_ptr->x(0);
    x_ctrl[1]=Control_node_ptr->x(1);
    Vector<double> u_exact(1);
    ExactSolnForUnsteadyHeat::get_exact_u(time_pt()->time(),x_ctrl,u_exact);
    trace_file << time_pt()->time() << " " << Control_node_ptr->value(0) << " " << u_exact[0] << " " << error   << " " << norm << endl;
}

int main() {
    UnsteadyHeatProblem<QUnsteadyHeatElement<2,4>> problem(&ExactSolnForUnsteadyHeat::get_source);

    DocInfo doc_info;
    doc_info.set_directory("RESLT");
    doc_info.number() = 0;

    ofstream trace_file;
    char filename[100];
    sprintf(filename, "%s/trace.dat", doc_info.directory().c_str());
    trace_file.open(filename);

    trace_file << "VARIABLES=\"times\", \"u<SUB>FE</SUB>\",\"u<SUB>exact</SUB>\",\"norm of error\",\"norm of solutions\"" << endl;

    double t_max = 0.5;
    double dt = 0.01;

    problem.initialise_dt(dt);
    problem.set_initial_condition();
    problem.doc_solution(doc_info, trace_file);
    doc_info.number()++;

    unsigned nstep = unsigned(t_max/dt);
    for (unsigned istep = 0; istep < nstep; istep++) {
        cout << " Timestep " << istep << endl;
        problem.unsteady_newton_solve(dt);
        problem.doc_solution(doc_info, trace_file);
        doc_info.number()++;
    }

    trace_file.close();
};