#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <EigenTypes.h>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.  
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename FORCE, typename STIFFNESS> 
inline void linearly_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            const Eigen::SparseMatrixd &mass,  FORCE &force, STIFFNESS &stiffness, 
                            Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness) {
    // (M - dt^2 K) qdot<t+1> = M qdot<t> + dt f(q<t>)
    //                 q<t+1> = q<t> + dt qdot<t+1>
    //  A = M - dt^2 K
    //  b = M qdot<t> + dt f(q<t>)
    //  Solve A qdot = b
    //  Then q = q + dt * qdot

    stiffness(tmp_stiffness, q, qdot);
    force(tmp_force, q, qdot);
    Eigen::SparseMatrixd A = mass - pow(dt, 2) * tmp_stiffness;
    Eigen::VectorXd b = mass * qdot + dt * tmp_force;

    // Use Simplicial LDLT solver to solve the SPD system
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        std::cerr << "ERROR: Solving Ax = b failed!" << std::endl;
        return;
    }
    qdot = solver.solve(b);
    q = q + dt * qdot; 
}
