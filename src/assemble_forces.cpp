#include <assemble_forces.h>
#include <dV_cloth_gravity_dq.h>
#include <dV_membrane_corotational_dq.h>
#include <mass_matrix_mesh.h>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                     double mu, double lambda) { 
    f.resize(q.size());
    f.setZero();
    // all triangles
    for (int r = 0; r < F.rows(); r ++) {
        Eigen::RowVector3i element = F.row(r);
        Eigen::Vector9d dV;
        Eigen::Matrix<double,1,9> tmp_row;
        tmp_row = dX.row(r);
        dV_membrane_corotational_dq(dV, q, Eigen::Map<const Eigen::Matrix3d>(tmp_row.data()), V, element, a0(r), mu, lambda);
        // f = - \partial{V}{q} = - \sum \partial{V_r}{q}
        //   = - \sum (E^T {dV})_r
        //  E^T(3 * idx[i] + j, 3i + j) = 1 for i = 0..2, j = 0..2
        //  dV(x) for x = 0..8
        //  (E^T dV)(3 * idx[i] + j) += -dV(3i + j) for i = 0..2, j = 0..2
        for (int i = 0; i < 3; i ++) {
            for (int j = 0; j < 3; j ++) {
                f(3*element(i)+j) -= dV(3*i+j);
            }
        }
    }
};
