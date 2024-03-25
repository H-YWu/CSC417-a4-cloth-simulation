#include <dV_cloth_gravity_dq.h>

void dV_cloth_gravity_dq(Eigen::VectorXd &fg, Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::Vector3d> g) {
    int n = M.rows() / 3;
    Eigen::VectorXd G = g.replicate(n, 1);
    fg = M * G;
}
