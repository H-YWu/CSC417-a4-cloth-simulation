#include <dF_cloth_triangle_dq.h>
#include <matrix_B_triangle.h>
#include <deformation_gradient.h>

void skew_matrix(Eigen::Matrix3d &skew, Eigen::Vector3d v) {
    skew << 0.0, -v(2), v(1),
            v(2), 0.0, -v(0),
            -v(1), v(0), 0.0;
}

void dF_cloth_triangle_dq(Eigen::Matrix99d &dFdq, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX,
                        Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element) {
    //  F = \sum_i x_i \diff{phi_i}{X} + n * N^T
    //  [flattened] F = B_j * q_j + N_mat * n_;
    //  dF/dq_j = B_j + N * ((1/|n_|)(I - n * n^T)([Dq1][-I O I] - [Dq2][-I I O]))
    Eigen::Matrix99d Bj;
    matrix_B_triangle(Bj, dX);
    std::vector<Eigen::Vector3d> Xj(3);
    std::vector<Eigen::Vector3d> qj(3);
    for (int i = 0; i < 3; i ++) {
        Xj[i] = V.row(element(i));
        qj[i] = q.segment(element(i)*3, 3);
    }
    Eigen::Vector3d N_ = triangle_normal(Xj[0], Xj[1], Xj[2]).normalized();
    Eigen::Matrix93d N; N.setZero();
    for (int i = 0; i < 3; i ++) {
        for (int k = 0; k < 3; k ++) {
            N(i + 3 * k, k) = N_(i);
        }
    }
    Eigen::Vector3d n_ = triangle_normal(qj[0], qj[1], qj[2]);
    double nnorm = n_.norm();
    Eigen::Vector3d n = n_.normalized();
    Eigen::Matrix3d Dq1, Dq2;
    Eigen::Vector3d dq1 = qj[1] - qj[0], dq2 = qj[2] - qj[0];
    skew_matrix(Dq1, dq1); 
    skew_matrix(Dq2, dq2); 
    Eigen::Matrix39d sel1, sel2;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d O = Eigen::Matrix3d::Zero();
    sel1 << -I, I, O; 
    sel2 << -I, O, I;
    Eigen::Matrix39d dn_dq = (1.0/nnorm) * (I - n * n.transpose()) * (Dq1 * sel2 - Dq2 * sel1);
    dFdq = Bj + N * dn_dq;
}