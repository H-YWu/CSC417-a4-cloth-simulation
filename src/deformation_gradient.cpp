#include <deformation_gradient.h>

Eigen::Vector3d triangle_normal(Eigen::Vector3d X0,
                     Eigen::Vector3d X1,
                     Eigen::Vector3d X2) {
    Eigen::Vector3d V01 = X1 - X0;
    Eigen::Vector3d V02 = X2 - X0;
    return V01.cross(V02);
}

void deformation_gradient(Eigen::Matrix3d &F,
                        Eigen::Ref<const Eigen::VectorXd> q, 
                        Eigen::Ref<const Eigen::MatrixXd> V,
                        Eigen::Ref<const Eigen::RowVectorXi> element,
                        Eigen::Ref<const Eigen::Matrix3d> dphidX){
    // the 3x3 delta X matrix of a single triangle
    std::vector<Eigen::Vector3d> x, X;
    Eigen::Matrix34d left;
    for (int i = 0; i < 3; i ++) {
        X[i] = V.row(element(i));
        x[i] = q.segment(element(i)*3, 3);
        left.col(i) = x[i]; 
    }
    left.col(3) = triangle_normal(x[0], x[1], x[2]);
    Eigen::Matrix43d right;
    right.block(0, 0, 3, 3) = dphidX;
    right.row(3) = triangle_normal(X[0], X[1], X[2]).transpose();;
    // F = \sum_i x_i \diff{phi_i}{X} + n * N^T
    F = left * right;
}