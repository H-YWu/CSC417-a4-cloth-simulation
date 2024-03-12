#include <dphi_cloth_triangle_dX.h>

//compute 3x3 deformation gradient 
void dphi_cloth_triangle_dX(Eigen::Matrix3d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    // lambda = (T^T * T)^{-1} * T^T * (X - X0) 
    // \partial{lambda}{X} = \partial{(T^T * T)^{-1} * T^T * (X - X0)}{X}
    //                     = \partial{(T^T * T)^{-1} * T^T * (X - X0)}{(X - X0)} \diff{(X - X0)}{X}
    //                     = (T^T * T)^{-1} * T^T * I
    //                     = (T^T * T)^{-1} * T^T
    // phi = [1 - lambda[+01]; lambda]
    // \partial{phi}{X} = [-\partial{lambda}{X}[+01]; \partial{lambda}{X}]
    std::vector<Eigen::Vector3d> XX(3);
    for (int i = 0; i < 3; i ++) {
        XX[i] = V.row(element(i));
    }
    Eigen::MatrixXd T(3, 2);
    for (int i = 0; i < 2; i ++) {
        T.col(i) = XX[i+1] - XX[0];
    }
    Eigen::MatrixXd TT = T.transpose();
    Eigen::MatrixXd TTT = T.transpose() * T;
    Eigen::MatrixXd TTTinv = TTT.inverse();
    Eigen::MatrixXd TTTinvTT = TTT.inverse() * TT;
    dphi.row(0) << 0., 0., 0.;
    for (int i = 0; i < 2; i ++) {
        dphi.row(i+1) = TTTinvTT.row(i);
        dphi.row(0) -= TTTinvTT.row(i);
    }
}