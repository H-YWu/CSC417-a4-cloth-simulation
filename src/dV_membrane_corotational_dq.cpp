#include <dV_membrane_corotational_dq.h>
#include <deformation_gradient.h>
#include <dF_cloth_triangle_dq.h>

void dV_membrane_corotational_dq(Eigen::Vector9d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {

    //Deformation Gradient
    Eigen::Matrix3d dx; //deformed tangent matrix 
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 

    deformation_gradient(dx, q, V, element, dX);
    Eigen::JacobiSVD<Eigen::Matrix3d> SVD(dx);
    S = SVD.singularValues();
    U = SVD.matrixU();
    W = SVD.matrixV();

    //Fix for inverted elements (thanks to Danny Kaufman)
    double det = S[0]*S[1];
    
     if(det <= -1e-10)
    {
        if(S[0] < 0) S[0] *= -1;
        if(S[1] < 0) S[1] *= -1;
        if(S[2] < 0) S[2] *= -1;
    }
    
    if(U.determinant() <= 0)
    {
        U(0, 2) *= -1;
        U(1, 2) *= -1;
        U(2, 2) *= -1;
    }
    
    if(W.determinant() <= 0)
    {
        W(0, 2) *= -1;
        W(1, 2) *= -1;
        W(2, 2) *= -1;
    }
    
    // dV/ds
    Eigen::Matrix3d dVds;
    double tmp = lambda * (S(0) + S(1) + S(2) - 3.0);
    for (int i = 0; i < 3; i ++) {
        dVds.coeffRef(i, i) = 2.0 * mu * (S(i) - 1.0) + tmp;
    }
    // dV/dF
    Eigen::Matrix3d dV_mat = U * dVds * W.transpose();
    for (int i = 0; i < 3; i ++) {
        for (int k = 0; k < 3; k ++) {
            dV(i * 3 + k) = dV_mat(i, k);
        }
    }
    // dV/dq_j
    Eigen::Matrix99d dFdq;
    dF_cloth_triangle_dq(dFdq, q, dX, V, element);
    //  dV/dq_j = (dF/dq_j)^T * dV/dF
    dV = area * dFdq.transpose() * dV;
}