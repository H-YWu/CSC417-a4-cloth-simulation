#include <d2V_membrane_corotational_dq2.h>
#include <dF_cloth_triangle_dq.h>
#include <deformation_gradient.h>
#include <dsvd.h>

void d2V_membrane_corotational_dq2(Eigen::Matrix99d &H, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    

    //SVD = USW^T
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 
    Eigen::Matrix3d F; //deformation gradient
    
    double tol = 1e-5;
    
    //Compute SVD of F here
    deformation_gradient(F, q, V, element, dX);
    Eigen::JacobiSVD<Eigen::Matrix3d> SVD(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    S = SVD.singularValues();
    U = SVD.matrixU();
    W = SVD.matrixV();
    
    //deal with singularity in the svd gradient
    if(std::fabs(S[0] - S[1]) < tol || std::fabs(S[1] - S[2]) < tol || std::fabs(S[0] - S[2]) < tol) {
        F += Eigen::Matrix3d::Random()*tol;
        Eigen::JacobiSVD<Eigen::Matrix3d> svd2(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd2.matrixU();
        W = svd2.matrixV();
        S = svd2.singularValues();
    }
    
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

    // Compute H, the hessian of the corotational energy
    Eigen::Matrix3d VT = W.transpose();
    Eigen::Matrix99d dFdq;
    // compute dFdq
    dF_cloth_triangle_dq(dFdq, q, dX, V, element);
    // compute SVD partial derivatives
    Eigen::Tensor3333d dU, dV;
    Eigen::Tensor333d dS;
    dsvd(dU, dS, dV, F);
    // d(dphi/dF)/dF_ij = dU/dF_ij S V^T + U dS/dF_ij V^T + U S dV^T/dF_ij
    Eigen::Tensor3333d d2phi;
    Eigen::Matrix3d delS, del2S;
    delS.setZero(); del2S.setZero();
    double tmp = lambda * (S(0) + S(1) + S(2) - 3.0);
    for (int i = 0; i < 3; i ++) {
        delS.coeffRef(i, i) = 2.0 * mu * (S(i) - 1.0) + tmp;
        for (int j = 0; j < 3; j ++) {
            del2S.coeffRef(i, j) = lambda * S(j);
            if (i == j) del2S.coeffRef(i, j) += 2.0 * mu;
        }
    }
    for (int i = 0; i < 3; i ++) {
        for (int j = 0; j < 3; j ++) {
            Eigen::Vector3d dsij = del2S * dS[i][j];
            Eigen::DiagonalMatrix<double, 3> dsij_mat(dsij);
            d2phi[i][j] = dU[i][j] * delS * VT
                        + U * dsij_mat * VT
                        + U * delS * dV[i][j].transpose();
        }
    }
    H.setZero();
    // Flatten the 3x3x3x3 tensor to a 9x9 matrix
    for (int i1 = 0; i1 < 3; i1 ++) {
        for (int j1 = 0; j1 < 3; j1 ++) {
            for (int i2 = 0; i2 < 3; i2 ++) {
                for (int j2 = 0; j2 < 3; j2 ++) {
                    //H(3 * i1 + j1, 3 * i2 + j2) = d2phi[i1][j1](i2, j2);
                    H(3 * i1 + j1, 3 * i2 + j2) = d2phi[i2][j2](i1, j1);
                }
            }
        }
    }
    // H_j = (dF/dq_j)^T * (d2phi/F2) * (dF/dq_j)
    H = area * dFdq.transpose() * H * dFdq;

    //fix errant eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix99d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 9; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }
    
    H = Evec * DiagEval * Evec.transpose();
    
}
