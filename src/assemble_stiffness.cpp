#include <assemble_stiffness.h>
#include <d2V_membrane_corotational_dq2.h>
#include <mass_matrix_mesh.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0, 
                     double mu, double lambda) { 
    K.resize(qdot.size(), qdot.size());
    // all triangles 
    std::vector<Eigen::Triplet<double>> triples;
    for (int r = 0; r < F.rows(); r ++) {
        Eigen::RowVector3i element = F.row(r);
        Eigen::Matrix99d H;
        H.setZero();
        Eigen::Matrix<double,1,9> tmp_row;
        tmp_row = dX.row(r);
        d2V_membrane_corotational_dq2(H, q, Eigen::Map<const Eigen::Matrix3d>(tmp_row.data()), V, element, a0(r), mu, lambda);
        // K = - \partial^2{V}{q} = - \sum \partial^2{V_r}{q}
        //   = - \sum (E^T {H} E)_r
        //  E(3i + j, 3 * idx[i] + j)   = 1 for i = 0..2, j = 0..2
        //  E^T(3 * idx[i] + j, 3i + j) = 1 for i = 0..2, j = 0..2
        //  H(x, y) for x,y = 0..8
        //  (H E)(x, 3 * idx[i] + j) -= H(x, 3i + j) for x = 0..11, i = 0..2, j = 0..2
        //  (E^T H E)(3 * idx[I] + J, 3 * idx[i] + j) -= H(3I + J, 3i + j) for I = 0..2, J = 0..2, i = 0..2, j = 0..2
        for (int I = 0; I < 3; I ++) {
            for (int J = 0; J < 3; J ++) {
                for (int i = 0; i < 3; i ++) {
                    for (int j = 0; j < 3; j ++) {
                        triples.push_back(Eigen::Triplet<double>(3*element(I)+J, 3*element(i)+j, -H(3*I+J,3*i+j)));
                    }
                }
            }
        }
    }
    K.setFromTriplets(triples.begin(), triples.end());    
};
