#include <mass_matrix_mesh.h>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, 
                         Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                         double density, Eigen::Ref<const Eigen::VectorXd> areas) {
    M.resize(q.size(), q.size());
    Eigen::Matrix3d Me;
    // generated using MATLAB symbolic toolkit
    Me(0, 0) = 1.0/2.0;
    Me(0, 1) = 1.0/4.0;
    Me(0, 2) = 1.0/4.0;
    Me(1, 0) = 1.0/4.0;
    Me(1, 1) = 1.0/2.0;
    Me(1, 2) = 1.0/4.0;
    Me(2, 0) = 1.0/4.0;
    Me(2, 1) = 1.0/4.0;
    Me(2, 2) = 1.0/2.0;
    Me *= density;
    // Assume h = 1.0
    Eigen::Matrix99d m;
    for (int i = 0; i < 3; i ++) {
        for (int j = 0; j < 3; j ++) {
            for (int k = 0; k < 3; k ++) {
                m(i * 3 + k, j * 3 + k) = Me(i, j); 
            }
        }
    }  

    // all triangles
    std::vector<Eigen::Triplet<double>> triples;
    for (int r = 0; r < F.rows(); r ++) {
        // E(3i + j, 3 * idx[i] + j)   = 1 for i = 0..2, j = 0..2
        // E^T(3 * idx[i] + j, 3i + j) = 1 for i = 0..2, j = 0..2
        // m(x, y) for x,y = 0..11
        //  (m E)(x, 3 * idx[i] + j) += m(x, 3i + j) for x = 0..11, i = 0..2, j = 0..2
        // M: (E^T m E) (3 * idx[I] + J, 3 * idx[i] + j) += m(3I + J, 3i + j) for I = 0..2, J = 0..2, i = 0..2, j = 0..2
        // area of this triangle
        double a = areas(r);
        for (int I = 0; I < 3; I ++) {
            for (int J = 0; J < 3; J ++) {
                for (int i = 0; i < 3; i ++) {
                    for (int j = 0; j < 3; j ++) {
                        triples.push_back(Eigen::Triplet<double>(3*element(I)+J, 3*element(i)+j, m(3*I+J,3*i+j)*a));
                    }
                }
            }
        }
    }
    M.setFromTriplets(triples.begin(), triples.end());
}
 