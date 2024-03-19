#include <matrix_B_triangle.h>

void matrix_B_triangle(Eigen::Matrix99d &B, Eigen::Ref<const Eigen::Matrix3d> dX) {
    B.setZero(9, 9);
    // auxilary matrix D = dX
    // B(j + 3k, 3i + k) = D(i, j) for i = 0..2, j = 0..2, k = 0..2
    std::vector<Eigen::Triplet<double>> triples;
    for (int i = 0; i < 3; i ++) {
        for (int j = 0; j < 3; j ++) {
            for (int k = 0; k < 3; k ++) {
                B(j+k*3, i*3+k) = dX(i, j);
            }
        }
    }
}