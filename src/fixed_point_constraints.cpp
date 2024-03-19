#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
    // q = P^T \hat{q} + q_{fixed}
    //  P(k * 3 + j, i * 3 + j) = 1 for k = all indices and i = not fixed indices in q; j = 0, 1, 2
    int n = q_size / 3;
    P.resize(3 * (n - indices.size()), q_size);
    int id = 0, k = 0;
    // Check all indices
    for (unsigned int i = 0; i < n; i ++) {
        if (i == indices[id]) { // jump fixed
            id ++;
        } else {
            for (int j = 0; j < 3; j ++) {
                P.insert(k + j, i * 3 + j) = 1;
            }
            k += 3;
        }
    }
}