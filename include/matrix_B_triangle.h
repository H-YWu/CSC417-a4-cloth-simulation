#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  dX - the 3x3 matrix containing dphi/dX
//Output:
// B - the 9x9 B auxilary matrix of a single triangle to calculate stacked F
void matrix_B_triangle(Eigen::Matrix99d &B, Eigen::Ref<const Eigen::Matrix3d> dX);