#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//Output:
//  dFdq - the per-triangle partial derivative of deformation gradient w.r.t. generalized coordinates q
void dF_cloth_triangle_dq(Eigen::Matrix99d &dFdq, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX,
                        Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element);