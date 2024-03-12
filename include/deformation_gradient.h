#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
// q - generalized coordinates of FEM system
// V - vertex matrix for the mesh
// element - vertex indices of the element
// dphidX - the 3x3 matrix containing dphi/dX
//Output:
// F - the 3x3 deformation gradient of a single triangle 
void deformation_gradient(Eigen::Matrix3d &F,
                        Eigen::Ref<const Eigen::VectorXd> q, 
                        Eigen::Ref<const Eigen::MatrixXd> V,
                        Eigen::Ref<const Eigen::RowVectorXi> element,
                        Eigen::Ref<const Eigen::Matrix3d> dphidX);