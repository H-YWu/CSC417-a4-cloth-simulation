#include <V_membrane_corotational.h>
#include <deformation_gradient.h>

//Allowed to use libigl SVD or Eigen SVD for this part
void V_membrane_corotational(double &energy, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    Eigen::Matrix3d F;
    deformation_gradient(F, q, V, element, dX);
    Eigen::JacobiSVD<Eigen::Matrix3d> SVD(F);
    Eigen::Vector3d sigma = SVD.singularValues();
    energy = 0.0;
    double tmp = -3.0;
    for (int i = 0; i < 3; i ++) {
        energy += pow(sigma(i) - 1.0, 2);
        tmp += sigma(i);
    }
    energy *= mu;
    energy += 0.5 * lambda * pow(tmp, 2);
}
