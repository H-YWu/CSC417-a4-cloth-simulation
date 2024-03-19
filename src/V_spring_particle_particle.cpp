#include <V_spring_particle_particle.h>

void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    // q = [q0; q1]
    // dx = Bq
    //  B = [-I, I]
    Eigen::Vector3d dx = q1 - q0;
    // V = 1/2 k (\sqrt{dx^T dx} - l0)^2
    V = 0.5 * stiffness * pow(sqrt(dx.dot(dx)) - l0, 2);
}