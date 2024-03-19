#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    // q = [q0; q1]
    // dx = Bq (q0 -> q1)
    //  B = [-I, I]
    // V = 1/2 k (\sqrt{dx^T dx} - l0)^2
    Eigen::Vector3d dx = q1 - q0;
    // transpose
    Eigen::Vector6d BTdx;
    BTdx << -dx, dx;
    double dxdot = dx.dot(dx);
    // here f is not force, force = - grad
    // f := grad = \partial{V}{q} = \partial{V}{dx} \diff{dx}{q}
    //                 = \partial{V}{dx} B
    //                 = 1/2 k (\sqrt{dx^T dx} - l0) \frac{dx^T}{\sqrt{dx^T dx}} B
    // f = [f0; f1]
    // however here f is a column vector
    //  so we use the B^Tdx = (dx^TB)^T
    f = 0.5 * stiffness * (sqrt(dxdot) - l0) * pow(dxdot, -0.5) * BTdx;
}