#include <velocity_filter_cloth_sphere.h>

void velocity_filter_cloth_sphere(Eigen::VectorXd &qdot, const std::vector<unsigned int> &indices, 
                                  const std::vector<Eigen::Vector3d> &normals) {
    for (int i = 0; i < indices.size(); i ++) {
        int id = indices[i];
        Eigen::Vector3d n = normals[i];
        Eigen::Vector3d vel = qdot.segment(i*3, 3);
        vel = vel - n.dot(vel) * n;
        qdot.segment<3>(i) = vel;
    }
}