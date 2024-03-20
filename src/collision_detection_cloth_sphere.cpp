#include <collision_detection_cloth_sphere.h>
#include <iostream>
void collision_detection_cloth_sphere(std::vector<unsigned int> &cloth_index, std::vector<Eigen::Vector3d> &normals, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Vector3d> center, double radius) {

    cloth_index.clear();
    normals.clear();

    double r2 = radius * radius;
    int n = q.size() / 3;

    for (int i = 0; i < n; i ++) {
        Eigen::Vector3d p = q.segment(i*3, 3);
        Eigen::Vector3d r2p = p - center;
        if (r2p.dot(r2p) <= r2) {
            cloth_index.push_back(i);
            normals.push_back(r2p.normalized());
        }
    }
}