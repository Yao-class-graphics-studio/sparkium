#include "sparks/assets/accelerated_mesh.h"

#include "algorithm"

namespace sparks {
AcceleratedMesh::AcceleratedMesh(const Mesh &mesh) : Mesh(mesh) {
  BuildAccelerationStructure();
}

AcceleratedMesh::AcceleratedMesh(const std::vector<Vertex> &vertices,
                                 const std::vector<uint32_t> &indices)
    : Mesh(vertices, indices) {
  BuildAccelerationStructure();
}

float AcceleratedMesh::TraceRay(const glm::vec3 &origin,
                                const glm::vec3 &direction,
                                float t_min,
                                HitRecord *hit_record) const {
  return Mesh::TraceRay(origin, direction, t_min, hit_record);
}

void AcceleratedMesh::BuildAccelerationStructure() {
}

}  // namespace sparks
