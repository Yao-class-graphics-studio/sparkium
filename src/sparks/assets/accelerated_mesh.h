#pragma once
#include "sparks/assets/aabb.h"
#include "sparks/assets/mesh.h"
#include <vector>

namespace sparks {

struct TreeNode {
  AxisAlignedBoundingBox aabb{};
  int child[2]{-1, -1};
  std::vector<uint32_t> indices{};
  bool IsLeaf() const { return child[0] == -1 && child[1] == -1; }
};

class AcceleratedMesh : public Mesh {
 public:
  AcceleratedMesh() = default;
  explicit AcceleratedMesh(const Mesh &mesh);
  AcceleratedMesh(const std::vector<Vertex> &vertices,
                  const std::vector<uint32_t> &indices);
  float TraceRay(const glm::vec3 &origin,
                 const glm::vec3 &direction,
                 float t_min,
                 HitRecord *hit_record) const override;
  void BuildAccelerationStructure();

 private:
  std::vector<TreeNode> nodes;
  int root;
  int BuildSubtree(const AxisAlignedBoundingBox &aabb,
                                         const std::vector<uint32_t> &indices);
  float TraceRay(const glm::vec3 &origin,
                 const glm::vec3 &direction,
                 float t_min,
                 HitRecord *hit_record,
                 int node_idx) const;
};
}  // namespace sparks
