#pragma once
#include "sparks/assets/aabb.h"
#include "sparks/assets/mesh.h"

namespace sparks {

namespace {
struct TreeNode {
  AxisAlignedBoundingBox aabb{};
  int child[2]{-1, -1};
};
}  // namespace

// Acceleration structure: Shengquan Du.

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
  float TraceRayBVH(int root,
                    int i_begin,
                    int i_end,
                    const glm::vec3 &origin,
                    const glm::vec3 &direction,
                    float t_min,
                    HitRecord *hit_record) const;
  void BuildBVH(int root, int i_begin, int i_end, int partition_dimension);
  std::vector<std::pair<AxisAlignedBoundingBox, int>> origin_bounding_boxes_;
  std::vector<AxisAlignedBoundingBox> treenode_bounding_boxes_;
  constexpr static int LEAFSIZE = 16;
  /*
   * You can add your acceleration structure contents here.
   * */
};
}  // namespace sparks
