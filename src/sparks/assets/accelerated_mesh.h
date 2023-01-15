#pragma once
#include "sparks/assets/aabb.h"
#include "sparks/assets/mesh.h"

namespace sparks {

struct KdAccelNode {
  // KdAccelNode Methods
  int s[2]{-1,-1};
  int offset=-1;
  int start=-1;
  AxisAlignedBoundingBox box;
  int split=-1;
  float split_pos=-1;
};
struct TreeSetting {
  int isectCost = 80;
  int traversalCost = 1;
  float emptyBonus = 0.5;
  int maxPrims = 1;
  int maxDepth = -1;
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
  void build_tree(TreeSetting &setting,
                  std::vector<int> &face_ids,
                  std::vector<AxisAlignedBoundingBox> all_aabb,
                  KdAccelNode &node,
                  int depth,
                  int failed_time);
 void Init(KdAccelNode &node, std::vector<int> &face_id);
 private:
  std::vector<KdAccelNode> tree_node;
  std::vector<int> ordered;

};
}  // namespace sparks
