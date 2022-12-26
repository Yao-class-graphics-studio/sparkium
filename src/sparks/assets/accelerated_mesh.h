#pragma once
#include "sparks/assets/aabb.h"
#include "sparks/assets/mesh.h"

namespace sparks {

namespace {
struct TreeNode {
  AxisAlignedBoundingBox aabb{};
  int child[8]{};
  glm::vec3 center;
};
}  // namespace

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
  /*
   * You can add your acceleration structure contents here.
   * */
  int treeRoot;
  int intersectCnt;
  std::vector<TreeNode> t;
  const int Ids[8]{0, 1, 2, 4, 6, 5, 3, 7};
  void buildTree(int &x, std::vector<int> faceIds);
  int relativeId(glm::vec3 a, glm::vec3 b) const;
  bool treeIntersect(int x,
                     const glm::vec3 &origin,
                     const glm::vec3 &direction,
                     float t_min,
                     float &t_max,
                     HitRecord *hit_record) const;
  bool triangleIntersect(int x,
                          const glm::vec3 &origin,
                          const glm::vec3 &direction,
                          float t_min,
                          float &t_max,
                          HitRecord *hit_record) const;
};
}  // namespace sparks
