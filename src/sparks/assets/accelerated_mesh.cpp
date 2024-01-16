#include "sparks/assets/accelerated_mesh.h"

#include "algorithm"

// Acceleration structure: Shengquan Du
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
  return TraceRayBVH(1,0,origin_bounding_boxes_.size(),origin, direction, t_min, hit_record);
    //return Mesh::TraceRay(origin, direction, t_min, hit_record);
}

float AcceleratedMesh::TraceRayBVH(int root, 
                                int i_begin, 
                                int i_end, 
                                const glm::vec3 &origin,
                                const glm::vec3 &direction,
                                float t_min,
                                HitRecord *hit_record) const {
  if (!treenode_bounding_boxes_[root].IsIntersect(origin, direction, t_min,
                                                  1e4f))
    return -1;
  if (i_end - i_begin <= LEAFSIZE)
  {
    float result = -1;
    for (int i = i_begin; i < i_end; i++)
    {
      const auto &v0 = vertices_[indices_[origin_bounding_boxes_[i].second * 3]];
      const auto &v1 =
          vertices_[indices_[origin_bounding_boxes_[i].second * 3 +1]];
      const auto &v2 =
          vertices_[indices_[origin_bounding_boxes_[i].second * 3 + 2]];
      TraceRayTriangleUpdate(origin, direction, t_min, hit_record, v0, v1, v2,
                             result);
    }
    return result;
  }
  HitRecord local_hit_record_left;
  HitRecord local_hit_record_right;
  int mid = (i_begin + i_end + 1) >> 1;
  float result_left =
      TraceRayBVH(root << 1, i_begin, mid, origin, direction, t_min,
                  hit_record ? &local_hit_record_left : nullptr);
  float result_right =
      TraceRayBVH(root << 1 | 1, mid, i_end, origin, direction, t_min,
                  hit_record ? &local_hit_record_right : nullptr);
  if (result_right <= t_min || (result_left > t_min && result_left < result_right))
  {
    if (hit_record)
        *hit_record = local_hit_record_left;
    return result_left;
  }
  else if (result_right > t_min)
  {
    if (hit_record)
        *hit_record = local_hit_record_right;
    return result_right;
  }
  return -1;
}

void AcceleratedMesh::BuildBVH(int root, int i_begin,
                               int i_end,
                               int partition_dimension) {
  if (i_end - i_begin <= LEAFSIZE) {
    treenode_bounding_boxes_[root] = origin_bounding_boxes_[i_begin].first;
    for (int i = i_begin + 1; i < i_end; i++)
      treenode_bounding_boxes_[root] |= origin_bounding_boxes_[i].first;
      return;
  }
  auto get_box_high = [partition_dimension](const AxisAlignedBoundingBox &a) -> float {
    switch (partition_dimension) {
      case 0:
        return a.x_high;
      case 1:
        return a.y_high;
      case 2:
        return a.z_high;
    }
    return a.x_high;
  };
  using T = std::pair<AxisAlignedBoundingBox, int>;
  sort(origin_bounding_boxes_.begin() + i_begin, origin_bounding_boxes_.begin() + i_end,
       [&get_box_high](const T &a, const T &b) {
         return std::make_pair(get_box_high(a.first), a.second) < std::make_pair(get_box_high(b.first), b.second);
       });
  int mid = (i_begin + i_end + 1) >> 1;
  int next_partition_dimension = (partition_dimension + 1) % 3;
  BuildBVH(root << 1, i_begin, mid, next_partition_dimension);
  BuildBVH(root << 1 | 1, mid, i_end, next_partition_dimension);
  treenode_bounding_boxes_[root] = treenode_bounding_boxes_[root << 1] |
                                   treenode_bounding_boxes_[root << 1 | 1];
}

void AcceleratedMesh::BuildAccelerationStructure() {
  origin_bounding_boxes_.resize(indices_.size() / 3);
  int leaf_lowbound = (LEAFSIZE >> 1) + 1;
  //int leaf_lowbound = 1;
  treenode_bounding_boxes_.resize((origin_bounding_boxes_.size() + leaf_lowbound -1)/leaf_lowbound * 4);
  for (int i = 0; i * 3 < indices_.size(); i++){
    auto triangle_box = (
        AxisAlignedBoundingBox(vertices_[indices_[i * 3]].position) |
        AxisAlignedBoundingBox(vertices_[indices_[i * 3 + 1]].position) |
        AxisAlignedBoundingBox(vertices_[indices_[i * 3 + 2]].position)
    );
    origin_bounding_boxes_[i] = std::make_pair(triangle_box, i); 
  }
  
  BuildBVH(1, 0, origin_bounding_boxes_.size(), 0);
  
}

}  // namespace sparks
