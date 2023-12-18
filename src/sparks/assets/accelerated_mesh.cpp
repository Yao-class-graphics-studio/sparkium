#include "sparks/assets/accelerated_mesh.h"

#include "algorithm"

const int MAX_TRIANGLE = 5;
const float MIN_VOLUME = 1e-9;
const float eps = 1e-3;

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
  return TraceRay(origin, direction, t_min, hit_record, root);
}

void AcceleratedMesh::BuildAccelerationStructure(){ root = BuildSubtree(GetAABB(glm::mat4{1.0f}), indices_); }

int AcceleratedMesh::BuildSubtree(
    const AxisAlignedBoundingBox &aabb,
    const std::vector<uint32_t> &indices) {
  TreeNode node{};
  node.aabb = aabb;
  if (indices.size() <= MAX_TRIANGLE * 3 || aabb.Volume() <= MIN_VOLUME) {
    node.indices = indices;
    nodes.push_back(node);
		return nodes.size() - 1;
	}

  // find the longest axis and split
  int axis = -1;
  float length = 0;
  if (aabb.x_high - aabb.x_low > length) {
    axis = 0;
    length = aabb.x_high - aabb.x_low;
  }
  if (aabb.y_high - aabb.y_low > length) {
		axis = 1;
		length = aabb.y_high - aabb.y_low;
	}
  if (aabb.z_high - aabb.z_low > length) {
    axis = 2;
    length = aabb.z_high - aabb.z_low;
	}
  AxisAlignedBoundingBox aabb_1 = aabb;
  AxisAlignedBoundingBox aabb_2 = aabb;
  if (axis == 0) {
    float mid = (aabb.x_high + aabb.x_low) / 2;
    aabb_1.x_high = mid;
    aabb_2.x_low = mid;
  } else if (axis == 1) {
		float mid = (aabb.y_high + aabb.y_low) / 2;
		aabb_1.y_high = mid;
		aabb_2.y_low = mid;
	} else if (axis == 2) {
		float mid = (aabb.z_high + aabb.z_low) / 2;
		aabb_1.z_high = mid;
		aabb_2.z_low = mid;
	}

  // distribute faces
  std::vector<uint32_t> indices_1;
  std::vector<uint32_t> indices_2;
  for (int i = 0; i < indices.size(); i += 3) {
    int j = i + 1, k = i + 2;
    const auto &v0 = vertices_[indices[i]];
    const auto &v1 = vertices_[indices[j]];
    const auto &v2 = vertices_[indices[k]];
    AxisAlignedBoundingBox aabb_triangle{v0.position, v1.position, v2.position};
    if (aabb_triangle.OverlapWith(aabb_1)) {
			indices_1.push_back(indices[i]);
			indices_1.push_back(indices[j]);
			indices_1.push_back(indices[k]);
		}
    if (aabb_triangle.OverlapWith(aabb_2)) {
			indices_2.push_back(indices[i]);
			indices_2.push_back(indices[j]);
			indices_2.push_back(indices[k]);
		}
  }

  // create children
  node.child[0] = BuildSubtree(aabb_1, indices_1);
  node.child[1] = BuildSubtree(aabb_2, indices_2);

  nodes.push_back(node);
  return nodes.size() - 1;
}

float AcceleratedMesh::TraceRay(const glm::vec3 &origin,
                                const glm::vec3 &direction,
                                float t_min,
                                HitRecord *hit_record,
                                int node_idx) const {
  auto node = &nodes[node_idx];
  auto ret = node->aabb.Intersect(origin, direction, t_min, 1e6f);
  if (!ret) return -1.0f;

  if (node->IsLeaf()) {
    t_min = ret->first;
    float t_max = ret->second;
    float result = -1.0f;
    auto &indices = node->indices;
    for (int i = 0; i < indices.size(); i += 3) {
      int j = i + 1, k = i + 2;
      const auto &v0 = vertices_[indices[i]];
      const auto &v1 = vertices_[indices[j]];
      const auto &v2 = vertices_[indices[k]];

      glm::mat3 A = glm::mat3(v1.position - v0.position,
                              v2.position - v0.position, -direction);
      if (std::abs(glm::determinant(A)) < 1e-9f) {
        continue;
      }
      A = glm::inverse(A);
      auto uvt = A * (origin - v0.position);
      auto &t = uvt.z;
      if (t < t_min - eps || t > t_max + eps || (result > 0.0f && t > result)) {
        continue;
      }
      auto &u = uvt.x;
      auto &v = uvt.y;
      auto w = 1.0f - u - v;
      auto position = origin + t * direction;
      if (u >= 0.0f && v >= 0.0f && u + v <= 1.0f) {
        result = t;
        if (hit_record) {
          auto geometry_normal = glm::normalize(
              glm::cross(v2.position - v0.position, v1.position - v0.position));
          if (glm::dot(geometry_normal, direction) < 0.0f) {
            hit_record->position = position;
            hit_record->geometry_normal = geometry_normal;
            hit_record->normal = v0.normal * w + v1.normal * u + v2.normal * v;
            hit_record->tangent =
                v0.tangent * w + v1.tangent * u + v2.tangent * v;
            hit_record->tex_coord =
                v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
            hit_record->front_face = true;
          } else {
            hit_record->position = position;
            hit_record->geometry_normal = -geometry_normal;
            hit_record->normal =
                -(v0.normal * w + v1.normal * u + v2.normal * v);
            hit_record->tangent =
                -(v0.tangent * w + v1.tangent * u + v2.tangent * v);
            hit_record->tex_coord =
                v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
            hit_record->front_face = false;
          }
        }
      }
		}
    return result;
  }

  // internal node
  auto ret_1 = nodes[node->child[0]].aabb.Intersect(origin, direction, t_min, 1e6f);
  auto ret_2 = nodes[node->child[1]].aabb.Intersect(origin, direction, t_min, 1e6f);
  if (!ret_1) {
    return TraceRay(origin, direction, t_min, hit_record, node->child[1]);
  } else if (!ret_2) {
    return TraceRay(origin, direction, t_min, hit_record, node->child[0]);
  } else {
    int first = ret_1->first < ret_2->first ? 0 : 1;
    int second = 1 - first;
    float result = TraceRay(origin, direction, t_min, hit_record, node->child[first]);
    if (result < 0)
      return TraceRay(origin, direction, t_min, hit_record, node->child[second]);
    return result;
  }
}

}  // namespace sparks
