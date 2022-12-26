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
  float result = -1.0f;
  treeIntersect(treeRoot, origin, direction, t_min, result, hit_record);
  return result;
  //return Mesh::TraceRay(origin, direction, t_min, hit_record);
}

void AcceleratedMesh::BuildAccelerationStructure() {
  t.resize(indices_.size() / 3);
  for (int i = 0; i < indices_.size() / 3; i++) {
    for (int j = 0; j < 3; j++)
      t[i].center += vertices_[indices_[3 * i + j]].position;
    t[i].center /= 3.0f;
    t[i].aabb = AxisAlignedBoundingBox(vertices_[indices_[3 * i]].position);
    for (int j = 1; j < 3; j++) {
      t[i].aabb |= AxisAlignedBoundingBox(vertices_[indices_[3 * i + j]].position);
    }
  }
  std::vector<int> faceIds;
  for (int i = 0; i < t.size(); i++)
    faceIds.push_back(i);
  buildTree(treeRoot, faceIds);
  std::cout << "# of triangles: " << t.size() << std::endl;
}

bool AcceleratedMesh::triangleIntersect(int x, const glm::vec3 &origin, const glm::vec3 &direction, float t_min, float &t_max, HitRecord *hit_record) const {
  int i = 3 * x, j = 3 * x + 1, k = 3 * x + 2;
  const auto &v0 = vertices_[indices_[i]];
  const auto &v1 = vertices_[indices_[j]];
  const auto &v2 = vertices_[indices_[k]];

  glm::mat3 A = glm::mat3(v1.position - v0.position, v2.position - v0.position,
                          -direction);
  if (std::abs(glm::determinant(A)) < 1e-9f) {
    return false;
  }
  A = glm::inverse(A);
  auto uvt = A * (origin - v0.position);
  auto &t = uvt.z;
  if (t < t_min || (t_max > 0.0f && t > t_max)) {
    return false;
  }
  auto &u = uvt.x;
  auto &v = uvt.y;
  auto w = 1.0f - u - v;
  auto position = origin + t * direction;
  if (u >= 0.0f && v >= 0.0f && u + v <= 1.0f) {
    t_max = t;
    if (hit_record) {
      auto geometry_normal = glm::normalize(
          glm::cross(v2.position - v0.position, v1.position - v0.position));
      if (glm::dot(geometry_normal, direction) < 0.0f) {
        hit_record->position = position;
        hit_record->geometry_normal = geometry_normal;
        hit_record->normal = v0.normal * w + v1.normal * u + v2.normal * v;
        hit_record->tangent = v0.tangent * w + v1.tangent * u + v2.tangent * v;
        hit_record->tex_coord =
            v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
        hit_record->front_face = true;
      } else {
        hit_record->position = position;
        hit_record->geometry_normal = -geometry_normal;
        hit_record->normal = -(v0.normal * w + v1.normal * u + v2.normal * v);
        hit_record->tangent =
            -(v0.tangent * w + v1.tangent * u + v2.tangent * v);
        hit_record->tex_coord =
            v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
        hit_record->front_face = false;
      }
    }
    return true;
  } else {
    return false;
  }
}

bool AcceleratedMesh::treeIntersect(int x, const glm::vec3 &origin, const glm::vec3 &direction, float t_min, float &t_max, HitRecord *hit_record) const {
  if (x == -1)
    return false;
  float tt_max = t_max > t_min ? t_max : 1e3f;
  if (!t[x].aabb.IsIntersect(origin, direction, t_min, tt_max))
    return false;

  const TreeNode &triIndex = t[x];
  bool ret = triangleIntersect(x, origin, direction, t_min, t_max, hit_record);

  bool chret = false;
  int id = relativeId(t[x].center, origin);
  for (int i = 0; i < 8; i++) {
    chret |= treeIntersect(t[x].child[id ^ Ids[i]], origin, direction,
                                t_min, t_max, hit_record);
  }
  if (chret == true)
    return chret;
  return ret;
}

int AcceleratedMesh::relativeId(glm::vec3 a, glm::vec3 b) const {
  int ret = 0;
  if (a[0] < b[0])
    ret |= 1;
  if (a[1] < b[1])
    ret |= 2;
  if (a[2] < b[2])
    ret |= 4;
  return ret;
}

void AcceleratedMesh::buildTree(int &x, std::vector<int> faceIds) {
  x = -1;
  if (faceIds.size() == 0)
    return;
  std::vector<int> divide[8];
  std::vector<std::pair<int, int>> score;
  score.resize(faceIds.size());

  for (int i = 0; i < faceIds.size(); i++)
    score[i] = std::make_pair(0, faceIds[i]);
  int center = faceIds.size() / 2;
  sort(score.begin(), score.end(),
       [this](std::pair<int, int> u, std::pair<int, int> v) -> bool {
         return this->t[u.second].center[0] < this->t[v.second].center[0];
       });
  for (int i = 0; i < faceIds.size(); i++)
    score[i].first += abs(i - center);
  sort(score.begin(), score.end(),
       [this](std::pair<int, int> u, std::pair<int, int> v) -> bool {
         return this->t[u.second].center[1] < this->t[v.second].center[1];
       });
  for (int i = 0; i < faceIds.size(); i++)
    score[i].first += abs(i - center);
  sort(score.begin(), score.end(),
       [this](std::pair<int, int> u, std::pair<int, int> v) -> bool {
         return this->t[u.second].center[2] < this->t[v.second].center[2];
       });
  for (int i = 0; i < faceIds.size(); i++)
    score[i].first += abs(i - center);
  sort(score.begin(), score.end());
  x = score[0].second;

  for (int i = 1; i < faceIds.size(); i++) {
    int id = score[i].second;
    divide[relativeId(t[x].center, t[id].center)].push_back(id);
  }
  for (int i = 0; i < 8; i++) {
    buildTree(t[x].child[i], divide[i]);
    if (t[x].child[i] != -1)
      t[x].aabb |= t[t[x].child[i]].aabb;
  }
}

}  // namespace sparks


