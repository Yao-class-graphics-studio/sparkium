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
struct KdToDo {
  const KdAccelNode *node;
  float tMin, tMax;
};
float AcceleratedMesh::TraceRay(const glm::vec3 &origin,
                                const glm::vec3 &direction,
                                float t_min,
                                HitRecord *hit_record) const {
  //return Mesh::TraceRay(origin, direction, t_min, hit_record);
  float sectMin, sectMax;
  float RayMax = 1e5f;
  bool hit = true;
  if (!tree_node[0].box.IsIntersect(origin, direction, t_min, RayMax, sectMin,
                              sectMax))
         return false;
    glm::vec3 invDir(1 / direction.x, 1 / direction.y, 1 / direction.z);
  constexpr int maxTodo = 64;
  KdToDo todo[maxTodo];
  int todoPos = 0;
  const KdAccelNode *node = &tree_node[0];
  while (node != nullptr) {
    if (RayMax < sectMin)
      break;
    if (node->start == -1) {
      int axis = node->split;
      float tPlane = (node->split_pos - origin[axis]) * invDir[axis];
      const KdAccelNode *firstChild, *secondChild;
      int belowFirst = (origin[axis] < node->split_pos) ||
                       (origin[axis] == node->split_pos && direction[axis] <= 0);
      if (belowFirst) {
        firstChild = &tree_node[node->s[0]];
        secondChild = &tree_node[node->s[1]];
      } else {
        firstChild = &tree_node[node->s[1]];
        secondChild = &tree_node[node->s[0]];
      }
      if (tPlane > sectMax || tPlane <= 0)
        node = firstChild;
      else if (tPlane < sectMin)
        node = secondChild;
      else {
        todo[todoPos].node = secondChild;
        todo[todoPos].tMin = tPlane;
        todo[todoPos].tMax = sectMax;
        ++todoPos;

        node = firstChild;
        sectMax = tPlane;
      }
    } else {
      int n_face = node->offset;
      for (int i = 0; i < n_face; ++i) {
        int index = ordered[node->start + i];
        const auto &v0 = vertices_[indices_[index*3]];
        const auto &v1 = vertices_[indices_[index * 3+1]];
        const auto &v2 = vertices_[indices_[index * 3+2]];
        glm::mat3 A = glm::mat3(v1.position - v0.position,
                                v2.position - v0.position, -direction);
        if (std::abs(glm::determinant(A)) < 1e-9f) {
          continue;
        }
        A = glm::inverse(A);
        auto uvt = A * (origin - v0.position);
        auto &t = uvt.z;
        if (t < t_min || (RayMax > 0.0f && t > RayMax)) {
          continue;
        }
        auto &u = uvt.x;
        auto &v = uvt.y;
        auto w = 1.0f - u - v;
        auto position = origin + t * direction;
        if (u >= 0.0f && v >= 0.0f && u + v <= 1.0f) {
          RayMax = t;
          hit = true;
          if (hit_record) {
            auto geometry_normal = glm::normalize(glm::cross(
                v2.position - v0.position, v1.position - v0.position));
            if (glm::dot(geometry_normal, direction) < 0.0f) {
              hit_record->position = position;
              hit_record->geometry_normal = geometry_normal;
              hit_record->normal =
                  v0.normal * w + v1.normal * u + v2.normal * v;
              hit_record->tangent =
                  v0.tangent * w + v1.tangent * u + v2.tangent * v;
              hit_record->tex_coord =
                  v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
              hit_record->front_face = true;
              hit_record->index_id = i / 3;
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
              hit_record->index_id = i;
            }
          }
        }
      }
      if (todoPos > 0) {
        --todoPos;
        node = todo[todoPos].node;
        sectMin = todo[todoPos].tMin;
        sectMax = todo[todoPos].tMax;
      }
      else break;
    }
  }
  if (!hit)
    RayMax = -1;
  return RayMax;
}

void AcceleratedMesh::Init(KdAccelNode &node, std::vector<int> &face_id) {
  node.start = ordered.size();
  node.offset = face_id.size();
  //printf("has");
  for (auto u : face_id) {
    //printf("%d", u);
    ordered.push_back(u);
  }
  //printf("\n");
}
void AcceleratedMesh::build_tree(TreeSetting &setting, std::vector<int> &face_ids,std::vector<AxisAlignedBoundingBox> all_aabb,int node_id,int depth,int failed_time) {
  int n_face = face_ids.size();
  if (n_face <= setting.maxPrims || depth == setting.maxDepth) {
    Init(tree_node[node_id], face_ids);
    return;
  }
  float oldCost = setting.isectCost * n_face;
  int bestAxis = -1, bestOffset = -1;
  float bestCost = 1e30f;
  float totalSA = tree_node[node_id].box.SurfaceArea();
  float invTotalSA = 1 / totalSA;
  glm::vec3 d = tree_node[node_id].box.Diagnal();
  for (int axis = 0; axis < 3; ++axis) {
    for (int i = 0; i < n_face; ++i) {
      int pid = face_ids[i];
      AxisAlignedBoundingBox &bounds = all_aabb[pid];
      edges[axis][2 * i] = BoundEdge(bounds.pMin()[axis], pid, true);
      edges[axis][2 * i + 1] = BoundEdge(bounds.pMax()[axis], pid, false);
    }
    std::sort(&edges[axis][0], &edges[axis][2 * n_face],
              [](const BoundEdge &e0, const BoundEdge &e1) -> bool {
                if (e0.t == e1.t)
                  return (int)e0.type < (int)e1.type;
                else
                  return e0.t < e1.t;
              });
    int nBelow = 0, nAbove = n_face;
    for (int i = 0; i < 2 * n_face; ++i) {
      if (edges[axis][i].type == EdgeType::End)
        --nAbove;
      float edgeT = edges[axis][i].t;
      if (edgeT > tree_node[node_id].box.pMin()[axis] &&
          edgeT < tree_node[node_id].box.pMax()[axis]) {
        int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
        float belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                             (edgeT - tree_node[node_id].box.pMin()[axis]) *
                                 (d[otherAxis0] + d[otherAxis1]));
        float aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                             (tree_node[node_id].box.pMax()[axis] - edgeT) *
                                 (d[otherAxis0] + d[otherAxis1]));
        float pBelow = belowSA * invTotalSA;
        float pAbove = aboveSA * invTotalSA;
        float eb = (nAbove == 0 || nBelow == 0) ? setting.emptyBonus : 0;
        float cost = setting.traversalCost +
                     setting.isectCost * (1 - eb) * (pBelow * nBelow + pAbove * nAbove);
        if (cost < bestCost) {
          bestCost = cost;
          bestAxis = axis;
          bestOffset = i;
        }
      }
      if (edges[axis][i].type == EdgeType::Start)
        ++nBelow;
    }
  }
  if (bestCost > oldCost)
    failed_time++;
  if (bestAxis == -1 || (bestCost > oldCost && n_face<16)||failed_time>=3) {
    Init(tree_node[node_id], face_ids);
    return;
  }
  std::vector<int> Lid, Rid;
  for (int i = 0; i < bestOffset; ++i)
    if (edges[bestAxis][i].type == EdgeType::Start)
      Lid.push_back(edges[bestAxis][i].primNum);
  for (int i = bestOffset + 1; i < 2 * n_face; ++i)
    if (edges[bestAxis][i].type == EdgeType::End)
      Rid.push_back(edges[bestAxis][i].primNum);
  float tSplit = edges[bestAxis][bestOffset].t;
  AxisAlignedBoundingBox bounds0 = tree_node[node_id].box,
                         bounds1 = tree_node[node_id].box;
  if (bestAxis == 0) {
    bounds0.x_high = bounds1.x_low = tSplit;
  }else if (bestAxis == 1) {
    bounds0.y_high = bounds1.y_low = tSplit;
  } else if (bestAxis == 2) {
    bounds0.z_high = bounds1.z_low = tSplit;
  }
  KdAccelNode lson, rson;
  lson.box = bounds0;
  rson.box = bounds1;
  tree_node[node_id].split = bestAxis;
  tree_node[node_id].split_pos = tSplit;
  tree_node[node_id].start = -1;
  tree_node[node_id].s[0] = tree_node.size();
  tree_node.push_back(lson);
  //printf("to %d\n", tree_node[node_id].s[0]);
  build_tree(setting, Lid, all_aabb, tree_node[node_id].s[0], depth + 1,
             failed_time);
  tree_node[node_id].s[1] = tree_node.size();
  tree_node.push_back(rson);
  //printf("to %d\n", tree_node[node_id].s[1]);
  build_tree(setting, Rid, all_aabb, tree_node[node_id].s[1], depth + 1,
             failed_time);
}
void AcceleratedMesh::BuildAccelerationStructure() {
  TreeSetting setting;
  std::vector<int> face_ids;
  std::vector<AxisAlignedBoundingBox> all_aabb;
  AxisAlignedBoundingBox test(vertices_[indices_[0]].position);
  for (int i = 0,j=0; i < indices_.size(); i += 3,j+=1) {
    face_ids.push_back(j);
    AxisAlignedBoundingBox tp =
        AxisAlignedBoundingBox(vertices_[indices_[i + 0]].position) |
        AxisAlignedBoundingBox(vertices_[indices_[i + 1]].position) |
        AxisAlignedBoundingBox(vertices_[indices_[i + 2]].position);
    all_aabb.push_back(tp);
    test |= tp;
  }
  setting.maxDepth = std::round(8 + 1.3f * log2(face_ids.size()));
  tree_node.push_back(KdAccelNode{});
  tree_node[0].box = test;
  int tot = face_ids.size();
  printf("tot %d\n", tot);
  //printf("to 0\n");
  for (int i = 0; i < 3; ++i)
    edges[i] = new BoundEdge[2*tot];
  build_tree(setting, face_ids, all_aabb, 0,0,0);
  for (int i = 0; i < 3; ++i)
    delete edges[i];
  return;
}

}  // namespace sparks
