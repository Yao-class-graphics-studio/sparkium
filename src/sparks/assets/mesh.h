#pragma once
#include "sparks/assets/model.h"
#include "sparks/assets/util.h"
#include "sparks/assets/vertex.h"
#include "vector"

namespace sparks {
class Mesh : public Model {
 public:
  Mesh() = default;
  Mesh(const Mesh &mesh);
  Mesh(const std::vector<Vertex> &vertices,
       const std::vector<uint32_t> &indices);
  explicit Mesh(const tinyxml2::XMLElement *element);
  ~Mesh() override = default;
  [[nodiscard]] float TraceRay(const glm::vec3 &origin,
                               const glm::vec3 &direction,
                               float t_min,
                               HitRecord *hit_record) const override;
  const char *GetDefaultEntityName() override;
  [[nodiscard]] AxisAlignedBoundingBox GetAABB(
      const glm::mat4 &transform) const override;
  [[nodiscard]] std::vector<Vertex> GetVertices() const override;
  [[nodiscard]] std::vector<uint32_t> GetIndices() const override;
  static Mesh Cube(const glm::vec3 &center, const glm::vec3 &size);
  static Mesh Sphere(const glm::vec3 &center = glm::vec3{0.0f},
                     float radius = 1.0f);
  static bool LoadObjFile(const std::string &obj_file_path, Mesh &mesh);
  void WriteObjFile(const std::string &file_path) const;
  void MergeVertices();
  double area(int idx) const{
    Vertex a = vertices_[indices_[3 * idx]], b = vertices_[indices_[3 * idx + 1]],
           c = vertices_[indices_[3 * idx + 2]];
    auto pa = a.position, pb = b.position, pc = c.position;
    auto v_ab = pb - pa;
    auto v_ac = pc - pa;
    auto cross = glm::cross(v_ab, v_ac);
    return glm::length(cross) / 2;
  };
  Vertex vertice(int idx) {
    return vertices_[indices_[idx]];
  }
  const std::vector<Vertex> &Get_vertices() const{
    return vertices_;
  };
  const std::vector<uint32_t> &Get_indices() const{
    return indices_;
  };

 protected:
  std::vector<Vertex> vertices_;
  std::vector<uint32_t> indices_;
};
}  // namespace sparks
