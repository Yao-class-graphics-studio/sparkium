#pragma once
#include "glm/glm.hpp"

namespace sparks {
struct AxisAlignedBoundingBox {
  float x_low{};
  float x_high{};
  float y_low{};
  float y_high{};
  float z_low{};
  float z_high{};
  AxisAlignedBoundingBox(float x_low,
                         float x_high,
                         float y_low,
                         float y_high,
                         float z_low,
                         float z_high);
  AxisAlignedBoundingBox(const glm::vec3 &position = glm::vec3{0.0f});
  [[nodiscard]] bool IsIntersect(const glm::vec3 &origin,
                                 const glm::vec3 &direction,
                                 float t_min,
                                 float t_max,
                                 float &intersection_range_low,
                                 float &intersection_range_high) const;
  AxisAlignedBoundingBox operator&(const AxisAlignedBoundingBox &aabb) const;
  AxisAlignedBoundingBox operator|(const AxisAlignedBoundingBox &aabb) const;
  AxisAlignedBoundingBox &operator&=(const AxisAlignedBoundingBox &aabb);
  AxisAlignedBoundingBox &operator|=(const AxisAlignedBoundingBox &aabb);
  glm::vec3 Diagnal() {
    return glm::vec3(x_high - x_low, y_high - y_low, z_high - z_low);
  }
  float SurfaceArea() {
    glm::vec3 d = Diagnal();
    return 2 * (d.x * d.y + d.y * d.z + d.z * d.x);
  }
  glm::vec3 pMin() {
    return glm::vec3(x_low, y_low, z_low);
  }
  glm::vec3 pMax() {
    return glm::vec3(x_high, y_high, z_high);
  }
};
}  // namespace sparks
