#pragma once
#include <optional>
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
  AxisAlignedBoundingBox(const glm::vec3 &pos1,
                         const glm::vec3 &pos2,
                         const glm::vec3 &pos3);
  AxisAlignedBoundingBox(const AxisAlignedBoundingBox &aabb);
  [[nodiscard]] bool IsIntersect(const glm::vec3 &origin,
                                 const glm::vec3 &direction,
                                 float t_min,
                                 float t_max) const;
  bool OverlapWith(const AxisAlignedBoundingBox &aabb) const;
  [[nodiscard]] std::optional<std::pair<float, float>> Intersect(const glm::vec3 &origin,
                                                                 const glm::vec3 &direction,
                                                                 float t_min,
                                                                 float t_max) const;
  float Volume() const { return (x_high - x_low) * (y_high - y_low) * (z_high - z_low); }
  AxisAlignedBoundingBox operator&(const AxisAlignedBoundingBox &aabb) const;
  AxisAlignedBoundingBox operator|(const AxisAlignedBoundingBox &aabb) const;
  AxisAlignedBoundingBox &operator&=(const AxisAlignedBoundingBox &aabb);
  AxisAlignedBoundingBox &operator|=(const AxisAlignedBoundingBox &aabb);
};
}  // namespace sparks
