#include "sparks/assets/aabb.h"

#include "algorithm"
#include "grassland/grassland.h"

namespace sparks {

AxisAlignedBoundingBox::AxisAlignedBoundingBox(float x_low,
                                               float x_high,
                                               float y_low,
                                               float y_high,
                                               float z_low,
                                               float z_high) {
  this->x_low = x_low;
  this->x_high = x_high;
  this->y_low = y_low;
  this->y_high = y_high;
  this->z_low = z_low;
  this->z_high = z_high;
}

AxisAlignedBoundingBox::AxisAlignedBoundingBox(const glm::vec3 &position) {
  x_low = position.x;
  x_high = position.x;
  y_low = position.y;
  y_high = position.y;
  z_low = position.z;
  z_high = position.z;
}

AxisAlignedBoundingBox::AxisAlignedBoundingBox(const glm::vec3& pos1,
                                               const glm::vec3& pos2,
                                               const glm::vec3& pos3) {
  x_low = std::min(std::min(pos1.x, pos2.x), pos3.x);
  x_high = std::max(std::max(pos1.x, pos2.x), pos3.x);
  y_low = std::min(std::min(pos1.y, pos2.y), pos3.y);
  y_high = std::max(std::max(pos1.y, pos2.y), pos3.y);
  z_low = std::min(std::min(pos1.z, pos2.z), pos3.z);
  z_high = std::max(std::max(pos1.z, pos2.z), pos3.z);
}

AxisAlignedBoundingBox::AxisAlignedBoundingBox(
    const AxisAlignedBoundingBox &aabb) {
  x_low = aabb.x_low;
	x_high = aabb.x_high;
	y_low = aabb.y_low;
	y_high = aabb.y_high;
	z_low = aabb.z_low;
	z_high = aabb.z_high;
}

bool AxisAlignedBoundingBox::IsIntersect(const glm::vec3 &origin,
                                         const glm::vec3 &direction,
                                         float t_min,
                                         float t_max) const {
  if (x_low <= origin.x && origin.x <= x_high && y_low <= origin.y &&
      origin.y <= y_high && z_low <= origin.z && origin.z <= z_high) {
    return true;
  }
  float intersection_range_low = t_max * (1.0f + t_min);
  float intersection_range_high = 0.0f;
  float t;
  glm::vec3 intersection;
#define TestIntersection(x, y, z)                                     \
  if (std::abs(direction.x) > 1e-5) {                                 \
    float inv_d = 1.0f / direction.x;                                 \
    t = (x##_low - origin.x) * inv_d;                                 \
    intersection = origin + direction * t;                            \
    if (y##_low <= intersection.y && intersection.y <= y##_high &&    \
        z##_low <= intersection.z && intersection.z <= z##_high) {    \
      intersection_range_low = std::min(intersection_range_low, t);   \
      intersection_range_high = std::max(intersection_range_high, t); \
    }                                                                 \
    t = (x##_high - origin.x) * inv_d;                                \
    intersection = origin + direction * t;                            \
    if (y##_low <= intersection.y && intersection.y <= y##_high &&    \
        z##_low <= intersection.z && intersection.z <= z##_high) {    \
      intersection_range_low = std::min(intersection_range_low, t);   \
      intersection_range_high = std::max(intersection_range_high, t); \
    }                                                                 \
  }
  TestIntersection(x, y, z);
  TestIntersection(z, x, y);
  TestIntersection(y, z, x);
  return intersection_range_high >= t_min && intersection_range_low <= t_max;
}

bool AxisAlignedBoundingBox::OverlapWith(
    const AxisAlignedBoundingBox &aabb) const {
  return x_low <= aabb.x_high && aabb.x_low <= x_high && y_low <= aabb.y_high &&
         aabb.y_low <= y_high && z_low <= aabb.z_high && aabb.z_low <= z_high;
}

std::optional<std::pair<float, float>> AxisAlignedBoundingBox::Intersect(
    const glm::vec3 &origin,
    const glm::vec3 &direction,
    float t_min,
    float t_max) const {
  float intersection_range_low = t_max * (1.0f + t_min);
  float intersection_range_high = 0.0f;
  float t;
  glm::vec3 intersection;
#define TestIntersection(x, y, z)                                     \
  if (std::abs(direction.x) > 1e-5) {                                 \
    float inv_d = 1.0f / direction.x;                                 \
    t = (x##_low - origin.x) * inv_d;                                 \
    intersection = origin + direction * t;                            \
    if (y##_low <= intersection.y && intersection.y <= y##_high &&    \
        z##_low <= intersection.z && intersection.z <= z##_high) {    \
      intersection_range_low = std::min(intersection_range_low, t);   \
      intersection_range_high = std::max(intersection_range_high, t); \
    }                                                                 \
    t = (x##_high - origin.x) * inv_d;                                \
    intersection = origin + direction * t;                            \
    if (y##_low <= intersection.y && intersection.y <= y##_high &&    \
        z##_low <= intersection.z && intersection.z <= z##_high) {    \
      intersection_range_low = std::min(intersection_range_low, t);   \
      intersection_range_high = std::max(intersection_range_high, t); \
    }                                                                 \
  }
  TestIntersection(x, y, z);
  TestIntersection(z, x, y);
  TestIntersection(y, z, x);
  if (intersection_range_high >= t_min && intersection_range_low <= t_max) {
    return std::make_pair(std::max(intersection_range_low, t_min),
                          std::min(intersection_range_high, t_max));
  } else {
    return {};
  }
}

AxisAlignedBoundingBox AxisAlignedBoundingBox::operator&(
    const AxisAlignedBoundingBox &aabb) const {
  return {std::max(x_low, aabb.x_low), std::min(x_high, aabb.x_high),
          std::max(y_low, aabb.y_low), std::min(y_high, aabb.y_high),
          std::max(z_low, aabb.z_low), std::min(z_high, aabb.z_high)};
}

AxisAlignedBoundingBox AxisAlignedBoundingBox::operator|(
    const AxisAlignedBoundingBox &aabb) const {
  return {std::min(x_low, aabb.x_low), std::max(x_high, aabb.x_high),
          std::min(y_low, aabb.y_low), std::max(y_high, aabb.y_high),
          std::min(z_low, aabb.z_low), std::max(z_high, aabb.z_high)};
}

AxisAlignedBoundingBox &AxisAlignedBoundingBox::operator&=(
    const AxisAlignedBoundingBox &aabb) {
  x_low = std::max(x_low, aabb.x_low);
  x_high = std::min(x_high, aabb.x_high);
  y_low = std::max(y_low, aabb.y_low);
  y_high = std::min(y_high, aabb.y_high);
  z_low = std::max(z_low, aabb.z_low);
  z_high = std::min(z_high, aabb.z_high);
  return *this;
}

AxisAlignedBoundingBox &AxisAlignedBoundingBox::operator|=(
    const AxisAlignedBoundingBox &aabb) {
  x_low = std::min(x_low, aabb.x_low);
  x_high = std::max(x_high, aabb.x_high);
  y_low = std::min(y_low, aabb.y_low);
  y_high = std::max(y_high, aabb.y_high);
  z_low = std::min(z_low, aabb.z_low);
  z_high = std::max(z_high, aabb.z_high);
  return *this;
}

}  // namespace sparks
