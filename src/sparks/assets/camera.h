#pragma once
#include "glm/glm.hpp"

namespace sparks {
class Camera {
 public:
  Camera(float fov = 60.0f, float aperture = 0.0f, float focal_length = 3.0f) : 
  fov_(fov), aperture_(aperture), focal_length_(focal_length) {}
  [[nodiscard]] glm::mat4 GetProjectionMatrix(float aspect) const;
  void GenerateRay(float aspect,
                   glm::vec2 range_low,
                   glm::vec2 range_high,
                   glm::vec3 &origin,
                   glm::vec3 &direction,
                   float rand_u = 0.0f,
                   float rand_v = 0.0f,
                   float rand_w = 0.0f,
                   float rand_r = 0.0f) const;
  bool ImGuiItems();
  void UpdateFov(float delta);

 private:
  float fov_{60.0f};
  float aperture_{0.0f};
  float focal_length_{3.0f};
};
}  // namespace sparks
