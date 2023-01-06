#pragma once
#include "cstdint"
#include "glm/glm.hpp"

namespace sparks {

enum MaterialType : uint32_t {
  MATERIAL_TYPE_LAMBERTIAN = 0,
  MATERIAL_TYPE_SPECULAR = 1,
  MATERIAL_TYPE_TRANSMISSIVE = 2,
  MATERIAL_TYPE_PRINCIPLED = 3,
  MATERIAL_TYPE_EMISSION = 4
};

struct Material {
  glm::vec3 ambient{0.0f};
  int ambient_texture_id{0};
  glm::vec3 diffuse{0.8f};
  int diffuse_texture_id{0};
  glm::vec3 specular{0.0f};
  int specular_texture_id{0};
  glm::vec3 transmittance{0.0f};
  float roughness{0.0f};
  glm::vec3 emission{0.0f};
  float ior{1.0f};

  MaterialType material_type{MATERIAL_TYPE_LAMBERTIAN};
  int reserve[3];
};
}  // namespace sparks
