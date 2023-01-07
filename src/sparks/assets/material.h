#pragma once
#include "cstdint"
#include "glm/glm.hpp"
#include "sparks/assets/util.h"

namespace sparks {

enum MaterialType : uint32_t {
  MATERIAL_TYPE_LAMBERTIAN = 0,
  MATERIAL_TYPE_SPECULAR = 1,
  MATERIAL_TYPE_TRANSMISSIVE = 2,
  MATERIAL_TYPE_PRINCIPLED = 3,
  MATERIAL_TYPE_EMISSION = 4
};

class Scene;

struct Material {
  // glm::vec3 ambient{0.0f};
  // int ambient_texture_id{0};
  glm::vec3 diffuse{1.0f};
  int diffuse_texture_id{0};

  glm::vec3 specular{0.0f};
  int specular_texture_id{0};

  glm::vec3 opacity{0.0f};
  int opacity_texture_id{1};

  glm::vec3 emission{0.0f};
  int emission_texture_id{0};

  glm::vec3 transmittance{0.0f};
  float ior{1.0f};

  float roughness{0.0f};
  int roughness_texture_id{0};
  float metallic{0.0f};
  int metallic_texture_id{0};

  float sheen{0.0f};
  int sheen_texture_id{0};
  float clearcoat_thickness{0.0f};
  float clearcoat_roughness{0.0f};

  float anisotropy{0.0f};
  float anisotropy_rotation{0.0f};
  int normal_texture_id{0};
  MaterialType material_type{MATERIAL_TYPE_LAMBERTIAN};
  // float reserve[2]{};
  Material() = default;
  explicit Material(const glm::vec3 &albedo);
  Material(Scene *scene, const tinyxml2::XMLElement *material_element);
};
}  // namespace sparks
