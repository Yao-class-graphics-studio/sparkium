#pragma once
#include "cstdint"
#include "glm/glm.hpp"
#include "sparks/assets/util.h"
#include "sparks/assets/hit_record.h"
#include "sparks/assets/bxdf.h"

namespace sparks {

enum MaterialType : int {
  MATERIAL_TYPE_LAMBERTIAN = 0,
  MATERIAL_TYPE_SPECULAR = 1,
  MATERIAL_TYPE_TRANSMISSIVE = 2,
  MATERIAL_TYPE_PRINCIPLED = 3,
  MATERIAL_TYPE_EMISSION = 4
};

class Scene;

struct Material {
  glm::vec3 albedo_color{0.8f};
  int albedo_texture_id{0};
  glm::vec3 emission{0.0f};
  float emission_strength{1.0f};
  float alpha{1.0f};
  MaterialType material_type{MATERIAL_TYPE_LAMBERTIAN};
  float metallic{0.0f}, eta{1.0f};
  float roughness{0.0f};
  float specularTint{0.0f};
  float anisotropic{0.0f};
  float sheen{0.0f}, sheenTint{0.0f};
  float clearcoat{0.0f}, clearcoatGloss{0.0f};
  float specTrans{0.0f};
  float flatness{0.0f}, diffTrans{0.0f}; //used by thin surface
  bool thin{false};
  glm::vec3 scatterDistance{0.0f};
  float reserve[2]{};
  Material() = default;
  explicit Material(const glm::vec3 &albedo);
  Material(Scene *scene, const tinyxml2::XMLElement *material_element);
  BSDF* ComputeBSDF(const HitRecord &hit,
					const Scene* scene) const;
};
}  // namespace sparks
