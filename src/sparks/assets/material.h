#pragma once
#include "cstdint"
#include "glm/glm.hpp"
#include "sparks/assets/util.h"
#include "sparks/assets/hit_record.h"
#include "sparks/assets/bxdf.h"
#include "sparks/assets/bssrdf.h"
#include "sparks/assets/medium.h"
#include "sparks/assets/subsurface_lib.h"

namespace sparks {

enum MaterialType : int {
  MATERIAL_TYPE_LAMBERTIAN = 0,
  MATERIAL_TYPE_SPECULAR = 1,
  MATERIAL_TYPE_TRANSMISSIVE = 2,
  MATERIAL_TYPE_PRINCIPLED = 3,
  MATERIAL_TYPE_EMISSION = 4,
  MATERIAL_TYPE_SUBSURFACE = 5,
  MATERIAL_TYPE_KDSUBSURFACE = 6,
  MATERIAL_TYPE_MEDIUM = 7,
  MATERIAL_TYPE_GRID_MEDIUM = 8
};

class Scene;

// New features appended: Jingyi Lyu + Shengquan Du.
class Material {
public:
  glm::vec3 albedo_color{0.8f};
  int albedo_texture_id{0};
  bool use_normal_texture{false};
  int normal_texture_id{0};
  glm::vec3 emission{0.0f};
  float emission_strength{1.0f};
  float cone{0.0f};
  bool use_alpha_texture{false};
  int alpha_texture_id{0};
  float alpha{1.0f};
  MaterialType material_type{MATERIAL_TYPE_LAMBERTIAN};
  float metallic{0.0f}, eta{1.45f};
  float roughness{0.5f};
  float specularTint{0.0f};
  float anisotropic{0.0f};
  float sheen{0.0f}, sheenTint{0.5f};
  float clearcoat{0.0f}, clearcoatGloss{0.0f};
  float specTrans{0.0f};
  float flatness{0.0f}, diffTrans{0.0f}; //used by thin surface
  bool thin{false};
  glm::vec3 scatterDistance{0.0f};
  bool use_ss_texture{false};
  int sigma_a_texture_id{0};
  int sigma_s_texture_id{0};
  PresetSSType ss_type{SS_NONE};
  glm::vec3 volumetric_emission{0.0f};
  Medium *medium{nullptr};
  bool false_surface{false};
  glm::vec3 sigma_a{0.0f}, sigma_s{0.0f};
  float g{0.0f};
  glm::vec3 mfp{0.0f};
  BSSRDFTable *table{nullptr};
  float reserve[2]{};
  Material() = default;
  Material(const Material &) = default;
  Material& operator = (const Material &) = default;
  explicit Material(const glm::vec3 &albedo);
  Material(Scene *scene, const tinyxml2::XMLElement *material_element);
  glm::vec3 GetAlbedoColor(const HitRecord &hit, const Scene *scene) const;
  float GetAlpha(const HitRecord &hit, const Scene *scene) const;
  glm::vec3 GetShaderNormal(const HitRecord &hit, const Scene *scene) const;
  HitRecord GetShaderHit(const HitRecord &hit, const Scene *scene) const;
  BSDF* ComputeBSDF(const HitRecord &hit, const Scene* scene) const;
  BSSRDF* ComputeBSSRDF(const HitRecord &hit, const glm::vec3 direction, const Scene* scene);
};
}  // namespace sparks
