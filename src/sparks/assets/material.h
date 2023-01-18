#pragma once
#include "cstdint"
#include "glm/glm.hpp"
#include "sparks/assets/util.h"
#include "sparks/assets/hit_record.h"
#include "random"
#include "glm/gtc/constants.hpp"
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
  float etaA{1.00029f};
  float etaB{1.5f};
  Material() = default;
  explicit Material(const glm::vec3 &albedo);
  Material(Scene *scene, const tinyxml2::XMLElement *material_element);
  glm::vec3 WorldToLocal(HitRecord &hit_record, glm::vec3 v) const{
    glm::vec3 n = hit_record.normal;
    int index = abs(n.x) < abs(n.y) ? (abs(n.x) < abs(n.z) ? 0 : 2)
                                    : (abs(n.y) < abs(n.z) ? 1 : 2);
    glm::vec3 s(0.f);
    s[index] = 1;
    s = normalize(cross(n, s));
    glm::vec3 t = normalize(cross(n, s));
    return glm::vec3(dot(v, s), dot(v, t), dot(v, n));
  }
  glm::vec3 LocalToWorld(HitRecord &hit_record, glm::vec3 v) const{
    glm::vec3 n = hit_record.normal;
    int index = abs(n.x) < abs(n.y)
        ? (abs(n.x) < abs(n.z) ? 0 : 2)
           : (abs(n.y) < abs(n.z) ? 1 : 2);
    glm::vec3 s(0.f);
    s[index] = 1;
    s = normalize(cross(n, s));
    glm::vec3 t = normalize(cross(n, s));
    return glm::vec3(s.x * v.x + t.x * v.y + n.x * v.z,
                     s.y * v.x + t.y * v.y + n.y * v.z,
                     s.z * v.x + t.z * v.y + n.z * v.z);
  }
  bool IsLambertian() const{
    return material_type == MATERIAL_TYPE_LAMBERTIAN;
  }
  bool IsSpecular() const {
    return material_type == MATERIAL_TYPE_SPECULAR||material_type==MATERIAL_TYPE_TRANSMISSIVE;
  }
  bool IsMirror() const {
    return material_type == MATERIAL_TYPE_SPECULAR;
  }
  bool IsTransmissive() const{
    return material_type == MATERIAL_TYPE_TRANSMISSIVE;
  }
  bool IsEmission() const{
    return material_type == MATERIAL_TYPE_EMISSION;
  }
  glm::vec3 local_f_Lambertian(glm::vec3 wo, glm::vec3 wi) const{
    return (wo.z * wi.z > 0) ?glm::one_over_pi<float>() * albedo_color:glm::vec3(0.f);
  }
  glm::vec3 local_f_Transmissive(glm::vec3 wo, glm::vec3 wi) const{
    return glm::vec3(0.f);
  }
  glm::vec3 local_f_Specular(glm::vec3 wo, glm::vec3 wi) const{
    return glm::vec3(0.f);
  }
  float local_pdf(glm::vec3 wo, glm::vec3 wi) const{
    if (IsLambertian()) {
      //return 1;
      return (wo.z * wi.z > 0) ? abs(wi.z) * glm::one_over_pi<float>() : 0;
    }
    if (IsSpecular())
      return 0.0f;
    return 0.0f;
  }
  glm::vec3 Sample_HemiSphere(std::mt19937 &rd) const{
    glm::vec2 u(std::uniform_real_distribution<float>(0.0f, 1.0f)(rd),
                std::uniform_real_distribution<float>(0.0f, 1.0f)(rd));
    glm::vec2 offset=2.f * u - glm::vec2(1, 1);
    glm::vec2 disk(0.0f);
    float r=0, theta=0;
    if (offset.x == 0 && offset.y == 0) {
      disk = glm::vec2(0.0f);
    } else{
      if (abs(offset.x) > abs(offset.y)) {
        r = offset.x;
        theta = glm::quarter_pi<float>() * (offset.y / offset.x);
      } else {
        r = offset.y;
        theta = glm::half_pi<float>() -
                glm::quarter_pi<float>() * (offset.x / offset.y);
      }
      disk=r * glm::vec2(std::cos(theta), std::sin(theta));
    }
    float z = sqrt(std::max(0.f, 1 - disk.x * disk.x - disk.y * disk.y));
    return glm::vec3(disk.x, disk.y, z);
  }
  float FrDielectric(float cosThetaI, float etaI, float etaT) const{
    cosThetaI = std::max(-1.0f,std::min(1.0f,cosThetaI));
    bool entering =
        cosThetaI > 0.f;
    if (!entering) {
      std::swap(etaI, etaT);
      cosThetaI = std::abs(cosThetaI);
    }

    float sinThetaI =
        std::sqrt(std::max((float)0, 1 - cosThetaI * cosThetaI));
    float sinThetaT = etaI / etaT * sinThetaI;
    float cosThetaT =
        std::sqrt(std::max((float)0, 1 - sinThetaT * sinThetaT));

    float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                  ((etaT * cosThetaI) + (etaI * cosThetaT));
    float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
                  ((etaI * cosThetaI) + (etaT * cosThetaT));
    return (Rparl * Rparl + Rperp * Rperp) / 2;
  }
  inline bool Refract(glm::vec3 &wi,
                      glm::vec3 &n,
                      float eta,
                      glm::vec3 *wt) const{
    float cosThetaI = dot(n, wi);
    float sin2ThetaI = std::max(0.f, 1.f - cosThetaI * cosThetaI);
    float sin2ThetaT = eta * eta * sin2ThetaI;
    if (sin2ThetaT >= 1) return false;
    float cosThetaT = std::sqrt(1 - sin2ThetaT);
    *wt = eta * -wi + (eta * cosThetaI - cosThetaT) * n;
    return true;
  }
  glm::vec3 Sample_f(HitRecord &hit_record,
                     std::mt19937 &rd,
                     glm::vec3 &woW,
                     glm::vec3 *wiW,
                     float *pdf) const{//note: we are sampling reflection functions
    glm::vec3 wo = WorldToLocal(hit_record, woW);
    //return -s;
    glm::vec3 f(0.f);

    if (IsLambertian()) {
      glm::vec3 wi = Sample_HemiSphere(rd);
      if (wo.z < 0)
        wi.z *= -1;
      *pdf = local_pdf(wo, wi);
      *wiW = LocalToWorld(hit_record, wi);
      f += local_f_Lambertian(wo, wi);
    } else if (IsMirror()) {
      *wiW = LocalToWorld(hit_record, glm::vec3(-wo.x, -wo.y, wo.z));
      *pdf = abs(dot(*wiW, hit_record.normal));  
      f += albedo_color;
    } else if (IsTransmissive()) {
      float F = FrDielectric(wo.z, etaA, etaB);
      if (std::uniform_real_distribution<float>(0.0f, 1.0f)(rd) < F) {
        glm::vec3 wi(-wo.x, -wo.y, wo.z);
        *wiW = LocalToWorld(hit_record, glm::vec3(-wo.x, -wo.y, wo.z));
        *pdf = F;
        f+= F * albedo_color / abs(wi.z);

      } else {
        bool entering = hit_record.front_face;
        float etaI = entering ? etaA : etaB;
        float etaT = entering ? etaB : etaA;
        glm::vec3 wi;
        glm::vec3 wo_forward = wo.z < 0 ? glm::vec3(0,0,-1) : glm::vec3(0,0,1);
        if (Refract(wo, wo_forward, etaI / etaT, &wi)) {
            glm::vec3 ft = albedo_color * (1 - F);
            //ft *= (etaI * etaI) / (etaT * etaT);
            *wiW = LocalToWorld(hit_record, wi);
            *pdf = 1 - F;
            f+=ft / abs(wi.z);
        }
      }
    }
    return f;
  }
  float pdf(HitRecord &hit_record,
                     glm::vec3 &woW, glm::vec3 &wiW) const{
    glm::vec3 wi = WorldToLocal(hit_record, wiW),
              wo = WorldToLocal(hit_record, woW);
    return local_pdf(wo, wi);

  }
  glm::vec3 f(HitRecord &hit_record, glm::vec3 &woW, glm::vec3 &wiW) const {
    glm::vec3 wi = WorldToLocal(hit_record, wiW),
              wo = WorldToLocal(hit_record, woW);
    bool reflect = dot(wiW, hit_record.geometry_normal) *
                       dot(woW, hit_record.geometry_normal) >
                   0;
    glm::vec3 f(0.0f);
    if (reflect && IsLambertian()) {
      f += local_f_Lambertian(wo, wi);
    }
    if (IsTransmissive()) {
      f += local_f_Transmissive(wo, wi);
    }
    if (reflect && IsMirror()) {
      f += local_f_Specular(wo, wi);
    }
    return f;
  }
};
}  // namespace sparks
