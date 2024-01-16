#pragma once
#include "cstdint"
#include "glm/glm.hpp"
#include "sparks/assets/util.h"
#include "sparks/assets/material.h"
#include "sparks/assets/bxdf.h"

// reference:
// https://github.com/mmp/pbrt-v3/blob/master/src/materials/disney.cpp
//  https://www.pbr-book.org/3ed-2018/contents

namespace sparks {

// Disney principled BSDF: Shengquan Du.

class Scene;

inline float Pow5(float x) {
  return x * x * x * x * x;
}
inline float SchlickWeight(float cosTheta) {
  float m = glm::clamp(1 - cosTheta, 0.0f, 1.0f);
  return Pow5(m);
}
inline float FrSchlick(float R0, float cosTheta) {
  return glm::mix(R0, 1.0f, SchlickWeight(cosTheta));
}
inline glm::vec3 FrSchlick(const glm::vec3 &R0, float cosTheta) {
  return glm::mix(R0, glm::vec3(1.0f), SchlickWeight(cosTheta));
}
inline float SchlickR0FromEta(float eta) {
  return sqr(eta - 1) / sqr(eta + 1);
}

class DisneyDiffuse : public BxDF {
 public:
  DisneyDiffuse(const glm::vec3 &R)
      : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R_(R) {
  }
  glm::vec3 f(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    float Fo = SchlickWeight(AbsCosTheta(wo));
    float Fi = SchlickWeight(AbsCosTheta(wi));
    
    return R_ * INV_PI * (1 - Fo / 2) * (1 - Fi / 2);
  }

 private:
  glm::vec3 R_;
};

class DisneyFakeSS : public BxDF {
 public:
  DisneyFakeSS(const glm::vec3 &R, float roughness)
      : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)),
        R_(R),
        roughness_(roughness) {
  }
  glm::vec3 f(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    glm::vec3 wh = wi + wo;
    if (wh == glm::vec3{0.0f})
      return glm::vec3{0.0f};
    wh = glm::normalize(wh);
    float cosThetaD = glm::dot(wi, wh);

    float Fss90 = sqr(cosThetaD) * roughness_;
    float Fo = SchlickWeight(AbsCosTheta(wo));
    float Fi = SchlickWeight(AbsCosTheta(wi));
    float Fss = glm::mix(1.0f, Fss90, Fo) * glm::mix(1.0f, Fss90, Fi);
    float ss =
        1.25f * (Fss * (1 / (AbsCosTheta(wo) + AbsCosTheta(wi)) - 0.5f) + 0.5f);
    return R_ * INV_PI * ss;
  }

 private:
  glm::vec3 R_;
  float roughness_;
};

class DisneyRetro : public BxDF {
 public:
  DisneyRetro(const glm::vec3 &R, float roughness)
      : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)),
        R_(R),
        roughness_(roughness) {
  }
  glm::vec3 f(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    glm::vec3 wh = wi + wo;
    if (wh == glm::vec3{0.0f})
      return glm::vec3{0.0f};
    wh = glm::normalize(wh);
    float cosThetaD = glm::dot(wi, wh);

    float Fo = SchlickWeight(AbsCosTheta(wo));
    float Fi = SchlickWeight(AbsCosTheta(wi));
    float Rr = 2 * roughness_ * sqr(cosThetaD);
    return R_ * INV_PI * Rr * (Fo + Fi + Fo * Fi * (Rr - 1));
  }

 private:
  glm::vec3 R_;
  float roughness_;
};

class DisneySheen : public BxDF {
 public:
  DisneySheen(const glm::vec3 &R)
      : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R_(R) {
  }
  glm::vec3 f(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    glm::vec3 wh = wi + wo;
    if (wh == glm::vec3{0.0f})
      return glm::vec3{0.0f};
    wh = glm::normalize(wh);
    float cosThetaD = glm::dot(wi, wh);

    return R_ * SchlickWeight(cosThetaD);
  }

 private:
  glm::vec3 R_;
};

inline float GTR1(float cosTheta, float alpha) {
  float alpha2 = alpha * alpha;
  return (alpha2 - 1) /
         (PI * std::log(alpha2) * (1 + (alpha2 - 1) * cosTheta * cosTheta));
}

// Smith masking/shadowing term.
inline float smithG_GGX(float cosTheta, float alpha) {
  float alpha2 = alpha * alpha;
  float cosTheta2 = cosTheta * cosTheta;
  return 1 / (cosTheta + std::sqrt(alpha2 + cosTheta2 - alpha2 * cosTheta2));
}


class DisneyClearcoat : public BxDF {
 public:
  DisneyClearcoat(float weight, float gloss)
      : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
        weight_(weight),
        gloss_(gloss) {
  }
  glm::vec3 f(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    glm::vec3 wh = wi + wo;
    if (wh == glm::vec3{0.0f})
      return glm::vec3{0.0f};
    wh = glm::normalize(wh);

    float Dr = GTR1(AbsCosTheta(wh), gloss_);
    float Fr = FrSchlick(.04, glm::dot(wo, wh));
    float Gr =
        smithG_GGX(AbsCosTheta(wo), .25) * smithG_GGX(AbsCosTheta(wi), .25);

    return glm::vec3{weight_ * Gr * Fr * Dr / 4};
  }
  glm::vec3 Sample_f(const glm::vec3 &wo,
                     glm::vec3 &wi,
                     const glm::vec2 &args,
                     float &pdf) const override {
    if (wo.z == 0)
      return glm::vec3{0.0f};

    float alpha2 = gloss_ * gloss_;
    float cosTheta = std::sqrt(
        std::max(0.0f, (1 - std::pow(alpha2, 1 - args[0])) / (1 - alpha2)));
    float sinTheta = std::sqrt(std::max(0.0f, 1 - cosTheta * cosTheta));
    float phi = 2 * PI * args[1];
    glm::vec3 wh{sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta};
    if (!SameHemiSphere(wo, wh))
      wh = -wh;

    wi = -glm::reflect(wo, wh);
    if (!SameHemiSphere(wo, wi))
      return glm::vec3{0.0f};

    pdf = Pdf(wo, wi);
    return f(wo, wi);
  }
  float Pdf(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    if (!SameHemiSphere(wo, wi))
      return 0;

    glm::vec3 wh = wi + wo;
    if (wh == glm::vec3{0.0f})
      return 0;
    wh = glm::normalize(wh);

    // The sampling routine samples wh exactly from the GTR1 distribution.
    // Thus, the final value of the PDF is just the value of the
    // distribution for wh converted to a mesure with respect to the
    // surface normal.
    float Dr = GTR1(AbsCosTheta(wh), gloss_);
    return Dr * AbsCosTheta(wh) / (4 * glm::dot(wo, wh));
  }

 private:
  float weight_, gloss_;
};

class DisneyFresnel : public Fresnel {
 public:
  DisneyFresnel(const glm::vec3 &R0, float metallic, float eta)
      : R0_(R0), metallic_(metallic), eta_(eta) {
  }
  glm::vec3 Evaluate(float cosI) const override {
    return glm::mix(FrDielectric(cosI, 1, eta_), FrSchlick(R0_, cosI), metallic_);
  }

 private:
  const glm::vec3 R0_;
  const float metallic_, eta_;
};

class DisneyMicrofacetDistribution : public TrowbridgeReitzDistribution {
 public:
  DisneyMicrofacetDistribution(float alphax, float alphay)
      : TrowbridgeReitzDistribution(alphax, alphay) {}
  float G(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    return G1(wo) * G1(wi);
  }
};



}  // namespace sparks
