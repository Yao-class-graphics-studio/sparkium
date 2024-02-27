#pragma once
#include "cstdint"
#include "glm/glm.hpp"
#include "sparks/assets/hit_record.h"
#include "sparks/util/util.h"
#include "vector"
#include "cassert"
#include "iostream"

namespace sparks {

// BxDF (except BSSRDF): Shengquan Du.

inline float sqr(float x) {
  return x * x;
}

inline glm::vec3 SampleFromCosine(const glm::vec2 &args) {
  float t1 = args[0], t2 = args[1];
  float theta = acos(sqrt(1 - t1)), phi = 2 * PI * t2;
  glm::vec3 localRes(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
  return localRes;
}
inline bool SameHemiSphere(const glm::vec3 &wo, const glm::vec3 &wi) {
  return wo.z * wi.z > 0.0f;
}
inline float CosTheta(const glm::vec3 &w) {
  return w.z;
}
inline float AbsCosTheta(const glm::vec3 &w) {
  return std::fabs(w.z);
}
inline float Cos2Theta(const glm::vec3 &w) {
  return w.z * w.z;
}
inline float Sin2Theta(const glm::vec3 &w) {
  return std::max(0.0f, 1 - Cos2Theta(w));
}
inline float Tan2Theta(const glm::vec3 &w) {
  return Sin2Theta(w) / Cos2Theta(w);
}
inline float SinTheta(const glm::vec3 &w) {
  return std::sqrt(Sin2Theta(w));
}

inline float CosPhi(const glm::vec3 &w) {
  float sinTheta = SinTheta(w);
  return sinTheta == 0 ? 1 : glm::clamp(w.x / sinTheta, -1.0f, 1.0f);
}

inline float SinPhi(const glm::vec3 &w) {
  float sinTheta = SinTheta(w);
  return sinTheta == 0 ? 0 : glm::clamp(w.y / sinTheta, -1.0f, 1.0f);
}

inline float Cos2Phi(const glm::vec3 &w) {
  float sin2Theta = Sin2Theta(w);
  return sin2Theta == 0 ? 1 : std::min(w.x * w.x / sin2Theta, 1.0f);  
}
inline float Sin2Phi(const glm::vec3 &w) {
  float sin2Theta = Sin2Theta(w);
  return sin2Theta == 0 ? 0 : std::min(w.y * w.y / sin2Theta, 1.0f);
}

inline bool Refract(const glm::vec3 &wi,
                    const glm::vec3 &n,
                    float eta,
                    glm::vec3 &wt) {
  float cosThetaI = glm::dot(n, wi);
  float sin2ThetaI = std::max(0.0f, 1 - cosThetaI * cosThetaI);
  float sin2ThetaT = eta * eta * sin2ThetaI;

  if (sin2ThetaT >= 1)
    return false;
  float cosThetaT = std::sqrt(1 - sin2ThetaT);
  wt = eta * (-wi) + (eta * cosThetaI - cosThetaT) * n;
  return true;
}

inline glm::vec3 FrDielectric(float cosThetaI, float etaI, float etaT) {
  cosThetaI = glm::clamp(cosThetaI, -1.0f, 1.0f);
  bool entering = cosThetaI > 0.0f;
  if (!entering) {
    std::swap(etaI, etaT);
    cosThetaI = std::fabs(cosThetaI);
  }

  float sinThetaI = std::sqrt(std::max(0.0f, 1 - cosThetaI * cosThetaI));
  float sinThetaT = etaI / etaT * sinThetaI;
  if (sinThetaT >= 1)
    return glm::vec3{1.0f};
  float cosThetaT = std::sqrt(std::max(0.0f, 1 - sinThetaT * sinThetaT));
  float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                ((etaT * cosThetaI) + (etaI * cosThetaT));
  float Rperp = ((etaI * cosThetaI) - etaT * (cosThetaT)) /
                ((etaI * cosThetaI) + (etaT * cosThetaT));
  return glm::vec3((Rparl * Rparl + Rperp * Rperp) / 2);
}

// each BxDFType contains  one of reflection/transmission and one of diffuse/glossy/specular. Retro-reflection is treated as glossy
enum BxDFType : int {
  BSDF_REFLECTION = 1<<0,
  BSDF_TRANSMISSION = 1<<1,
  BSDF_DIFFUSE = 1<<2,
  BSDF_GLOSSY = 1<<3,
  BSDF_SPECULAR = 1<<4,
  BSDF_ALL = BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR,
};

//BxDF assuming normal = \hat{z} 
struct BxDF {
  BxDF() = default;
  virtual ~BxDF() {}
  explicit BxDF(const BxDFType &bxdf_type) : type_(bxdf_type) {
  }
  bool MatchesFlags(BxDFType t) const {
    return (type_ & t) == type_;
  }
  //calculate bxdf 
  virtual glm::vec3 f(const glm::vec3 &wo,
                      const glm::vec3 &wi)const=0;
  // Sample wi by bxdf and return the F result
  virtual glm::vec3 Sample_f(const glm::vec3 &wo, glm::vec3 &wi, const glm::vec2 &args, float &pdf) const;
  virtual float Pdf(const glm::vec3 &wo,
                    const glm::vec3 &wi) const;
  const BxDFType type_{BSDF_REFLECTION};
};

// scaled BxDF
class ScaledBxDF final : public BxDF {
 public:
  ScaledBxDF() = delete;
  ScaledBxDF(const BxDF *bxdf, const glm::vec3 &scale)
      : bxdf_(bxdf), scale_(scale) {
  }
  ~ScaledBxDF() {
    delete bxdf_;
  }
  // calculate bxdf
  glm::vec3 f(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    return scale_ * bxdf_->f(wo, wi);
  }
  // Sample wi by bxdf and return the F result
  glm::vec3 Sample_f(const glm::vec3 &wo,
                     glm::vec3 &wi,
                     const glm::vec2 &args,
                     float &pdf) const override {
    return scale_ * bxdf_->Sample_f(wo, wi, args, pdf);
  }
  float Pdf(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    return bxdf_->Pdf(wo, wi);
  }
 private:
  const BxDF *bxdf_;
  glm::vec3 scale_;
};

//a collection of BRDFs & BTDFs 
struct BSDF {
  BSDF() = delete;
  BSDF(const HitRecord &hit);
  void Add(BxDF *bxdf, float weight) {
    if (weight > 0.0f)
        bxdfs_.push_back(std::make_pair(bxdf, weight));
  }
  int NumComponents(BxDFType flags) const {
    int num = 0;
    for (auto [bxdf, weight] : bxdfs_) {
      if (bxdf->MatchesFlags(flags))
        ++num;
    }
    return num;
  }
  glm::vec3 LocalToWorld(const glm::vec3 &w) const {
    return local_to_world_ * w;
  }
  glm::vec3 WorldToLocal(const glm::vec3 &w) const {
    return glm::transpose(local_to_world_) * w;
  }

  glm::vec3 f(const glm::vec3 &woWorld,
              const glm::vec3 &wiWorld,
              BxDFType flags) const;
  //select bxdf according to weight, arg and type, and then remap arg to [0,1)
  BxDF *SelectBxDF(float &arg, BxDFType type, int &matchingComps) const;
  
  glm::vec3 Sample_f(const glm::vec3 &woWorld,
                     glm::vec3 &wiWorld,
                     const glm::vec2 &args,
                     float &pdf,
                     BxDFType type,
                     BxDFType &sampledType) const;
  float Pdf(const glm::vec3 &woWorld,
            const glm::vec3 &wiWorld,
            BxDFType flags = BxDFType(BSDF_ALL)) const;
  virtual ~BSDF() {
    for (auto [bxdf, weight] : bxdfs_)
      delete bxdf;
  }

  HitRecord hit_;
  glm::mat3 local_to_world_;
  std::vector<std::pair<BxDF *, float>> bxdfs_;
};

class MicrofacetDistribution {
 public:
  virtual ~MicrofacetDistribution() {}
  virtual float D(const glm::vec3 &wh) const = 0;
  virtual float Lambda(const glm::vec3 &w) const = 0;
  float G1(const glm::vec3 &w) const {
    return 1 / (1 + Lambda(w));
  }
  virtual float G(const glm::vec3 &wo, const glm::vec3 &wi) const {
    return 1 / (1 + Lambda(wo) + Lambda(wi));
  }
  virtual glm::vec3 Sample_wh(const glm::vec3 &wo,
                              const glm::vec2 &args) const = 0;
  float Pdf(const glm::vec3 &wo, const glm::vec3 &wh) const {
    assert(!std::isnan(wh[0]) && !std::isnan(wh[1]) && !std::isnan(wh[2]));
    if (sampleVisibleArea_)
      return D(wh) * G1(wo) * std::fabs(glm::dot(wo, wh)) / AbsCosTheta(wo);
    else
      return D(wh) * AbsCosTheta(wh);
  }

 protected:
  MicrofacetDistribution(bool sampleVisibleArea)
      : sampleVisibleArea_(sampleVisibleArea) {
  }
  const bool sampleVisibleArea_;
};

static void TrowbridgeReitzSample11(float cosTheta,
                                    float U1,
                                    float U2,
                                    float &slope_x,
                                    float &slope_y) {
  if (cosTheta > 1.0f - 1e-4f) {
    float r = std::sqrt(U1 / (1 - U1));
    float phi = 2 * PI * U2;
    slope_x = r * std::cos(phi);
    slope_y = r * std::sin(phi);
    return;
  }

  float sinTheta = std::sqrt(std::max(0.0f, 1 - cosTheta * cosTheta));
  float tanTheta = sinTheta / cosTheta;
  float a = 1 / tanTheta;
  float G1 = 2 / (1 + std::sqrt(1.0f + 1.0f / (a * a)));

  float A = 2 * U1 / G1 - 1;
  float tmp = std::min(1e10f, 1 / (A * A - 1));
  float B = tanTheta;
  float D =
      std::sqrt(std::max(B * B * tmp * tmp - (A * A - B * B) * tmp, 0.0f));
  float slope_x_1 = B * tmp - D;
  float slope_x_2 = B * tmp + D;
  slope_x = (A < 0 || slope_x_2 > 1 / tanTheta) ? slope_x_1 : slope_x_2;

  float S;
  if (U2 > 0.5f) {
    S = 1;
    U2 = 2 * (U2 - 0.5);
  } else {
    S = -1;
    U2 = 2 * (0.5 - U2);
  }
  float z = (U2 * (U2 * (U2 * 0.27385f - 0.73369f) + 0.46341f)) /
            (U2 * (U2 * (U2 * 0.093073f + 0.309420f) - 1.000000f) + 0.597999f);
  slope_y = S * z * std::sqrt(1.f + slope_x * slope_x);
  //std::cerr << slope_x << " " << slope_y << std::endl;
  assert(!std::isinf(slope_y));
  assert(!std::isnan(slope_y));
}

static glm::vec3 TrowbridgeReitzSample(const glm::vec3 &wi,
                                      float alpha_x,
                                      float alpha_y,
                                      float U1,
                                      float U2) {
  // 1. stretch wi
  glm::vec3 wiStretched =
      glm::normalize(glm::vec3(alpha_x * wi.x, alpha_y * wi.y, wi.z));

  // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
  float slope_x, slope_y;
  TrowbridgeReitzSample11(CosTheta(wiStretched), U1, U2, slope_x, slope_y);

  // 3. rotate
  float tmp = CosPhi(wiStretched) * slope_x - SinPhi(wiStretched) * slope_y;
  slope_y = SinPhi(wiStretched) * slope_x + CosPhi(wiStretched) * slope_y;
  slope_x = tmp;

  // 4. unstretch
  slope_x = alpha_x * slope_x;
  slope_y = alpha_y * slope_y;

  // 5. compute normal
  return glm::normalize(glm::vec3(-slope_x, -slope_y, 1.));
}

class TrowbridgeReitzDistribution : public MicrofacetDistribution {
 public:
  TrowbridgeReitzDistribution(float alphax, float alphay, bool sampleVis = true)
      : MicrofacetDistribution(sampleVis),
        alphax_(std::max(1e-3f, alphax)),
        alphay_(std::max(1e-3f, alphay)) {
  }
  float D(const glm::vec3 &wh) const override {
    float tan2Theta = Tan2Theta(wh);

    if (std::isinf(tan2Theta))
      return 0.0f;
    const float cos4Theta = sqr(Cos2Theta(wh));
    float e = (Cos2Phi(wh) / (alphax_ * alphax_) +
               Sin2Phi(wh) / (alphay_ * alphay_)) *
              tan2Theta;
    assert(!std::isnan(e));
    return 1 / (PI * alphax_ * alphay_ * cos4Theta * sqr(1 + e));
  }
  glm::vec3 Sample_wh(const glm::vec3 &wo,
                      const glm::vec2 &args) const override {
    glm::vec3 wh;
    if (!sampleVisibleArea_) {
      assert(0);
      float cosTheta = 0, phi = (2 * PI) * args[1];
      if (alphax_ == alphay_) {
        float tanTheta2 = alphax_ * alphax_ * args[0] / (1.0f - args[0]);
        cosTheta = 1 / std::sqrt(1 + tanTheta2);
      } else {
        float phi = std::atan(alphay_ / alphax_ * std::tan(2 * PI * args[1] + .5f * PI));
        if (args[1] > .5f)
          phi += PI;
        float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
        const float alphax2 = alphax_ * alphax_, alphay2 = alphay_ * alphay_;
        const float alpha2 =
            1 / (cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
        float tanTheta2 = alpha2 * args[0] / (1 - args[0]);
        cosTheta = 1 / std::sqrt(1 + tanTheta2);
      }
      float sinTheta =
          std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
      wh = glm::vec3(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
      if (!SameHemiSphere(wo, wh))
        wh = -wh;
    } else {
      bool flip = wo.z < 0;
      wh = TrowbridgeReitzSample(flip ? -wo : wo, alphax_, alphay_, args[0], args[1]);
      if (flip)
        wh = -wh;
    }
    return wh;
  }

 private:
  float Lambda(const glm::vec3 &w) const override {
    float absTanTheta = std::sqrt(Tan2Theta(w));
    if (std::isinf(absTanTheta))
      return 0;
    float alpha = std::sqrt(Cos2Phi(w) * alphax_ * alphax_ +
                            Sin2Phi(w) * alphay_ * alphay_);
    float alpha2Tan2Theta = sqr(alpha * absTanTheta);
    return (-1 + std::sqrt(1.0f + alpha2Tan2Theta)) / 2;
  }
  const float alphax_, alphay_;
};

class Fresnel {
 public:
  virtual ~Fresnel() {}
  virtual glm::vec3 Evaluate(float cosThetaI) const = 0;
};

class FresnelDielectric : public Fresnel {
 public:
  glm::vec3 Evaluate(float cosThetaI) const override;
  FresnelDielectric(float etaI, float etaT) : etaI_(etaI), etaT_(etaT) {}

 private:
  float etaI_, etaT_;
};

class MicrofacetReflection final : public BxDF {
 public:
  MicrofacetReflection(const glm::vec3 &R,
                       MicrofacetDistribution *distribution,
                       Fresnel *fresnel)
      : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
        R_(R),
        distribution_(distribution),
        fresnel_(fresnel) {
  }
  ~MicrofacetReflection() {
    delete distribution_;
    delete fresnel_;
  }
  glm::vec3 f(const glm::vec3 &wo, const glm::vec3 &wi) const override;
  glm::vec3 Sample_f(const glm::vec3 &wo,
                     glm::vec3 &wi,
                     const glm::vec2 &args,
                     float &pdf) const override;
  float Pdf(const glm::vec3 &wo, const glm::vec3 &wi) const override;

 private:
  const glm::vec3 R_;
  const MicrofacetDistribution *distribution_;
  const Fresnel *fresnel_;
};

class MicrofacetTransmission final : public BxDF {
 public:
  MicrofacetTransmission(const glm::vec3 &T,
                         MicrofacetDistribution *distribution,
                         float etaA,
                         float etaB)
      : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_GLOSSY)),
        T_(T),
        distribution_(distribution),
        etaA_(etaA),
        etaB_(etaB),
        fresnel_(etaA,etaB) {}
  glm::vec3 f(const glm::vec3 &wo, const glm::vec3 &wi) const override;
  glm::vec3 Sample_f(const glm::vec3 &wo,
                     glm::vec3 &wi,
                     const glm::vec2 &args,
                     float &pdf) const override;
  float Pdf(const glm::vec3 &wo, const glm::vec3 &wi) const override;
  ~MicrofacetTransmission() {
    delete distribution_;
  }

 private:
  const glm::vec3 T_;
  const MicrofacetDistribution *distribution_;
  const float etaA_, etaB_;
  const FresnelDielectric fresnel_;
};

class LambertianReflection : public BxDF {
 public:
  LambertianReflection(const glm::vec3 &R)
      : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R_(R) {
  }
  glm::vec3 f(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    return R_ * INV_PI;
  }
 private:
  const glm::vec3 R_;
};

class LambertianTransmission : public BxDF {
 public:
  LambertianTransmission(const glm::vec3 &T)
      : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_DIFFUSE)), T_(T) {
  }
  glm::vec3 f(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    return T_ * INV_PI;
  }
  glm::vec3 Sample_f(const glm::vec3 &wo,
                     glm::vec3 &wi,
                     const glm::vec2 &args,
                     float &pdf) const {
    wi = SampleFromCosine(args);
    if (wo.z > 0)
      wi.z *= -1;
    pdf = Pdf(wo, wi);
    return f(wo, wi);
  }
  float Pdf(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    return !SameHemiSphere(wo, wi) ? AbsCosTheta(wi) * INV_PI : 0;
  }

 private:
  glm::vec3 T_;
};

class SpecularReflection : public BxDF {
 public:
  ~SpecularReflection() {
    delete fresnel_;
  }
  SpecularReflection(const glm::vec3 &R, Fresnel *fresnel)
      : BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)),
        R_(R),
        fresnel_(fresnel) {
  }
  glm::vec3 f(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    return glm::vec3{0.0f};
  }
  glm::vec3 Sample_f(const glm::vec3 &wo, glm::vec3 &wi, const glm::vec2 &args, float& pdf) const override {
    wi = glm::vec3{-wo.x, -wo.y, wo.z};
    pdf = 1000.0f;
    return pdf * fresnel_->Evaluate(CosTheta(wi)) * R_ / AbsCosTheta(wi);
  }
  float Pdf(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    return 0;
  }
 private:
  glm::vec3 R_;
  Fresnel *fresnel_;
};
class SpecularTransmission : public BxDF {
 public:
  SpecularTransmission(const glm::vec3 &T, float etaA, float etaB)
      : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)),
        T_(T),
        etaA_(etaA),
        etaB_(etaB),
        fresnel_(etaA, etaB) {
  }
  glm::vec3 f(const glm::vec3 &wo, const glm::vec3 &wi) const {
    return glm::vec3{0.0f};
  }
  glm::vec3 Sample_f(const glm::vec3 &wo,
                     glm::vec3 &wi,
                     const glm::vec2 &args,
                     float &pdf) const override {
    bool entering = CosTheta(wo) > 0;
    float etaI = entering ? etaA_ : etaB_;
    float etaT = entering ? etaB_ : etaA_;
    if (!Refract(wo, glm::vec3(0.0f,0.0f,(wo.z<0.0f?-1.0f:1.0f)),etaI/etaT, wi))
        return glm::vec3{0.0f};
    pdf=1000.0f;
    glm::vec3 ft = T_ * (glm::vec3(1.0f) - fresnel_.Evaluate(CosTheta(wi)));
    ft *= (etaI * etaI) / (etaT * etaT);
    return pdf * ft / AbsCosTheta(wi);
  }
  float Pdf(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    return 0;
  }

 private:
  glm::vec3 T_;
  float etaA_, etaB_;
  FresnelDielectric fresnel_;
};

class AlphaTransmission : public BxDF {
 public:
  AlphaTransmission(const glm::vec3 &T)
      : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)),
        T_(T) {
  }
  glm::vec3 f(const glm::vec3 &wo, const glm::vec3 &wi) const {
    return glm::vec3{0.0f};
  }
  glm::vec3 Sample_f(const glm::vec3 &wo,
                     glm::vec3 &wi,
                     const glm::vec2 &args,
                     float &pdf) const override {
    wi = -wo;
    pdf = 1000.0f;
    return pdf * T_ / AbsCosTheta(wi);
  }
  float Pdf(const glm::vec3 &wo, const glm::vec3 &wi) const override {
    return 0;
  }

 private:
  glm::vec3 T_;
};
}  // namespace sparks
