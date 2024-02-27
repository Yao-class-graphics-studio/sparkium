#include "grassland/grassland.h"
#include "sparks/assets/bxdf.h"
#include "grassland/util/util.h"

namespace sparks {

// BxDF (except BSSRDF): Shengquan Du.

// BxDF method definition
glm::vec3 BxDF::Sample_f(const glm::vec3 &wo, glm::vec3 &wi, const glm::vec2 &args, float& pdf) const{
  wi = SampleFromCosine(args);
  if (wo.z < 0.0f)
    wi.z *= -1;
  pdf = AbsCosTheta(wi) * INV_PI;
  return f(wo, wi);
}
float BxDF::Pdf(const glm::vec3 &wo, const glm::vec3 &wi) const {
  return SameHemiSphere(wo,wi) ? AbsCosTheta(wi) * INV_PI : 0.0f;
}

// compute fresnel term for dielectric material 


glm::vec3 FresnelDielectric::Evaluate(float cosThetaI) const {
  return glm::vec3(FrDielectric(cosThetaI, etaI_, etaT_));
}

// MicrofacetReflection method definition
glm::vec3 MicrofacetReflection::f(const glm::vec3 &wo,
                                  const glm::vec3 &wi) const {
  float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
  if (cosThetaO == 0.0f || cosThetaI == 0.0f)
    return glm::vec3{0.0f};
  glm::vec3 wh = wi + wo;
  if (wh == glm::vec3(0.0f))
    return glm::vec3{0.0f};
  wh = glm::normalize(wh);

  glm::vec3 F = fresnel_->Evaluate(glm::dot(wi, (wh.z < 0.0f ? -wh : wh)));
  return R_ * distribution_->D(wh) * distribution_->G(wo, wi) * F /
         (4 * cosThetaI * cosThetaO);
}

glm::vec3 MicrofacetReflection::Sample_f(const glm::vec3 &wo,
                                         glm::vec3 &wi,
                                         const glm::vec2 &args,
                                         float &pdf) const {
  if (wo.z == 0)
    return glm::vec3{0.0f};
  glm::vec3 wh = distribution_->Sample_wh(wo, args);
  if (glm::dot(wo, wh) < 0)
    return glm::vec3{0.0f};
  wi = -glm::reflect(wo, wh);
  if (!SameHemiSphere(wo, wi))
    return glm::vec3{0.0f};

  pdf = distribution_->Pdf(wo, wh) / (4 * glm::dot(wo, wh));
  return f(wo, wi);
}

float MicrofacetReflection::Pdf(const glm::vec3 &wo,
                                const glm::vec3 &wi) const {
    if (!SameHemiSphere(wo, wi))
        return 0.0f;
    glm::vec3 wh = glm::normalize(wo + wi);
    return distribution_->Pdf(wo, wh) / (4 * glm::dot(wo, wh));
}

// MicrofacetTransmission method definition
glm::vec3 MicrofacetTransmission::f(const glm::vec3 &wo, const glm::vec3 &wi) const {
    if (SameHemiSphere(wo, wi))
        return glm::vec3{0.0f};
    float cosThetaO = CosTheta(wo);
    float cosThetaI = CosTheta(wi);
    if (cosThetaI == 0 || cosThetaO == 0)
        return glm::vec3{0.0f};
    float eta = cosThetaO > 0 ? (etaB_ / etaA_) : (etaA_ / etaB_);
    glm::vec3 wh = wo + wi * eta;
    if (glm::length(wh) > 1e-6f)
        wh = glm::normalize(wh);
    else
        wh = glm::vec3{0.0f, 0.0f, 1.0f};
    if (wh.z < 0)
        wh = -wh;
    float dotOH = glm::dot(wo, wh), dotIH = glm::dot(wi,wh);

    if (dotOH * dotIH > 0.0f)
        return glm::vec3{0.0f};
    glm::vec3 F = fresnel_.Evaluate(dotOH);
    float sqrtDenom = dotOH + eta * dotIH;
    float factor = 1 / eta;

    return (glm::vec3(1.0f) - F) * T_ *
           std::fabs(distribution_->D(wh) * distribution_->G(wo, wi) * eta *
                     eta * dotIH * dotOH * factor * factor /
                     (cosThetaI * cosThetaO * sqrtDenom * sqrtDenom));
}



glm::vec3 MicrofacetTransmission::Sample_f(const glm::vec3 &wo,
                                           glm::vec3 &wi,
                                           const glm::vec2 &args,
                                           float &pdf) const {
    if (wo.z == 0)
        return glm::vec3{0.0f};
    glm::vec3 wh = distribution_->Sample_wh(wo, args);
    if (glm::dot(wo, wh) < 0)
        return glm::vec3{0.0f};
    float eta = CosTheta(wo) > 0 ? (etaA_ / etaB_) : (etaB_ / etaA_);
    if (!Refract(wo, wh, eta, wi))
        return glm::vec3{0.0f};
    pdf = Pdf(wo, wi);
    assert(!std::isnan(pdf));
    assert(!std::isinf(pdf));
    return f(wo, wi);
}

float MicrofacetTransmission::Pdf(const glm::vec3 &wo,
                                  const glm::vec3 &wi) const {
    if (SameHemiSphere(wo, wi))
        return 0;
    float eta = CosTheta(wo) > 0 ? (etaB_ / etaA_) : (etaA_ / etaB_);
    glm::vec3 wh = wo + wi * eta;
    if (glm::length(wh) > 1e-6f)
        wh = glm::normalize(wh);
    else
        wh = glm::vec3{0.0f, 0.0f, 1.0f};
    assert(!std::isnan(-wh[0]) && !std::isnan(wh[1]) && !std::isnan(wh[2]));
    if (glm::dot(wo, wh) * glm::dot(wi, wh) > 0)
        return 0;
    float sqrtDenom = glm::dot(wo, wh) + eta * glm::dot(wi, wh);
    float dwh_dwi =
        std::fabs(eta * eta * glm::dot(wi, wh) / (sqrtDenom * sqrtDenom));
    return distribution_->Pdf(wo, wh) * dwh_dwi;
}

// BSDF method definition
BSDF::BSDF(const HitRecord &hit) : hit_(hit) {
    
  if (!hit_.front_face) {
      hit_.geometry_normal *= -1.0f;
      hit_.normal *= -1.0f;
      hit_.tangent *= -1.0f;
  }
  glm::vec3 norm = hit_.normal;
  if (std::fabs(norm.z) > 1.0f - 1e-6) {
    local_to_world_ = glm::mat3(1.0f);
  } else {
        glm::vec3 y = glm::normalize(glm::cross(norm, glm::vec3(0.0f,0.0f,1.0f)));
        glm::vec3 x = glm::normalize(glm::cross(y, norm));
        local_to_world_ = glm::mat3(x,y,norm);
  }
}

glm::vec3 BSDF::f(const glm::vec3 &woWorld,
        const glm::vec3 &wiWorld,
        BxDFType flags) const {
  glm::vec3 wo = WorldToLocal(woWorld), wi = WorldToLocal(wiWorld);
  if (wo.z == 0.0f)
        return glm::vec3{0.0f};
  bool reflect = glm::dot(wiWorld, hit_.geometry_normal) *
                     glm::dot(woWorld, hit_.geometry_normal) >
                 0.0f;
  glm::vec3 f(0.0f);
  for (auto [bxdf, weight] : bxdfs_) {
        if (bxdf -> MatchesFlags(flags) &&
            ((reflect && (bxdf->type_ & BSDF_REFLECTION)) ||
            (!reflect && (bxdf -> type_ & BSDF_TRANSMISSION)))) {
          f += bxdf->f(wo, wi);
        }
  }
  return f;
}

BxDF* BSDF::SelectBxDF(float& arg, BxDFType type, int& matchingComps) const {
  std::vector<std::pair<BxDF*, float> > candidates;
  float preweight = 0.0f;
  matchingComps = 0;
  for (auto [bxdf, weight] : bxdfs_) {
    if (bxdf->MatchesFlags(type)) {
          candidates.emplace_back(bxdf, preweight);
          preweight += weight;
          ++matchingComps;
    }
  }
  if (candidates.size() == 0)
    return nullptr;
  candidates.emplace_back(nullptr, preweight);
  arg *= preweight;
  auto it = --std::upper_bound(candidates.begin(), candidates.end(), arg,
                         [](const float &a, const std::pair<BxDF *, float> &b) { return a < b.second; });
  arg = std::min((arg - it->second) / ((it + 1)->second - it->second), 1.0f-1e-9f);
  return it->first;
}

glm::vec3 BSDF::Sample_f(const glm::vec3 &woWorld,
               glm::vec3 &wiWorld,
               const glm::vec2 &args,
               float &pdf,
               BxDFType type,
               BxDFType &sampledType) const {
  int matchingComps = 0;
  glm::vec2 argsRemapped = args;
  BxDF *chosenBxdf = SelectBxDF(argsRemapped[0], type, matchingComps);
  if (chosenBxdf == nullptr)
  {
    pdf = 0.0f;
    sampledType = BxDFType(0);
    return glm::vec3{0.0f};
  }
  glm::vec3 wo = WorldToLocal(woWorld), wi;
  if (wo.z == 0.0f)
    return glm::vec3{0.0f};
  pdf = 0.0f;
  sampledType = chosenBxdf->type_;
  glm::vec3 f = chosenBxdf->Sample_f(wo, wi, argsRemapped, pdf);
  if (pdf == 0.0f) {
    sampledType = BxDFType(0);
    return glm::vec3{0.0f};
  }
  wiWorld = LocalToWorld(wi);

  float matchingWeight = 0.0f;
  for (auto [bxdf, weight] : bxdfs_) {
    if (bxdf->MatchesFlags(type)) {
          matchingWeight += weight;
          if (bxdf == chosenBxdf)
            pdf *= weight;
    }
  }

  if (!(chosenBxdf->type_ & BSDF_SPECULAR) && matchingComps > 1) {
    bool reflect = dot(wiWorld, hit_.geometry_normal) *
                       glm::dot(woWorld, hit_.geometry_normal) >
                   0;
    for (auto [bxdf, weight] : bxdfs_) {
          if (bxdf->MatchesFlags(type)) {
            if (bxdf != chosenBxdf)
              pdf += weight * bxdf->Pdf(wo, wi);

            if ((reflect && (bxdf->type_ & BSDF_REFLECTION)) ||
                (!reflect && (bxdf->type_ & BSDF_TRANSMISSION))) {
                if (bxdf != chosenBxdf)
                  f += bxdf->f(wo, wi);
            }
          }
    }
  }
  pdf /= matchingWeight;
  return f;
}

float BSDF::Pdf(const glm::vec3 &woWorld,
                const glm::vec3 &wiWorld,
                BxDFType flags) const {
  glm::vec3 wo = WorldToLocal(woWorld), wi = WorldToLocal(wiWorld);
  if (wo.z == 0.0f)
    return 0.0f;
  float pdf = 0.0f;
  float matchingWeight = 0.0f;
  for (auto [bxdf, weight] : bxdfs_) {
    if (bxdf->MatchesFlags(flags)) {
          matchingWeight += weight;
          pdf += weight * bxdf->Pdf(wo, wi);
    }
  }
  return matchingWeight > 1e-6f ? pdf / matchingWeight : 0.0f;
}


}  // namespace sparks