#pragma once
#include "sparks/assets/bxdf.h"
#include "sparks/assets/medium.h"
#include <memory>
#include <random>
#include <functional>

// all w are local direction
// BSSRDF: Jingyi Lyu
namespace sparks {

class Scene;
class Material;
typedef std::function<float(const glm::vec3&, const glm::vec3&, float, float, HitRecord*)> TraceMethod;
// typedef std::function<float(const glm::vec3&, const glm::vec3&, float, float, HitRecord*)> TraceMethod;

enum BSSRDFType: int {
    BSSRDF_TABULATED = 1 << 0,
    BSSRDF_DISNEY = 1 << 1
};

class BSSRDF {
public:
    BSSRDFType type;
    BSSRDF() = default;
    virtual ~BSSRDF() {}
    BSSRDF(const glm::vec3 po_, const glm::vec3 wo_, float eta_, BSSRDFType type_): 
           po(po_), wo(wo_), eta(eta_), type(type_) {}
    virtual glm::vec3 S(const glm::vec3 pi, const glm::vec3 wi) const = 0;
    virtual glm::vec3 Sample_S(const TraceMethod &trace, glm::vec3 &pi, glm::vec3 &wi, glm::vec3 &wiNormal, float &rawPdf, float &fullPdf, 
                               const int entity, std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const = 0;
protected:
    glm::vec3 po, wo;
    float eta;
};

class SeparableBSSRDF: public BSSRDF {
public:
    SeparableBSSRDF(const glm::vec3 po_, const glm::vec3 wo_, float eta_, BSSRDFType type_, 
                    const glm::vec3 normal_, const Material *material_):
                    BSSRDF(po_, wo_, eta_, type_), normal(normal_), material(material_) {
        tx = normal == glm::vec3(0.0f, 0.0f, 1.0f) ? 
             glm::normalize(glm::cross(normal, glm::vec3(0.0f, 1.0f, 0.0f))) :
             glm::normalize(glm::cross(normal, glm::vec3(0.0f, 0.0f, 1.0f)));
        ty = glm::normalize(glm::cross(normal, tx));
        local2World = glm::mat3{tx, ty, normal};
        world2Local = glm::inverse(local2World);
    }
    glm::vec3 S(const glm::vec3 pi, const glm::vec3 wi) const override;
    glm::vec3 Sw(const glm::vec3 wi) const;
    glm::vec3 Sp(const glm::vec3 pi) const;
    glm::vec3 Sample_S(const TraceMethod &trace, glm::vec3 &pi, glm::vec3 &wi, glm::vec3 &wiNormal, float &rawPdf, float &fullPdf, 
                       const int entity, std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const;
    glm::vec3 Sample_Sp(const TraceMethod &trace, glm::vec3 &pi, glm::vec3 &wiNormal, float &pdf, 
                        const int entity, std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const;
    float Pdf_Sp(const HitRecord &pi) const;
    virtual glm::vec3 Sr(float r) const = 0;
    virtual float Sample_Sr(int ch, float distance) const = 0;
    virtual float Pdf_Sr(int ch, float sample) const = 0;
protected:
    glm::vec3 normal, tx, ty;
    glm::mat3 local2World, world2Local;
    const Material *material;
};

class BSSRDFTable {
   public:
    const int nRhoSamples, nRadiusSamples;
    float *rhoSamples, *radiusSamples;
    float *profile;
    float *rhoEff;
    float *profileCDF;

    BSSRDFTable() : nRhoSamples(0), nRadiusSamples(0) {}
    BSSRDFTable(int nRhoSamples, int nRadiusSamples)
        : nRhoSamples(nRhoSamples),
          nRadiusSamples(nRadiusSamples),
          rhoSamples(new float[nRhoSamples]),
          radiusSamples(new float[nRadiusSamples]),
          profile(new float[nRadiusSamples * nRhoSamples]),
          rhoEff(new float[nRhoSamples]),
          profileCDF(new float[nRadiusSamples * nRhoSamples]) {
    }
    ~BSSRDFTable() {
        delete[] rhoSamples;
        delete[] radiusSamples;
        delete[] profile;
        delete[] rhoEff;
        delete[] profileCDF;
    }
    // BSSRDFTable(const BSSRDFTable&) = default;
    // BSSRDFTable& operator = (const BSSRDFTable&) = default;
    float EvalProfile(int rhoIndex, int radiusIndex) const {
        return profile[rhoIndex * nRadiusSamples + radiusIndex];
    }
};

class TabulatedBSSRDF: public SeparableBSSRDF {
    glm::vec3 sigma_t, rho;
    const BSSRDFTable &table;
public:
    TabulatedBSSRDF(const glm::vec3 po_, const glm::vec3 wo_, float eta_, BSSRDFType type_, 
                    const glm::vec3 normal_, const Material *material_,
                    const glm::vec3 a_, const glm::vec3 s_, const BSSRDFTable &table_):
                    SeparableBSSRDF(po_, wo_, eta_, type_, normal_, material_), 
                    sigma_t(a_ + s_), table(table_) {
                        for(int i = 0; i <= 2; i++)
                        rho[i] = sigma_t[i] != 0 ? s_[i] / sigma_t[i] : 0;
                    }
    glm::vec3 Sr(float r) const override;
    float Sample_Sr(int ch, float distance) const override;
    float Pdf_Sr(int ch, float sample) const override;
};

class DisneyBSSRDF : public SeparableBSSRDF {
   public:
    DisneyBSSRDF(const glm::vec3 &R,
                 const glm::vec3 &d,
                 const glm::vec3 &po,
                 const glm::vec3 &wo, 
                 const glm::vec3 &normal,
                 float eta,
                 const Material *material)
        : SeparableBSSRDF(po,
                          wo,
                          eta,
                          BSSRDFType(BSSRDF_DISNEY),
                          normal,
                          material),
          R_(R),
          d_(glm::vec3{0.2f} * d) {
    }
    //glm::vec3 S(const glm::vec3 pi, const glm::vec3 wi) const override;
    glm::vec3 Sr(float d) const override;
    float Sample_Sr(int ch, float u) const override;
    float Pdf_Sr(int ch, float r) const override;
   private:
    glm::vec3 R_, d_;
};

void ComputeBeamDiffusionBSSRDF(float g, float eta, BSSRDFTable &table);
void SubsurfaceFromDiffuse(const BSSRDFTable &table, const glm::vec3 rhoEff, const glm::vec3 &mfp, 
                           glm::vec3 &sigma_a, glm::vec3 &sigma_s);

}