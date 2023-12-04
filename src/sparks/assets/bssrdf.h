#pragma once
#include "sparks/assets/bxdf.h"
#include "sparks/assets/scene.h"
#include "sparks/assets/medium.h"
#include <memory>
#include <random>

// all w are local direction

namespace sparks {

enum BSSRDFType: int {
    BSSRDF_TABULATED = 1 << 0,
};

class BSSRDF {
public:
    BSSRDFType type;
    BSSRDF() = default;
    virtual ~BSSRDF() {}
    BSSRDF(const glm::vec3 po_, const glm::vec3 wo_, float eta_, BSSRDFType type_): 
           po(po_), wo(wo_), eta(eta_), type(type_) {}
    virtual glm::vec3 S(const glm::vec3 pi, const glm::vec3 wi) const = 0;
    virtual glm::vec3 Sample_S(const Scene *scene, glm::vec3 &pi, glm::vec3 &wi, float &rawPdf, float &fullPdf, 
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
    }
    glm::vec3 S(const glm::vec3 pi, const glm::vec3 wi) const override;
    glm::vec3 Sw(const glm::vec3 wi) const;
    glm::vec3 Sp(const glm::vec3 pi) const;
    glm::vec3 Sample_S(const Scene *scene, glm::vec3 &pi, glm::vec3 &wi, float &rawPdf, float &fullPdf, 
                       const int entity, std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const;
    glm::vec3 Sample_Sp(const Scene *scene, glm::vec3 &pi, float &pdf, 
                        const int entity, std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const;
    float Pdf_Sp(const HitRecord &pi) const;
    virtual glm::vec3 Sr(float r) const = 0;
    virtual float Sample_Sr(int ch, float distance) const = 0;
    virtual float Pdf_Sr(int ch, float sample) const = 0;
protected:
    glm::vec3 normal, tx, ty;
    const Material *material;
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

class BSSRDFTable {
public:
    const int nRhoSamples, nRadiusSamples;
    std::unique_ptr<float[]> rhoSamples, radiusSamples;
    std::unique_ptr<float[]> profile;
    std::unique_ptr<float[]> rhoEff;
    std::unique_ptr<float[]> profileCDF;

    BSSRDFTable(): nRhoSamples(0), nRadiusSamples(0) {}
    BSSRDFTable(int nRhoSamples, int nRadiusSamples)
    : nRhoSamples(nRhoSamples),
      nRadiusSamples(nRadiusSamples),
      rhoSamples(new float[nRhoSamples]),
      radiusSamples(new float[nRadiusSamples]),
      profile(new float[nRadiusSamples * nRhoSamples]),
      rhoEff(new float[nRhoSamples]),
      profileCDF(new float[nRadiusSamples * nRhoSamples]) {}
    float EvalProfile(int rhoIndex, int radiusIndex) const {
        return profile[rhoIndex * nRadiusSamples + radiusIndex];
    }
};

void ComputeBeamDiffusionBSSRDF(float g, float eta, BSSRDFTable &table);

}