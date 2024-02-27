#pragma once
#include "sparks/assets/model.h"
#include "glm/glm.hpp"
#include <memory>
#include <cstring>
#include <cassert>
#include <random>

namespace sparks {

// Volumetric medium: Jingyi Lyu.

bool GetMediumScatteringProperties(const std::string &name,
                                   glm::vec3 *sigma_a,
                                   glm::vec3 *sigma_prime_s);

float PhaseHG(float cosTheta, float g);

class Medium {
public:
    virtual ~Medium() {}
    virtual glm::vec3 Tr(glm::vec3 origin, glm::vec3 direction, float tMax,
                         std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const = 0;
    virtual glm::vec3 Sample(const glm::vec3 origin, const glm::vec3 direction, glm::vec3 &sample, float &pdf, 
                             bool &insideMedium, float tMax,
                             std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const = 0;
    virtual float Sample_p(const glm::vec3 &wo, glm::vec3 &wi, 
                           std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const = 0;
    virtual float p(glm::vec3 wo, glm::vec3 wi) const = 0;
    virtual glm::vec3 getEmission() const = 0;
    virtual void Update(glm::vec3 a_, glm::vec3 s_, glm::vec3 emission_, float g_) = 0;
protected:
    const float MAX_FLOAT = 1e10f;
};

class HomogeneousMedium: public Medium {
    glm::vec3 sigma_a, sigma_s, sigma_t;
    glm::vec3 emission;
    float g;
public:
    HomogeneousMedium(glm::vec3 a_, glm::vec3 s_, glm::vec3 emission_, float g_):
        sigma_a(a_), sigma_s(s_), sigma_t(a_ + s_), emission(emission_), g(g_) {};
    glm::vec3 Tr(glm::vec3 origin, glm::vec3 direction, float tMax,
                 std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const;
    glm::vec3 Sample(const glm::vec3 origin, const glm::vec3 direction, glm::vec3 &sample, float &pdf, 
                             bool &insideMedium, float tMax,
                             std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const;
    float Sample_p(const glm::vec3 &wo, glm::vec3 &wi, 
                   std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const;
    float p(glm::vec3 wo, glm::vec3 wi) const {
        return PhaseHG(glm::dot(wo, wi), g);
    }
    glm::vec3 getEmission() const {
        return emission;
    }
    void Update(glm::vec3 a_, glm::vec3 s_, glm::vec3 emission_, float g_) {
        sigma_a = a_, sigma_s = s_, emission = emission_, g = g_;
        sigma_t = sigma_a + sigma_s;
    }
};

class GridDensityMedium: public Medium {
    glm::vec3 sigma_a, sigma_s;
    glm::vec3 emission;
    float g;
    int nx, ny, nz;
    glm::mat4 world2Medium, scale;
    float* density;
    float sigma_t;
    float invMaxDensity;
public:
    GridDensityMedium(glm::vec3 a_, glm::vec3 s_, glm::vec3 emission_, float g_, 
                      int nx_, int ny_, int nz_, glm::mat4 transform_, glm::mat4 linear_, float *d_):
        sigma_a(a_), sigma_s(s_), emission(emission_), g(g_), nx(nx_), ny(ny_), nz(nz_), 
        world2Medium(transform_), scale(linear_), density(new float[nx * ny * nz]) {
            memcpy(density, d_, sizeof(float) * nx * ny * nz);
            sigma_t = (sigma_a + sigma_s)[0];
            float maxDensity = 0;
            for(int i = 0; i < nx * ny * nz; i++) {
                // std::cerr << density[i] << std::endl;
                maxDensity = std::max(maxDensity, density[i]);
            }
            assert(maxDensity > 0);
            invMaxDensity = 1 / maxDensity;
        }
    ~GridDensityMedium() {
        delete[] density;
    }
    glm::vec3 Tr(glm::vec3 origin, glm::vec3 direction, float tMax,
                 std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const;
    glm::vec3 Sample(const glm::vec3 origin, const glm::vec3 direction, glm::vec3 &sample, float &pdf, 
                             bool &insideMedium, float tMax,
                             std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const;
    float Sample_p(const glm::vec3 &wo, glm::vec3 &wi, 
                           std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const;
    float Density(const glm::vec3 &p) const;
    float D(const glm::ivec3 &p) const {
        if(p[0] < 0 || p[0] >= nx || p[1] < 0 || p[1] >= ny || p[2] < 0 || p[2] >= nz)
            return 0.0f;
        return density[p[2] * nx * ny + p[1] * nx + p[0]];
    }
    float p(glm::vec3 wo, glm::vec3 wi) const {
        return PhaseHG(glm::dot(wo, wi), g);
    }
    glm::vec3 getEmission() const {
        return emission;
    }
    void Update(glm::vec3 a_, glm::vec3 s_, glm::vec3 emission_, float g_) {
        sigma_a = a_, sigma_s = s_, emission = emission_, g = g_;
        sigma_t = (sigma_a + sigma_s)[0];
    }
};


}