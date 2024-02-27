#include "sparks/assets/medium.h"

namespace sparks {

// Volumetric medium: Jingyi Lyu
// medium sampling utility

bool IntersectWithGrid(const glm::vec3 origin, const glm::vec3 direction, float &tMin, float &tMax) {
    tMin = 0, tMax = 1e10f;
    for(int i = 0; i <= 2; i++) {
        if(direction[i] == 0) {
            if(origin[i] >= 0 && origin[i] <= 1)
                continue;
            // std::cerr << "?" << std::endl;
            return false;    
        }
        float x = -origin[i] / direction[i];
        float y = (1 - origin[i]) / direction[i];
        tMax = std::min(tMax, std::max(x, y));
        tMin = std::max(tMin, std::min(x, y));
    }
    // std::cerr << tMax << std::endl;
    return tMin <= tMax;
}

float PhaseHG(float cosTheta, float g) {
    float den = 1 + g * g + 2 * g * cosTheta;
    return (1 - g * g) / (4 * PI * den * glm::sqrt(den));
}

float lerp(float x, float f0, float f1) {
    return (1 - x) * f0 + x * f1;
}

glm::vec3 HomogeneousMedium::Tr(glm::vec3 origin, glm::vec3 direction, float tMax,
                                std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const {
    return glm::exp(-sigma_t * std::min(tMax, MAX_FLOAT));
}

glm::vec3 HomogeneousMedium::Sample(const glm::vec3 origin, const glm::vec3 direction, glm::vec3 &sample, float &pdf, 
                                    bool &insideMedium, float tMax,
                                    std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const {
    int ch = std::min(2, (int)(uniform(rd) * 3));
    float d = -glm::log(1 - uniform(rd)) / sigma_t[ch];
    float t = std::min(d, tMax);
    insideMedium = d < tMax;
    if(insideMedium) {
        sample = origin + d * direction;
    }
    glm::vec3 tr;
    for(int i = 0; i <= 2; i++)
        tr[i] = glm::exp(-sigma_t[i] * std::min(t, MAX_FLOAT));
    glm::vec3 density = insideMedium ? sigma_t * tr : tr;
    pdf = 0;
    for(int i = 0; i <= 2; i++)
        pdf += density[i];
    pdf /= 3;
    return insideMedium ? tr * sigma_s / pdf : tr / pdf;
}

float HomogeneousMedium::Sample_p(const glm::vec3 &wo, glm::vec3 &wi, 
                                  std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const {
    float cosTheta;
    if(std::fabs(g) < 1e-3)
        cosTheta = 1 - 2 * uniform(rd);
    else {
        float sqrTerm = (1 - g * g) / (1 + g - 2 * g * uniform(rd));
        cosTheta = (1 + g * g - sqrTerm * sqrTerm) / (2 * g);
    }
    float sinTheta = glm::sqrt(std::max(0.0f, 1 - cosTheta * cosTheta));
    float phi = 2 * PI * uniform(rd);
    float cosPhi = glm::cos(phi), sinPhi = glm::sqrt(1 - cosPhi * cosPhi);
    glm::vec3 localRes = glm::vec3(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
    glm::vec3 v1 = wo == glm::vec3(0.0f, 0.0f, 1.0f) ? 
                   glm::normalize(glm::cross(wo, glm::vec3(0.0f, 1.0f, 0.0f))) :
                   glm::normalize(glm::cross(wo, glm::vec3(0.0f, 0.0f, 1.0f)));
    glm::vec3 v2 = glm::normalize(glm::cross(wo, v1));
    glm::mat3 local2World{v1, v2, wo};
    wi = local2World * wo;
    return PhaseHG(cosTheta, g);
}

float GridDensityMedium::Sample_p(const glm::vec3 &wo, glm::vec3 &wi, 
                                  std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const {
    float cosTheta;
    if(std::fabs(g) < 1e-3)
        cosTheta = 1 - 2 * uniform(rd);
    else {
        float sqrTerm = (1 - g * g) / (1 + g - 2 * g * uniform(rd));
        cosTheta = (1 + g * g - sqrTerm * sqrTerm) / (2 * g);
    }
    float sinTheta = glm::sqrt(std::max(0.0f, 1 - cosTheta * cosTheta));
    float phi = 2 * PI * uniform(rd);
    float cosPhi = glm::cos(phi), sinPhi = glm::sqrt(1 - cosPhi * cosPhi);
    glm::vec3 localRes = glm::vec3(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
    glm::vec3 v1 = wo == glm::vec3(0.0f, 0.0f, 1.0f) ? 
                   glm::normalize(glm::cross(wo, glm::vec3(0.0f, 1.0f, 0.0f))) :
                   glm::normalize(glm::cross(wo, glm::vec3(0.0f, 0.0f, 1.0f)));
    glm::vec3 v2 = glm::normalize(glm::cross(wo, v1));
    glm::mat3 local2World{v1, v2, wo};
    wi = local2World * wo;
    return PhaseHG(cosTheta, g);
}

glm::vec3 GridDensityMedium::Tr(glm::vec3 origin, glm::vec3 direction, float tMax,
                                std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const {
    glm::vec3 localOrigin = glm::vec3(world2Medium * glm::vec4{origin, 1.0f});
    glm::vec3 localDirection = glm::vec3(scale * glm::vec4{direction, 1.0f});
    // glm::vec3 tmp = glm::vec3(world2Medium * glm::vec4(0, 0, 0, 1.0f));
    // std::cerr << localOrigin[0] << std::endl;
    // if(origin[1] >= 248.7f && origin[1] <= 248.8f)
    //     std::cerr << localOrigin[0] << std::endl;
    float gridMin, gridMax;
    if(!IntersectWithGrid(localOrigin, localDirection, gridMin, gridMax))
        return glm::vec3{1.0f};
    float t = gridMin, res = 1.0f;
    gridMax = std::min(gridMax, tMax);
    while(true) {
        t -= glm::log(1 - uniform(rd)) * invMaxDensity / sigma_t;
        if(t >= gridMax)
            break;
        float nowDensity = Density(localOrigin + t * localDirection);
        // std::cerr << invMaxDensity << std::endl;
        res *= 1 - std::max(0.0f, nowDensity * invMaxDensity);
    }
    // std::cerr << res << std::endl;
    return glm::vec3{res};
}

glm::vec3 GridDensityMedium::Sample(const glm::vec3 origin, const glm::vec3 direction, glm::vec3 &sample, float &pdf, 
                                    bool &insideMedium, float tMax,
                                    std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const {
    glm::vec3 localOrigin = glm::vec3(world2Medium * glm::vec4{origin, 1.0f});
    glm::vec3 localDirection = glm::vec3(scale * glm::vec4{direction, 1.0f});
    float gridMin, gridMax;
    if(!IntersectWithGrid(localOrigin, localDirection, gridMin, gridMax))
        return glm::vec3{1.0f};
    float t = gridMin;
    gridMax = std::min(gridMax, tMax);
    while(true) {
        t -= glm::log(1 - uniform(rd)) * invMaxDensity / sigma_t;
        if(t >= gridMax)
            break;
        if(Density(localOrigin + t * localDirection) * invMaxDensity > uniform(rd)) {
            sample = origin + t * direction;
            insideMedium = true;
            return sigma_s / sigma_t;
        }
    }
    return glm::vec3{1.0f};
}

float GridDensityMedium::Density(const glm::vec3 &p) const {
    glm::vec3 pSamples(p[0] * nx - 0.5f, p[1] * ny - 0.5f, p[2] * nz - 0.5f);
    glm::ivec3 pi = glm::ivec3((int)pSamples[0], (int)pSamples[1], (int)pSamples[2]);
    glm::vec3 d = pSamples - glm::vec3(pi[0], pi[1], pi[2]);
    float d00 = lerp(d[0], D(pi), D(pi + glm::ivec3(1, 0, 0)));
    float d01 = lerp(d[0], D(pi + glm::ivec3(0, 1, 0)), D(pi + glm::ivec3(1, 1, 0)));
    float d10 = lerp(d[0], D(pi + glm::ivec3(0, 0, 1)), D(pi + glm::ivec3(1, 0, 1)));
    float d11 = lerp(d[0], D(pi + glm::ivec3(0, 1, 1)), D(pi + glm::ivec3(1, 1, 1)));
    float d0 = lerp(d[1], d00, d10);
    float d1 = lerp(d[1], d01, d11);
    return lerp(d[2], d0, d1);
}

}