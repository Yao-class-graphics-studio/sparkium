#include "sparks/assets/bssrdf.h"
#include <vector>
#include "glm/gtc/type_ptr.hpp"

namespace sparks {

// BSSRDF: Jingyi Lyu.

// BSSRDF Utility

template <typename Predicate>
int FindInterval(int n, const Predicate &f) {
  int first = 0, len = n;
    while(len > 0) {
        int half = len >> 1, middle = first + half;
        if (f(middle)) {
            first = middle + 1;
            len -= half + 1;
        } else
            len = half;
    }
    return std::min(std::max(first - 1, 0), n - 2);
}

float FrFirstMoment(float eta) {
    float eta2 = eta * eta, eta3 = eta2 * eta, eta4 = eta2 * eta2, eta5 = eta2 * eta3;
    if(eta < 1)
        return 0.45966f - 1.73965f * eta + 3.37668f * eta2 - 3.904945 * eta3 +
               2.49277f * eta4 - 0.68441f * eta5;
    else
        return -4.61686f + 11.1136f * eta - 10.4646f * eta2 + 5.11455f * eta3 -
               1.27198f * eta4 + 0.12746f * eta5;
}

float FrSecondMoment(float eta) {
    float eta2 = eta * eta, eta3 = eta2 * eta, eta4 = eta3 * eta, eta5 = eta4 * eta;
    if (eta < 1) {
        return 0.27614f - 0.87350f * eta + 1.12077f * eta2 - 0.65095f * eta3 +
               0.07883f * eta4 + 0.04860f * eta5;
    } else {
        float r_eta = 1 / eta, r_eta2 = r_eta * r_eta, r_eta3 = r_eta2 * r_eta;
        return -547.033f + 45.3087f * r_eta3 - 218.725f * r_eta2 +
               458.843f * r_eta + 404.557f * eta - 189.519f * eta2 +
               54.9327f * eta3 - 9.00603f * eta4 + 0.63942f * eta5;
    }
}


bool CatmullRomWeights(int size, const float *nodes, float x, int &offset, float *weights) {
    if(x < nodes[0] || x > nodes[size - 1])
        return false;
    int idx = FindInterval(size, [&](int i) { return nodes[i] <= x; });
    offset = idx - 1;
    float x0 = nodes[idx], x1 = nodes[idx + 1];
    float t = (x - x0) / (x1 - x0), t2 = t * t, t3 = t2 * t;
    float w0, w3;
    weights[0] = weights[3] = 0.0f;
    weights[1] = 2 * t3 - 3 * t2 + 1;
    weights[2] = -2 * t3 + 3 * t2;
    if(idx > 0) {
        w0 = (t3 - 2 * t2 + t) * (x1 - x0) / (x1 - nodes[idx - 1]);
        weights[0] = -w0, weights[2] += w0;
    } else {
        w0 = t3 - 2 * t2 + t;
        weights[0] = 0, weights[1] -= w0, weights[2] += w0;
    }
    if(idx + 2 < size) {
        w3 = (t3 - t2) * (x1 - x0) / (nodes[idx + 2] - x0);
        weights[1] -= w3;
        weights[3] = w3;
    } else {
        w3 = t3 - t2;
        weights[3] = 0.0f, weights[1] -= w3, weights[2] += w3;
    }
    return true;
}

float SampleCatmullRom2D(int size1, int size2, const float *nodes1, const float *nodes2, 
                         const float *values, const float *cdf, float alpha, float u) {
    int offset;
    float weights[4];
    if(!CatmullRomWeights(size1, nodes1, alpha, offset, weights))
        return 0.0f;
    auto interpolate = [&](const float *array, int idx) {
        float res = 0;
        for(int i = 0; i <= 3; i++)
            if(weights[i] != 0)
                res += array[(offset + i) * size2 + idx] * weights[i];
        return res;
    };
    float maximum = interpolate(cdf, size2 - 1);
    u *= maximum;
    int idx = FindInterval(size2, [&](int i) { return interpolate(cdf, i) <= u; });
    float f0 = interpolate(values, idx), f1 = interpolate(values, idx + 1);
    float x0 = nodes2[idx], x1 = nodes2[idx + 1];
    float w = x1 - x0;
    float d0, d1;
    u = (u - interpolate(cdf, idx)) / w;
    if(idx > 0)
        d0 = w * (f1 - interpolate(values, idx - 1)) / (x1 - nodes2[idx - 1]);
    else
        d0 = f1 - f0;
    if(idx < size2 - 2)
        d1 = w * (interpolate(values, idx + 2) - f0) / (nodes2[idx + 2] - x0);
    else
        d1 = f1 - f0;
    float t;
    if(f0 != f1)
        t = (f0 - glm::sqrt(std::max(0.0f, f0 * f0 + 2 * u * (f1 - f0)))) / (f0 - f1);
    else
        t = u / f0;
    float a = 0, b = 1, Fhat, fhat;
    while(true) {
        if(!(t >= a && t <= b)) 
            t = 0.5f * (a + b);
        Fhat = t * (f0 + t * (.5f * d0 + t * ((1.f / 3.f) * (-2 * d0 - d1) + f1 - f0 + t * (.25f * (d0 + d1) + .5f * (f0 - f1)))));
        fhat = f0 + t * (d0 + t * (-2 * d0 - d1 + 3 * (f1 - f0) + t * (d0 + d1 + 2 * (f0 - f1))));
        if(std::abs(Fhat - u) < 1e-6f || b - a < 1e-6f) 
            break;
        if(Fhat - u < 0)
            a = t;
        else
            b = t;
        t -= (Fhat - u) / fhat;
    }
    return x0 + w * t;
}

float IntegrateCatmullRom(int n, const float *x, const float *values, float *cdf) {
    float sum = 0;
    cdf[0] = 0;
    for(int i = 0; i < n - 1; i++) {
        float x0 = x[i], x1 = x[i + 1];
        float f0 = values[i], f1 = values[i + 1];
        float w = x1 - x0;
        float d0, d1;
        if(i > 0)
            d0 = w * (f1 - values[i - 1]) / (x1 - x[i - 1]);
        else
            d0 = f1 - f0;
        if(i + 2 < n)
            d1 = w * (values[i + 2] - f0) / (x[i + 2] - x0);
        else
            d1 = f1 - f0;
        sum += ((d0 - d1) * (1.0f / 12.0f) + (f0 + f1) * 0.5f) * w;
        cdf[i + 1] = sum;
    }
    return sum;
}

float InvertCatmullRom(int n, const float *x, const float *values, float u) {
    if (!(u > values[0]))
        return x[0];
    else if (!(u < values[n - 1]))
        return x[n - 1];
    int i = FindInterval(n, [&](int i) { return values[i] <= u; });
    float x0 = x[i], x1 = x[i + 1];
    float f0 = values[i], f1 = values[i + 1];
    float w = x1 - x0;
    float d0, d1;
    if (i > 0)
        d0 = w * (f1 - values[i - 1]) / (x1 - x[i - 1]);
    else
        d0 = f1 - f0;
    if (i + 2 < n)
        d1 = w * (values[i + 2] - f0) / (x[i + 2] - x0);
    else
        d1 = f1 - f0;
    float a = 0, b = 1, t = 0.5f;
    float Fhat, fhat;
    while(true) {
        if (!(t > a && t < b)) 
            t = 0.5f * (a + b);
        float t2 = t * t, t3 = t2 * t;
        Fhat = (2 * t3 - 3 * t2 + 1) * f0 + (-2 * t3 + 3 * t2) * f1 + (t3 - 2 * t2 + t) * d0 + (t3 - t2) * d1;
        fhat = (6 * t2 - 6 * t) * f0 + (-6 * t2 + 6 * t) * f1 + (3 * t2 - 4 * t + 1) * d0 + (3 * t2 - 2 * t) * d1;
        if (std::fabs(Fhat - u) < 1e-6f || b - a < 1e-6f) 
            break;
        if (Fhat - u < 0)
            a = t;
        else
            b = t;
        t -= (Fhat - u) / fhat;
    }
    return x0 + t * w;
}

// float PhaseHG(float cosTheta, float g) {
//     float den = 1 + g * g + 2 * g * cosTheta;
//     return (1 - g * g) / (4 * PI * den * glm::sqrt(den));
// }

float BeamDiffusionMS(float sigma_s, float sigma_a, float g, float eta, float r) {
    const int nSamples = 50;
    float Ed = 0.0f;
    float sigmap_s = sigma_s * (1 - g);
    float sigmap_t = sigma_a + sigmap_s;
    float rho_p = sigmap_s / sigmap_t;
    float Dg = (2 * sigma_a + sigmap_s) / (3 * sigmap_t * sigmap_t);
    float sigma_tr = glm::sqrt(sigma_a / Dg);
    float fm1 = FrFirstMoment(eta), fm2 = FrSecondMoment(eta);
    float ze = -2 * Dg * (1 + 3 * fm2) / (1 - 2 * fm1);
    float cPhi = 0.25f * (1 - 2 * fm1), cE = 0.5f * (1 - 3 * fm2);
    for(int i = 0; i < nSamples; i++) {
        float zr = -std::log(1 - (i + 0.5f) / nSamples) / sigmap_t;
        float zv = -zr + 2 * ze;
        float dr = glm::sqrt(r * r + zr * zr), dv = glm::sqrt(r * r + zv * zv);
        float phiD = INV_PI / (4 * Dg) * (glm::exp(-sigma_tr * dr) / dr - std::exp(-sigma_tr * dv) / dv);
        float EDn = INV_PI / 4 * (zr * (1 + sigma_tr * dr) * glm::exp(-sigma_tr * dr) / (dr * dr * dr) - 
                                  zv * (1 + sigma_tr * dv) * glm::exp(-sigma_tr * dv) / (dv * dv * dv));
        float E = phiD * cPhi + EDn * cE;
        float k = 1 - glm::exp(-2 * sigmap_t * (dr + zr));
        Ed += k * rho_p * rho_p * E; 
    }
    return Ed / nSamples;
}

float BeamDiffusionSS(float sigma_s, float sigma_a, float g, float eta, float r) {
    const int nSamples = 50;
    float sigma_t = sigma_a + sigma_s, rho = sigma_s / sigma_t;
    float tCrit = r * glm::sqrt(eta * eta - 1);
    float Ess = 0;
    for(int i = 0; i < nSamples; i++) {
        float ti = tCrit - glm::log(1 - (i + 0.5f) / nSamples) / sigma_t;
        float d = glm::sqrt(r * r + ti * ti);
        float cosTheta0 = ti / d;
        Ess += rho * glm::exp(-sigma_t * (d + tCrit)) / (d * d) * 
               PhaseHG(cosTheta0, g) * (1 - FrDielectric(-cosTheta0, 1, eta)[0]) * 
               std::fabs(cosTheta0);
    }
    return Ess / nSamples;
}

void ComputeBeamDiffusionBSSRDF(float g, float eta, BSSRDFTable &table) {
    table.radiusSamples[0] = 0;
    table.radiusSamples[1] = 2.5e-3f;
    for(int i = 2; i < table.nRadiusSamples; i++) 
        table.radiusSamples[i] = table.radiusSamples[i - 1] * 1.2f;
    for(int i = 0; i < table.nRhoSamples; i++)
        table.rhoSamples[i] = (1 - glm::exp(-8.0f * i / (table.nRhoSamples - 1))) / (1 - glm::exp(-8.0f));
    for(int i = 0; i < table.nRhoSamples; i++) {
        for(int j = 0; j < table.nRadiusSamples; j++) {
            float rho = table.rhoSamples[i], r = table.radiusSamples[j];
            table.profile[i * table.nRadiusSamples + j] = 2 * PI * r * 
                                                          (BeamDiffusionSS(rho, 1 - rho, g, eta, r) + 
                                                           BeamDiffusionMS(rho, 1 - rho, g, eta, r));
            // std::cerr << table.profile[i * table.nRadiusSamples + j] << std::endl;
            table.rhoEff[i] = IntegrateCatmullRom(table.nRadiusSamples, table.radiusSamples, 
                                                  &table.profile[i * table.nRadiusSamples], &table.profileCDF[i * table.nRadiusSamples]);
        }
    }
}

void SubsurfaceFromDiffuse(const BSSRDFTable &table, const glm::vec3 rhoEff, const glm::vec3 &mfp, glm::vec3 &sigma_a, glm::vec3 &sigma_s) {
    for(int ch = 0; ch <= 2; ch++) {
        float rho = InvertCatmullRom(table.nRhoSamples, table.rhoSamples, table.rhoEff, rhoEff[ch]);
        sigma_s[ch] = rho / mfp[ch];
        sigma_a[ch] = (1 - rho) / mfp[ch];
    }
}

// BSSRDF implementation

glm::vec3 SeparableBSSRDF::S(const glm::vec3 pi, const glm::vec3 wi) const {
    return (glm::vec3{1.0f} - FrDielectric(glm::dot(wo, normal), 1, eta)) * Sp(pi) * Sw(wi);
}

glm::vec3 SeparableBSSRDF::Sw(const glm::vec3 wi) const {
    glm::vec3 wwi = world2Local * wi;
    float c = 1 - 2 * FrFirstMoment(1 / eta);
    return (glm::vec3{1.0f} - FrDielectric(CosTheta(wwi), 1, eta)) / (c * PI);
}

glm::vec3 SeparableBSSRDF::Sp(const glm::vec3 pi) const {
    return Sr(glm::length(pi - po));
}

glm::vec3 SeparableBSSRDF::Sample_S(const TraceMethod &trace, glm::vec3 &pi, glm::vec3 &wi, glm::vec3 &wiNormal, float &rawPdf, float &fullPdf,  
                                    const int entity, std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const {
    glm::vec3 sp = Sample_Sp(trace, pi, wiNormal, rawPdf, entity, rd, uniform);
    if(sp != glm::vec3{0.0f}) {
        glm::vec2 args(uniform(rd), uniform(rd));
        wi = SampleFromCosine(args);
        fullPdf = rawPdf * SinTheta(wi) * AbsCosTheta(wi) * INV_PI;
        // std::cerr << rawPdf << " " << fullPdf << std::endl;
        glm::vec3 sx = wiNormal == glm::vec3(0.0f, 0.0f, 1.0f) ? 
                       glm::normalize(glm::cross(wiNormal, glm::vec3(0.0f, 1.0f, 0.0f))) :
                       glm::normalize(glm::cross(wiNormal, glm::vec3(0.0f, 0.0f, 1.0f)));
        glm::vec3 sy = glm::normalize(glm::cross(wiNormal, sx));
        glm::mat3 wiLocal2World = glm::mat3{sx, sy, wiNormal};
        //printf("dir before transform: %f %f %f\n", wi[0], wi[1], wi[2]);
        //printf("normal: %f %f %f\n", wiNormal[0], wiNormal[1], wiNormal[2]);
        wi = glm::normalize(wiLocal2World * wi);
        //printf("addr; dir after transform: %p %f %f %f\n", glm::value_ptr(wi), wi[0], wi[1], wi[2]);
    } 
    return sp;
}

glm::vec3 SeparableBSSRDF::Sample_Sp(const TraceMethod &trace, glm::vec3 &pi, glm::vec3 &wiNormal, float &pdf, 
                                     const int entity, std::mt19937 &rd, std::uniform_real_distribution<float> &uniform) const {
    float u1 = uniform(rd);
    glm::vec3 vx, vy, vz;
    if(u1 < 0.5f)
        vx = tx, vy = ty, vz = normal;
    else if(u1 < 0.75f)
        vx = ty, vy = normal, vz = tx;
    else
        vx = normal, vy = tx, vz = ty;
    int ch = std::min(2, (int)(uniform(rd) * 3));
    float r = Sample_Sr(ch, uniform(rd));
    if(r < 0)
        return glm::vec3{0.0f};
    float phi = 2 * PI * uniform(rd);
    float rMax = Sample_Sr(ch, 0.999f);
    if(r > rMax)
        return glm::vec3{0.0f};
    float l = 2 * glm::sqrt(rMax * rMax - r * r);
    std::vector<HitRecord> isects;
    glm::vec3 origin = po + r * (vx * glm::cos(phi) + vy * glm::sin(phi)) - l * vz / 2.0f;
    // std::cerr << r << std::endl;
    while(glm::distance(origin, po) <= rMax + 0.0001f) {
        HitRecord currentHit;
        float t = trace(origin, vz, 0.0001f, 1e10f, &currentHit);
        if(t < 0)
            break;
        if(currentHit.hit_entity_id == entity) {
            isects.emplace_back(currentHit);
            // if(currentHit.position[1] == 165.0f)
            //     std::cerr << origin[1] << std::endl;
        }
        origin = currentHit.position + vz * 0.0001f;
    }
    int isectCnt = isects.size();
    if(isectCnt == 0)
        return glm::vec3{0.0f};
    int selected = std::min(isectCnt - 1, (int)(uniform(rd) * isectCnt));
    pi = isects[selected].position;
    wiNormal = isects[selected].geometry_normal;
    if(!isects[selected].front_face)
        wiNormal *= -1.0f;
    pdf = Pdf_Sp(isects[selected]) / isectCnt;
    return Sp(isects[selected].position);
}

float SeparableBSSRDF::Pdf_Sp(const HitRecord &pi) const {
    glm::vec3 d = pi.position - po;
    glm::vec3 dLocal = glm::vec3(glm::dot(tx, d), glm::dot(ty, d), glm::dot(normal, d));
    glm::vec3 nLocal = glm::vec3(glm::dot(tx, pi.geometry_normal), 
                                 glm::dot(ty, pi.geometry_normal),
                                 glm::dot(normal, pi.geometry_normal));
    float rProj[] = {glm::sqrt(dLocal[1] * dLocal[1] + dLocal[2] * dLocal[2]),
                     glm::sqrt(dLocal[2] * dLocal[2] + dLocal[0] * dLocal[0]),
                     glm::sqrt(dLocal[0] * dLocal[0] + dLocal[1] * dLocal[1])};
    float res = 0, axisProb[] = {0.25f, 0.25f, 0.5f};
    for(int axis = 0; axis <= 2; axis++) 
        for(int ch = 0; ch <= 2; ch++)
            res += Pdf_Sr(ch, rProj[axis]) * std::fabs(nLocal[axis]) * axisProb[axis] / 3.0f;
    return res;
}

glm::vec3 TabulatedBSSRDF::Sr(const float r) const {
    glm::vec3 sr{0.0f};
    for(int ch = 0; ch <= 2; ch++) {
        float rOptical = r * sigma_t[ch];
        int rhoOffset, radiusOffset;
        float rhoWeights[4], radiusWeights[4];
        bool flag1 = CatmullRomWeights(table.nRhoSamples, table.rhoSamples, rho[ch], rhoOffset, rhoWeights);
        bool flag2 = CatmullRomWeights(table.nRadiusSamples, table.radiusSamples, rOptical, radiusOffset, radiusWeights);
        if(!flag1 || !flag2)
            continue;
        float res = 0;
        for(int i = 0; i < 4; i++)
            for(int j = 0; j < 4; j++) {
                float w = rhoWeights[i] * radiusWeights[j];
                // std::cerr << w << std::endl;
                if(w != 0)
                    res += w * table.EvalProfile(rhoOffset + i, radiusOffset + j);
            }
        if(rOptical != 0.0f)
            res /= 2 * PI * rOptical;
        sr[ch] = std::max(res * sigma_t[ch] * sigma_t[ch], 0.0f);
    }
    return sr;
}

float TabulatedBSSRDF::Sample_Sr(int ch, float u) const {
    if(sigma_t[ch] == 0)
        return -1;
    return SampleCatmullRom2D(table.nRhoSamples, table.nRadiusSamples, table.rhoSamples, table.radiusSamples,
                              table.profile, table.profileCDF, rho[ch], u) / sigma_t[ch];
}

float TabulatedBSSRDF::Pdf_Sr(int ch, float r) const {
    float rOptical = r * sigma_t[ch];
    int rhoOffset, radiusOffset;
    float rhoWeights[4], radiusWeights[4];
    bool flag1 = CatmullRomWeights(table.nRhoSamples, table.rhoSamples, rho[ch], rhoOffset, rhoWeights);
    bool flag2 = CatmullRomWeights(table.nRadiusSamples, table.radiusSamples, rOptical, radiusOffset, radiusWeights);
    if(!flag1 || !flag2)
        return 0.0f;
    float sr = 0, rhoEff = 0;
    for(int i = 0; i <= 3; i++) {
        if(rhoWeights[i] == 0)
            continue;
        rhoEff += table.rhoEff[rhoOffset + i] * rhoWeights[i];
        for(int j = 0; j <= 3; j++)
            sr += rhoWeights[i] * radiusWeights[j] * table.EvalProfile(rhoOffset + i, radiusOffset + j);
    }
    return std::max(sr * sigma_t[ch] * sigma_t[ch] / rhoEff, 0.0f);
}

//glm::vec3 DisneyBSSRDF::S(const glm::vec3 pi, const glm::vec3 wi) const {
//    glm::vec3 a = glm::normalize(pi - po);
//    float fade = 1;
//    glm::vec3 n = normal;
//    float cosTheta = glm::dot(a, n);
//    if (cosTheta > 0) {
//        float sinTheta = std::sqrt(std::max(0.0f, 1 - cosTheta * cosTheta));
//        glm::vec3 a2 = n * sinTheta - (a - n * cosTheta) * cosTheta / sinTheta;
//        fade = std::max(0.0f, glm::dot())
//    }
//}

glm::vec3 DisneyBSSRDF::Sr(float r) const {
    if (r < 1e-6f)
        r = 1e-6f;
    auto f = [r](float d) {
      if (d == 0.0f)
        return 0.0f;
      return (std::exp(-r / d) + std::exp(-r / (3 * d))) / (8 * PI * d * r);
    };
    return R_ *
        glm::vec3(f(d_[0]), f(d_[1]), f(d_[2]));
}

float DisneyBSSRDF::Sample_Sr(int ch, float u) const {
    if (u < 0.25f) {
        u = std::min(u * 4, 1 - 1e-6f);
        return d_[ch] * std::log(1 / (1 - u));
    } else {
        u = std::min((u - 0.25f) / 0.75f, 1-1e-6f);
        return 3 * d_[ch] * std::log(1 / (1 - u));
    }
}

float DisneyBSSRDF::Pdf_Sr(int ch, float r) const {
    if (r < 1e-6f)
        r = 1e-6f;
    if (d_[ch] == 0.0f)
        return 0.0f;
    return (0.25f * std::exp(-r / d_[ch]) / (2 * PI * d_[ch] * r) +
            0.75 * std::exp(-r / (3 * d_[ch])) / (6 * PI * d_[ch] * r));
}

}