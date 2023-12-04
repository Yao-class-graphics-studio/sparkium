#pragma once

#include "sparks/assets/material.h"
#include "sparks/assets/bssrdf.h"

namespace sparks {

class SubsurfaceMaterial: public Material {
    BSSRDFTable table;
public:
    bool use_ss_texture = false;
    int sigma_a_texture_id = 0;
    int sigma_s_texture_id = 0;
    Medium *medium;
    glm::vec3 sigma_a, sigma_s;
    SubsurfaceMaterial() = default;
    explicit SubsurfaceMaterial(const glm::vec3 &albedo, const glm::vec3 a_, const glm::vec3 s_);
    SubsurfaceMaterial(Scene *scene, const tinyxml2::XMLElement *material_element);
    BSSRDF* ComputeBSSRDF(const HitRecord &hit, const glm::vec3 direction, const Scene* scene) const;
};

}