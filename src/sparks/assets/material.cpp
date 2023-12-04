#include "sparks/assets/material.h"

#include "grassland/grassland.h"
#include "sparks/assets/scene.h"
#include "sparks/assets/texture.h"
#include "sparks/util/util.h"
#include "sparks/assets/bxdf.h"
#include "sparks/assets/disney_material.h"

namespace sparks {

namespace {
std::unordered_map<std::string, MaterialType> material_name_map{
    {"lambertian", MATERIAL_TYPE_LAMBERTIAN},
    {"specular", MATERIAL_TYPE_SPECULAR},
    {"transmissive", MATERIAL_TYPE_TRANSMISSIVE},
    {"principled", MATERIAL_TYPE_PRINCIPLED},
    {"emission", MATERIAL_TYPE_EMISSION},
    {"subsurface", MATERIAL_TYPE_SUBSURFACE},
    {"kdsubsurface", MATERIAL_TYPE_KDSUBSURFACE},
    {"medium", MATERIAL_TYPE_MEDIUM}};
}

Material::Material(Scene *scene, const tinyxml2::XMLElement *material_element)
    : table(100, 64)  {
  if (!material_element) {
    return;
  }

  albedo_color = glm::vec3{1.0f};

  auto child_element = material_element->FirstChildElement("albedo");
  if (child_element) {
    albedo_color = StringToVec3(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("albedo_texture");
  if (child_element) {
    std::string path = child_element->FindAttribute("value")->Value();
    Texture albedo_texture(1, 1);
    if (Texture::Load(path, albedo_texture)) {
      albedo_texture_id =
          scene->AddTexture(albedo_texture, PathToFilename(path));
    }
  }

  child_element = material_element->FirstChildElement("normal_texture");
  if (child_element) {
    std::string path = child_element->FindAttribute("value")->Value();
    Texture normal_texture(1, 1);
    if (Texture::Load(path, normal_texture)) {
      use_normal_texture = true;
      normal_texture_id =
          scene->AddTexture(normal_texture, PathToFilename(path));
    }
  }

  child_element = material_element->FirstChildElement("emission");
  if (child_element) {
    emission = StringToVec3(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("emission_strength");
  if (child_element) {
    emission_strength =
        std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("alpha");
  if (child_element) {
    alpha = std::stof(child_element->FindAttribute("value")->Value());
  }

    child_element = material_element->FirstChildElement("metallic");
  if (child_element) {
    metallic = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("eta");
  if (child_element) {
    eta = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("roughness");
  if (child_element) {
    roughness = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("specularTint");
  if (child_element) {
    specularTint = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("anisotropic");
  if (child_element) {
    anisotropic = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("sheen");
  if (child_element) {
    sheen = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("sheenTint");
  if (child_element) {
    sheenTint = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("clearcoat");
  if (child_element) {
    clearcoat = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("clearcoatGloss");
  if (child_element) {
    clearcoatGloss = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("specTrans");
  if (child_element) {
    specTrans = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("flatness");
  if (child_element) {
    flatness = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("diffTrans");
  if (child_element) {
    diffTrans = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("scatterDistance");
  if (child_element) {
    scatterDistance =
        StringToVec3(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("thin");
  if (child_element) {
    thin = StringToBool(child_element->FindAttribute("value")->Value());
  }

  float g = 0.0f;

  child_element = material_element->FirstChildElement("sigma_a");
  if (child_element) {
    sigma_a = StringToVec3(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("sigma_s");
  if (child_element) {
    sigma_s = StringToVec3(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("mfp");
  if (child_element) {
    mfp = StringToVec3(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("sigma_a_texture");
  if (child_element) {
    std::string path = child_element->FindAttribute("value")->Value();
    Texture sigma_a_texture(1, 1);
    if (Texture::Load(path, sigma_a_texture)) {
      sigma_a_texture =
          scene->AddTexture(sigma_a_texture, PathToFilename(path));
    }
    use_ss_texture = true;
  }

  child_element = material_element->FirstChildElement("sigma_s_texture");
  if (child_element) {
    std::string path = child_element->FindAttribute("value")->Value();
    Texture sigma_a_texture(1, 1);
    if (Texture::Load(path, sigma_a_texture)) {
      sigma_a_texture =
          scene->AddTexture(sigma_a_texture, PathToFilename(path));
    }
    use_ss_texture = true;
  }

  child_element = material_element->FirstChildElement("g");
  if (child_element) {
    g = std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("volumetric_emission");
  if (child_element) {
    volumetric_emission = StringToVec3(child_element->FindAttribute("value")->Value());
  }

  material_type =
      material_name_map[material_element->FindAttribute("type")->Value()];

  if (material_type == MATERIAL_TYPE_SUBSURFACE || material_type == MATERIAL_TYPE_KDSUBSURFACE) {
    ComputeBeamDiffusionBSSRDF(g, eta, table);
  }

  if (material_type == MATERIAL_TYPE_MEDIUM) {
    medium = new HomogeneousMedium(sigma_a, sigma_s, volumetric_emission, g);
  }
}

Material::Material(const glm::vec3 &albedo) {
  albedo_color = albedo;
}

BSDF* Material::ComputeBSDF(const HitRecord &hit,
                         const Scene* scene) const {
  glm::vec3 color = albedo_color;
  //if (albedo_texture_id != -1)
    color *=
        glm::vec3(scene->GetTexture(albedo_texture_id).Sample(hit.tex_coord));
  HitRecord textureHit = hit;
  if (use_normal_texture)
    textureHit.normal =
        glm::vec3(scene->GetTexture(normal_texture_id).Sample(hit.tex_coord));

  BSDF* bsdf = new BSDF(textureHit);
  switch (material_type) { 
    case MATERIAL_TYPE_LAMBERTIAN: {
      bsdf->Add(new LambertianReflection(color), 1);
      break;
    }
    case MATERIAL_TYPE_EMISSION:{
      bsdf->Add(new LambertianReflection(color),1);
      break;
    }
    case MATERIAL_TYPE_SPECULAR: {
      bsdf->Add(new SpecularReflection(glm::vec3{1.0f},new FresnelDielectric(1,1e9)),1);
      break;
    }
    case MATERIAL_TYPE_TRANSMISSIVE: {
      bsdf->Add(new SpecularTransmission(glm::vec3{1.0f}, 1.0f, 1.5f), 1);
      break;
    }
    case MATERIAL_TYPE_PRINCIPLED: {
      // Diffuse
      float diffuseWeight = (1 - metallic) * (1 - specTrans);
      float dt = diffTrans / 2; // 0: all diffuse is reflected -> 1, transmitted
      float lum = 0.212671f * color[0] + 0.715160f * color[1] + 0.072169f * color[2];
      // normalize lum. to isolate hue+sat
      glm::vec3 Ctint = lum > 0 ? (color / lum) : glm::vec3{1.0f} ;

      glm::vec3 Csheen;
      if (sheen > 0) {
        Csheen = glm::mix(glm::vec3{1.0f}, Ctint, sheenTint);        
      }

      if (diffuseWeight > 0) {
        if (thin) {
          // Blend between DisneyDiffuse and fake subsurface based on
          // flatness.  Additionally, weight using diffTrans.
          bsdf->Add(new DisneyDiffuse(diffuseWeight * (1 - flatness) *
                                      (1 - dt) * color),
                    diffuseWeight * (1 - flatness) * (1 - dt));
          bsdf->Add(new DisneyFakeSS(
                        diffuseWeight * flatness * (1 - dt) * color, roughness),
                    diffuseWeight * flatness * (1 - dt));
        } else {
            //sd = scatterDistance
          if (scatterDistance==glm::vec3{0.0f})
            // No subsurface scattering; use regular (Fresnel modified)
            // diffuse.
            bsdf->Add(new DisneyDiffuse(diffuseWeight * color), diffuseWeight);
          else {
              //TODO: BSSRDF
            // Use a BSSRDF instead.
            //si->bsdf->Add(
            //    ARENA_ALLOC(arena, SpecularTransmission)(1.f, 1.f, e, mode));
            //si->bssrdf = ARENA_ALLOC(arena, DisneyBSSRDF)(c * diffuseWeight, sd,
            //                                              *si, e, this, mode);
          }
        }

        // Retro-reflection.
        bsdf->Add(new DisneyRetro(diffuseWeight * color, roughness),
                  diffuseWeight);

        // Sheen (if enabled)
        if (sheen > 0)
          bsdf->Add(new DisneySheen(diffuseWeight * sheen * Csheen),
                    diffuseWeight*sheen);
      }

      // Create the microfacet distribution for metallic and/or specular
      // transmission.
      float aspect = std::sqrt(1 - anisotropic * 0.9f);
      float ax = std::max(1e-3f, sqr(roughness) / aspect);
      float ay = std::max(1e-3f, sqr(roughness) * aspect);
      MicrofacetDistribution *distrib =
          new DisneyMicrofacetDistribution(ax, ay);

      // Specular is Trowbridge-Reitz with a modified Fresnel function.
      glm::vec3 Cspec0 =
          glm::mix(SchlickR0FromEta(eta) *
                       glm::mix(glm::vec3{1.0f}, Ctint, specularTint),
                   color, metallic);
      Fresnel *fresnel =
          new DisneyFresnel(Cspec0, metallic, eta);
      bsdf->Add(new MicrofacetReflection(glm::vec3{1.0f}, distrib, fresnel),1.0f);

      // Clearcoat
      if (clearcoat > 0) {
        bsdf->Add(new DisneyClearcoat(clearcoat,
                                      glm::mix(0.1f, 0.001f, clearcoatGloss)),
                  clearcoat);
      }

      // BTDF
      if (specTrans > 0) {
        // Walter et al's model, with the provided transmissive term scaled
        // by sqrt(color), so that after two refractions, we're back to the
        // provided color.
        glm::vec3 T = specTrans * glm::sqrt(color);
        if (thin) {
          // Scale roughness based on IOR (Burley 2015, Figure 15).
          float rscaled = (0.65f * eta - 0.35f) * roughness;
          float ax = std::max(1e-3f, sqr(rscaled) / aspect);
          float ay = std::max(1e-3f, sqr(rscaled) * aspect);
          MicrofacetDistribution *scaledDistrib = new TrowbridgeReitzDistribution(ax, ay);
          bsdf->Add(new MicrofacetTransmission(T, scaledDistrib, 1., eta),
                    specTrans);
        } else {
          MicrofacetDistribution *distribTransmission =
              new DisneyMicrofacetDistribution(ax, ay);
          bsdf->Add(new MicrofacetTransmission(T, distribTransmission, 1., eta), specTrans);
        }
      }
      if (thin) {
        // Lambertian, weighted by (1 - diffTrans)
        bsdf->Add(new LambertianTransmission(dt * color),dt);
      }
      break;
    }
  }
  return bsdf;
}

BSSRDF* Material::ComputeBSSRDF(const HitRecord &hit, const glm::vec3 direction, const Scene* scene) const {
    if(material_type != MATERIAL_TYPE_SUBSURFACE && material_type != MATERIAL_TYPE_KDSUBSURFACE)
      return nullptr;
    BSSRDF *newBSSRDF = new TabulatedBSSRDF(hit.position, direction, eta, BSSRDF_TABULATED, hit.geometry_normal, 
                                            this, sigma_a, sigma_s, table);
    return newBSSRDF;
}

}  // namespace sparks