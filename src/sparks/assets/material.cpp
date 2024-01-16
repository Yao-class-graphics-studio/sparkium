#include "sparks/assets/material.h"

#include "grassland/grassland.h"
#include "sparks/assets/scene.h"
#include "sparks/assets/texture.h"
#include "sparks/util/util.h"
#include "sparks/assets/bxdf.h"
#include "sparks/assets/disney_material.h"
#include "glm/gtc/matrix_transform.hpp"
#include <fstream>

namespace sparks {

// utility

// Read density grid from file. Jingyi Lyu.
float* readDensityGrid(const std::string &path, int &nx, int &ny, int &nz) {
  std::ifstream in;
  in.open(path, std::ios::in);
  if(!in.is_open())
    return nullptr;
  in >> nx >> ny >> nz;
  std::cerr << nx << " " << ny << " " << nz << " " << path << std::endl;
  float *grid = new float[nx * ny * nz];
  for(int i = 0; i < nz; i++) 
    for(int j = 0; j < ny; j++)
      for(int k = 0; k < nx; k++)
        in >> grid[i * ny * nx + j * nx + k];
  in.close();
  return grid;
}
namespace {
std::unordered_map<std::string, MaterialType> material_name_map{
    {"lambertian", MATERIAL_TYPE_LAMBERTIAN},
    {"specular", MATERIAL_TYPE_SPECULAR},
    {"transmissive", MATERIAL_TYPE_TRANSMISSIVE},
    {"principled", MATERIAL_TYPE_PRINCIPLED},
    {"emission", MATERIAL_TYPE_EMISSION},
    {"subsurface", MATERIAL_TYPE_SUBSURFACE},
    {"kdsubsurface", MATERIAL_TYPE_KDSUBSURFACE},
    {"medium", MATERIAL_TYPE_MEDIUM},
    {"grid_medium", MATERIAL_TYPE_GRID_MEDIUM}
};

// Map for BSSRDF Type. Jingyi Lyu.
std::unordered_map<std::string, PresetSSType> ss_type_name_map{
    {"apple", SS_APPLE},
    {"chicken_1", SS_CHICKEN_1},
    {"chicken_2", SS_CHICKEN_2},
    {"cream", SS_CREAM},
    {"ketchup", SS_KETCHUP},
    {"marble", SS_MARBLE},
    {"potato", SS_POTATO},
    {"skimmilk", SS_SKIMMILK},
    {"skin_1", SS_SKIN_1},
    {"skin_2", SS_SKIN_2},
    {"spectralon", SS_SPECTRALON},
    {"wholemilk", SS_WHOLEMILK},
    {"lowfat_milk", SS_LOWFAT_MILK},
    {"reduced_milk", SS_REDUCED_MILK},
    {"regular_milk", SS_REGULAR_MILK},
    {"espresso", SS_ESPRESSO},
    {"mint_mocha_coffee", SS_MINT_MOCHA_COFFEE},
    {"lowfat_soy_milk", SS_LOWFAT_SOY_MILK},
    {"regular_soy_milk", SS_REGULAR_SOY_MILK},
    {"lowfat_chocolate_milk", SS_LOWFAT_CHOCOLATE_MILK},
    {"regular_chocolate_milk", SS_REGULAR_CHOCOLATE_MILK},
    {"coke", SS_COKE},
    {"pepsi", SS_PEPSI},
    {"sprite", SS_SPRITE},
    {"gatorade", SS_GATORADE},
    {"chardonnay", SS_CHARDONNAY},
    {"white_zinfandel", SS_WHITE_ZINFANDEL},
    {"merlot", SS_MERLOT},
    {"budweiser_beer", SS_BUDWEISER_BEER},
    {"coors_light_beer", SS_COORS_LIGHT_BEER},
    {"clorox", SS_CLOROX},
    {"apple_juice", SS_APPLE_JUICE},
    {"cranberry_juice", SS_CRANBERRY_JUICE},
    {"grape_juice", SS_GRAPE_JUICE},
    {"ruby_grape_juice", SS_RUBY_GRAPE_JUICE},
    {"white_grapefuite_juice", SS_WHITE_GRAPEfRUITE_JUICE},
    {"shampoo", SS_SHAMPOO},
    {"strawberry_shampoo", SS_STRAWBERRY_SHAMPOO},
    {"head_and_shoulders_shampoo", SS_HEAD_AND_SHOULDERS_SHAMPOO},
    {"lemon_tee_powder", SS_LEMON_TEE_POWDER},
    {"orange_powder", SS_ORANGE_POWDER},
    {"pink_lemonade_powder", SS_PINK_LEMONADE_POWDER},
    {"cappuccino_powder", SS_CAPPUCCINO_POWDER},
    {"salt_powder", SS_SALT_POWDER},
    {"sugar_powder", SS_SUGAR_POWDER},
    {"suisse_mocha_powder", SS_SUISSE_MOCHA_POWDER},
    {"pacific_ocean_surface_water", SS_PACIFIC_OCEAN_SURFACE_WATER},
    {"none", SS_NONE}
};

}

// Read material from XML. Jingyi Lyu + Shengquan Du. 
Material::Material(Scene *scene, const tinyxml2::XMLElement *material_element)
    : table(new BSSRDFTable(100, 64))  {
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
    std::string filename = PathToFilename(path);
    int tmp_id = scene->GetTextureId(filename);
    if (tmp_id != -1)
    {
      use_normal_texture = true;
      normal_texture_id = tmp_id;
    } else {
      Texture normal_texture(1, 1);
      if (Texture::Load(path, normal_texture)) {
        use_normal_texture = true;
        normal_texture_id =
            scene->AddTexture(normal_texture, PathToFilename(path));
      }
    }
  }

  child_element = material_element->FirstChildElement("alpha_texture");
  if (child_element) {
    std::string path = child_element->FindAttribute("value")->Value();
    std::string filename = PathToFilename(path);
    int tmp_id = scene->GetTextureId(filename);
    if (tmp_id != -1) {
      use_alpha_texture = true;
      alpha_texture_id = tmp_id;
    } else {
      Texture alpha_texture(1, 1);
      if (Texture::Load(path, alpha_texture)) {
        use_alpha_texture = true;
        alpha_texture_id =
            scene->AddTexture(alpha_texture, PathToFilename(path));
      }
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

  child_element = material_element->FirstChildElement("cone");
  if (child_element) {
    cone = std::stof(child_element->FindAttribute("value")->Value());
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

  child_element = material_element->FirstChildElement("sigma_a");
  if (child_element) {
    sigma_a = StringToVec3(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("sigma_s");
  if (child_element) {
    sigma_s = StringToVec3(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("ss_type");
  if (child_element) {
    std::string type = child_element->FindAttribute("value")->Value();
    ss_type = ss_type_name_map[type];
    GetPresetSS(ss_type, sigma_a, sigma_s);
  }

  child_element = material_element->FirstChildElement("false_surface");
  if (child_element) {
    false_surface = StringToBool(child_element->FindAttribute("value")->Value());
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
    ComputeBeamDiffusionBSSRDF(g, eta, *table);
  }

  if (material_type == MATERIAL_TYPE_MEDIUM) {
    medium = new HomogeneousMedium(sigma_a, sigma_s, volumetric_emission, g);
  }

  if (material_type == MATERIAL_TYPE_GRID_MEDIUM) {
    std::string path;
    child_element = material_element->FirstChildElement("density_grid_path");
    if (child_element) {
      path = child_element->FindAttribute("value")->Value();
    }

    glm::vec3 pos000, pos001, pos010, pos100;
    child_element = material_element->FirstChildElement("pos000");
    if (child_element) {
      pos000 = StringToVec3(child_element->FindAttribute("value")->Value());
    }
    child_element = material_element->FirstChildElement("pos001");
    if (child_element) {
      pos001 = StringToVec3(child_element->FindAttribute("value")->Value());
    }
    child_element = material_element->FirstChildElement("pos010");
    if (child_element) {
      pos010 = StringToVec3(child_element->FindAttribute("value")->Value());
    }
    child_element = material_element->FirstChildElement("pos100");
    if (child_element) {
      pos100 = StringToVec3(child_element->FindAttribute("value")->Value());
    }

    int nx, ny, nz;
    float *grid = readDensityGrid(path, nx, ny, nz);
    assert(grid != nullptr);
    glm::mat3 scale{pos100 - pos000, pos010 - pos000, pos001 - pos000};
    glm::mat4 transform{
        glm::vec4{scale[0], 0.0f},
        glm::vec4{scale[1], 0.0f},
        glm::vec4{scale[2], 0.0f},
        glm::vec4(pos000[0], pos000[1], pos000[2], 1.0f)
    };
    medium = new GridDensityMedium(sigma_a, sigma_s, volumetric_emission, g, 
                                   nx, ny, nz, glm::inverse(transform), glm::inverse(scale), grid);
  }
}

Material::Material(const glm::vec3 &albedo) {
  albedo_color = albedo;
}

// Shengquan Du
glm::vec3 Material::GetAlbedoColor(const HitRecord &hit,
                                   const Scene *scene) const {
  return albedo_color * glm::vec3(scene->GetTexture(albedo_texture_id).Sample(hit.tex_coord));
}

// Shengquan Du
float Material::GetAlpha(const HitRecord &hit, const Scene *scene) const {
  if (!use_alpha_texture)
    return alpha;
  return scene->GetTexture(alpha_texture_id).Sample(hit.tex_coord)[3];
}

// Return actual hit (considering normal texture etc.). Shengquan Du.
HitRecord Material::GetShaderHit(const HitRecord &hit, const Scene *scene) const{
  HitRecord textureHit = hit;
  if (!textureHit.front_face) {
    textureHit.front_face = true;
    textureHit.geometry_normal *= -1;
    textureHit.normal *= -1;
    textureHit.tangent *= -1;
  }
  if (use_normal_texture) {
      // seems not very correct, to be checked.....
    glm::mat3 local_to_world;
    glm::vec3 norm = textureHit.normal;
    if (std::fabs(norm.z) > 1.0f - 1e-6) {
      local_to_world = glm::mat3(1.0f);
    } else {
      glm::vec3 y =
          glm::normalize(glm::cross(norm, glm::vec3(0.0f, 0.0f, 1.0f)));
      glm::vec3 x = glm::normalize(glm::cross(y, norm));
      local_to_world = glm::mat3(x, y, norm);
    }
    glm::vec3 textureValue =
        scene->GetTexture(normal_texture_id).Sample(hit.tex_coord);
    textureValue = glm::normalize(textureValue - glm::vec3{0.5f});
    textureHit.normal = glm::normalize(local_to_world * textureValue);
  }
  return textureHit;
}

// Shengquan Du
glm::vec3 Material::GetShaderNormal(const HitRecord &hit,
                                    const Scene *scene) const {
  if (!use_normal_texture)
    return hit.front_face ? hit.normal : -hit.normal;
  else
    return GetShaderHit(hit, scene).normal;
}

// Shengquan Du.
BSDF* Material::ComputeBSDF(const HitRecord &hit,
                         const Scene* scene) const {
  glm::vec3 color = GetAlbedoColor(hit, scene);
  float real_alpha = GetAlpha(hit, scene);
  //float real_alpha = alpha; // It seems that there are still some bugs on alpha channel, so do not use it
  BSDF* bsdf = new BSDF(GetShaderHit(hit,scene));
  switch (material_type) { 
      case MATERIAL_TYPE_SUBSURFACE: {
          bsdf->Add(new SpecularTransmission(glm::vec3{1.0f}, 1.0f, 1.0f), 1.0f);
          break;
      }
      case MATERIAL_TYPE_KDSUBSURFACE: {
          bsdf->Add(new SpecularTransmission(glm::vec3{1.0f}, 1.0f, 1.0f), 1.0f);
          break;
      }
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
        bsdf->Add(new SpecularTransmission(glm::vec3{ 1.0f }, 1.0f, eta), 1);
        //bsdf->Add(new SpecularReflection(glm::vec3{ 1.0f },
        //    new FresnelDielectric(1.0f, eta)),
        //    1.0f);
        break;
    }
    case MATERIAL_TYPE_MEDIUM: {
        bsdf->Add(new SpecularTransmission(color, 1.0f, eta), 1);
        break;
    }
    case MATERIAL_TYPE_GRID_MEDIUM: {
        bsdf->Add(new SpecularTransmission(color, 1.0f, eta), 1);
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
            bsdf->Add(new SpecularTransmission(glm::vec3{1.0f}, 1.0f, eta),diffuseWeight);
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
      bsdf->Add(new MicrofacetReflection(color, distrib, fresnel),1.0f);

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
        if (eta > 1.0f) {
          if (thin) {
            // Scale roughness based on IOR (Burley 2015, Figure 15).
            float rscaled = (0.65f * eta - 0.35f) * roughness;
            float ax = std::max(1e-3f, sqr(rscaled) / aspect);
            float ay = std::max(1e-3f, sqr(rscaled) * aspect);
            MicrofacetDistribution *scaledDistrib =
                new TrowbridgeReitzDistribution(ax, ay);
            bsdf->Add(new MicrofacetTransmission(T, scaledDistrib, 1., eta),
                      specTrans);
          } else {
            MicrofacetDistribution *distribTransmission =
                new DisneyMicrofacetDistribution(ax, ay);
            bsdf->Add(
                new MicrofacetTransmission(T, distribTransmission, 1., eta),
                specTrans);
          }
        } else {
          bsdf->Add(new SpecularTransmission(T, 1.0f, 1.0f), specTrans);
        }
      }
      if (thin) {
        // Lambertian, weighted by (1 - diffTrans)
        bsdf->Add(new LambertianTransmission(dt * color),dt);
      }
      break;
    }
  }
  if (real_alpha != 1.0f) {
    float sum_weight = 0;
    for (auto &[bxdf, weight] : bsdf->bxdfs_) {
      bxdf = new ScaledBxDF(bxdf, glm::vec3{real_alpha});
      sum_weight += weight;
      weight *= real_alpha;
    }
    bsdf->Add(new AlphaTransmission(glm::vec3{1.0f-real_alpha}), sum_weight * (1 - real_alpha));
  }
  return bsdf;
}

// Jingyi Lyu.
BSSRDF* Material::ComputeBSSRDF(const HitRecord &hit, const glm::vec3 direction, const Scene* scene) {
  if (material_type != MATERIAL_TYPE_SUBSURFACE &&
      material_type != MATERIAL_TYPE_KDSUBSURFACE &&
      (material_type != MATERIAL_TYPE_PRINCIPLED ||
       scatterDistance == glm::vec3{0.0f}))
      return nullptr;
  HitRecord textureHit = GetShaderHit(hit,scene);
  if (material_type == MATERIAL_TYPE_KDSUBSURFACE || material_type ==MATERIAL_TYPE_SUBSURFACE) {
    if (material_type == MATERIAL_TYPE_KDSUBSURFACE)
      SubsurfaceFromDiffuse(*table, GetAlbedoColor(hit,scene), mfp, sigma_a, sigma_s);
    BSSRDF *newBSSRDF = new TabulatedBSSRDF(
        hit.position, -direction, eta, BSSRDF_TABULATED, textureHit.normal,
        this, sigma_a, sigma_s, *table);
    return newBSSRDF;
  } else if (material_type == MATERIAL_TYPE_PRINCIPLED &&
      scatterDistance != glm::vec3{0.0f} && !thin) {
    float diffuseWeight = (1 - metallic) * (1 - specTrans);
    if (diffuseWeight == 0.0f)
      return nullptr;
    glm::vec3 color = GetAlbedoColor(hit,scene);
    return new DisneyBSSRDF(color * diffuseWeight, scatterDistance,
                            textureHit.position, -direction, textureHit.normal,
                            eta, this);
  }
  return nullptr;
}

}  // namespace sparks