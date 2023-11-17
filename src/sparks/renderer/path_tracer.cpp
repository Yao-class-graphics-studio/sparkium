#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

#include <cassert>

namespace sparks {
PathTracer::PathTracer(const RendererSettings *render_settings,
                       const Scene *scene) {
  render_settings_ = render_settings;
  scene_ = scene;
  rd = std::mt19937(std::random_device()());
  uniform = std::uniform_real_distribution<float>(0, 1);
}


// glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
//                                 glm::vec3 direction,
//                                 int x,
//                                 int y,
//                                 int sample,
//                                 bool moreBounces,
//                                 float weight) {
//   glm::vec3 throughput{1.0f};
//   glm::vec3 radiance{0.0f};
//   HitRecord hit_record;
//   const int max_bounce = render_settings_->num_bounces;
//   std::mt19937 rd(sample ^ x ^ y);
//   for (int i = 0; i < max_bounce; i++) {
//     auto t = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit_record);
//     if (t > 0.0f) {
//       auto &material =
//           scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
//       if (material.material_type == MATERIAL_TYPE_EMISSION) {
//         radiance += throughput * material.emission * material.emission_strength;
//         break;
//       } else {
//         throughput *=
//             material.albedo_color *
//             glm::vec3{scene_->GetTextures()[material.albedo_texture_id].Sample(
//                 hit_record.tex_coord)};
//         origin = hit_record.position;
//         direction = scene_->GetEnvmapLightDirection();
//         radiance += throughput * scene_->GetEnvmapMinorColor();
//         throughput *=
//             std::max(glm::dot(direction, hit_record.normal), 0.0f) * 2.0f;
//         if (scene_->TraceRay(origin, direction, 1e-3f, 1e4f, nullptr) < 0.0f) {
//           radiance += throughput * scene_->GetEnvmapMajorColor();
//         }
//         break;
//       }
//     } else {
//       radiance += throughput * glm::vec3{scene_->SampleEnvmap(direction)};
//       break;
//     }
//   }
//   return radiance;
// }


// sample from light source, returns sample result(position on light source), pdf, normal vector and light
void PathTracer::SampleFromLight(glm::vec3 &res, glm::vec3 &norm, float &area, int except) {
  const std::vector<Entity> &lightSources = scene_->GetEntities(); // all of them shall be meshes
  std::vector<float> prefixSum;
  std::vector<std::tuple<glm::vec3, glm::vec3, glm::vec3>> triangles; 
  std::vector<Material> materials;
  float totalArea = 0;
  // calculate total area
  for (int i = 0; i < lightSources.size(); i++) {
    const Entity &object = lightSources[i];
    if (object.GetMaterial().material_type != MATERIAL_TYPE_EMISSION || i == except)
      continue;
    const std::vector<Vertex> &vertices = object.GetModel()->GetVertices();
    const std::vector<uint32_t> &indices = object.GetModel()->GetIndices();
    glm::mat4 transform = object.GetTransformMatrix();
    for(int i = 0; i < indices.size(); i += 3) {
      int j = i + 1, k = i + 2;
      glm::vec4 v0(vertices[i].position, 1.0f);
      glm::vec4 v1(vertices[j].position, 1.0f);
      glm::vec4 v2(vertices[k].position, 1.0f);
      glm::vec3 u0 = glm::vec3(transform * v0);
      glm::vec3 u1 = glm::vec3(transform * v1);
      glm::vec3 u2 = glm::vec3(transform * v2);
      float nowArea = glm::length(glm::cross(u2 - u0, u1 - u0)) / 2;
      totalArea += nowArea;
      prefixSum.emplace_back(totalArea);
      triangles.emplace_back(std::make_tuple(u0, u1, u2));
      materials.emplace_back(object.GetMaterial());
    }
  }
  if(totalArea == 0) {
    area = 0;
    return;
  }
  // uniformly sample
  float sample = uniform(rd) * totalArea;
  int totalTriangles = triangles.size();
  int triangleFound = std::lower_bound(prefixSum.begin(), prefixSum.end(), sample) - prefixSum.begin();
  float t = uniform(rd);
  float x = 1 - sqrt(1 - t); 
  float y = uniform(rd) * (1 - x);
  glm::vec3 v[3];
  std::tie(v[0], v[1], v[2]) = triangles[triangleFound];
  res = x * v[0] + y * v[1] + (1 - x - y) * v[2];
  norm = glm::normalize(glm::cross(v[2] - v[0], v[1] - v[0]));
  area = totalArea;
}

// sample from cosine(among the hemisphere)
void PathTracer::SampleFromCosine(glm::vec3 &res, float &pdf, glm::vec3 localNorm) {
  float t1 = uniform(rd), t2 = uniform(rd);
  float theta = acos(sqrt(1 - t1)), phi = 2 * PI * t2;
  glm::vec3 localRes(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
  float dotProduct = glm::dot(localNorm, glm::vec3(0, 0, 1));
  if(std::fabs(dotProduct) == 1.0f)
    res = dotProduct > 0 ? localRes : -1.0f * localRes;
  else {
    glm::mat4 local2World = glm::rotate(glm::mat4(1.0f), acos(glm::dot(glm::vec3(0, 0, 1), localNorm)), glm::normalize(glm::cross(glm::vec3(0, 0, 1), localNorm)));
    res = glm::vec3(local2World * glm::vec4(localRes, 1.0f));
  }
  res = glm::normalize(res);
  pdf = cos(theta) * sin(theta) / PI;
}

void PathTracer::Sample(glm::vec3 &res, float &pdf, glm::vec3 in, glm::vec3 localNorm, Material &material) {
  switch(material.material_type) {
    case MATERIAL_TYPE_LAMBERTIAN: 
      SampleFromCosine(res, pdf, localNorm);
      break;
    case MATERIAL_TYPE_EMISSION:
      SampleFromCosine(res, pdf, localNorm);
      break;
    case MATERIAL_TYPE_SPECULAR:
      res = -2.0f * dot(localNorm, in) * localNorm + in, pdf = 1000.0f;
      res = glm::normalize(res);
      break;
  }
}

float PathTracer::getPdfByMaterial(glm::vec3 sample, glm::vec3 in, glm::vec3 norm, Material &material) {
  float theta = 0;
  switch(material.material_type) {
    case MATERIAL_TYPE_LAMBERTIAN: 
      theta = std::max(glm::dot(sample, norm), 0.0f);
      return cos(theta) * sin(theta) / PI;
    case MATERIAL_TYPE_EMISSION:
      theta = std::max(glm::dot(sample, norm), 0.0f);
      return cos(theta) * sin(theta) / PI;
    case MATERIAL_TYPE_SPECULAR:
      return sample == glm::normalize(-2.0f * dot(norm, in) * norm + in) ? 1000.0f : 0;
    default:
      return 0.0f;
  }
}

float PathTracer::getPdfByLight(glm::vec3 pos, glm::vec3 sample, float area) {
  HitRecord tmp;
  float t = scene_->TraceRay(pos, sample, 1e-3f, 1e4f, &tmp);
  if(t > 0.0f && scene_->GetEntity(tmp.hit_entity_id).GetMaterial().material_type == MATERIAL_TYPE_EMISSION)
    return 1 / area * t * t / glm::dot(tmp.geometry_normal, -sample);
  else  
    return 0.0f;
}

glm::vec3 PathTracer::getBRDF(Material &material, HitRecord &hit, glm::vec3 norm, glm::vec3 in, glm::vec3 out) {
  glm::vec3 color = material.albedo_color * glm::vec3(scene_->GetTexture(material.albedo_texture_id).Sample(hit.tex_coord));
  // std::cerr << color.x << std::endl;
  switch(material.material_type) {
    case MATERIAL_TYPE_LAMBERTIAN: 
      return color / PI;
    case MATERIAL_TYPE_EMISSION:
      return color / PI;
    case MATERIAL_TYPE_SPECULAR:
      return out == glm::normalize(-2.0f * dot(norm, in) * norm + in) ? glm::vec3{1.0f} * 1000.0f / glm::dot(out, norm) : glm::vec3{0.0f};
    default:
      return glm::vec3{0.0f}; 
  }
}

// path trace main algorithm
glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample,
                                bool moreBounces,
                                float weight,
                                int bounces) {           
  glm::vec3 emission{0.0f}, direct{0.0f}, env{0.0f}, incident{0.0f};
  HitRecord hit;
  float intersection = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit);
  if(intersection <= 0.0f)
    return bounces == 0 ? glm::vec3(scene_->SampleEnvmap(direction)) : glm::vec3{0.0f};
  glm::vec3 pos = hit.position + 3e-5f * hit.geometry_normal;
  glm::vec3 norm = hit.geometry_normal;
  Material material = scene_->GetEntity(hit.hit_entity_id).GetMaterial();
  // emission
  emission = material.emission * material.emission_strength * weight;
  if(!moreBounces)
    return emission; // only direct illumination is effected by MIS
  // direct illumination (from environment map)
  glm::vec3 throughput = material.albedo_color * glm::vec3(scene_->GetTexture(material.albedo_texture_id).Sample(hit.tex_coord));
  glm::vec3 envDirection = scene_->GetEnvmapLightDirection();
  env += throughput * scene_->GetEnvmapMinorColor();
  if(scene_->TraceRay(pos, envDirection, 1e-3f, 1e4f, nullptr) < 0.0f)
    env += throughput * scene_->GetEnvmapMajorColor() * std::max(glm::dot(envDirection, norm), 0.0f) * 2.0f;
  // direct illumination (from emitting entities)
  glm::vec3 directSample, directDir = glm::vec3(0, 0, 1);
  float directPdf, lightArea;
  glm::vec3 directNorm;
  SampleFromLight(directSample, directNorm, lightArea, hit.hit_entity_id);
  if(lightArea != 0) {
    directDir = glm::normalize(directSample - pos);
    HitRecord tmpHit;
    float check = scene_->TraceRay(pos, directDir, 1e-3f, 1e4f, &tmpHit);
    if(glm::length(pos + check * directDir - directSample) < 0.01f) {
      directPdf = 1 / lightArea * glm::dot(directSample - pos, directSample - pos) / std::fabs(glm::dot(directNorm, -directDir));
      Material lightMaterial = scene_->GetEntity(tmpHit.hit_entity_id).GetMaterial();
      direct = lightMaterial.emission * lightMaterial.emission_strength *
              getBRDF(material, hit, norm, direction, directDir) *
              std::max(glm::dot(norm, directDir), 0.0f) /
              directPdf;
    }
  } else {
    directPdf = 0;
  }
  float russianRoulette = uniform(rd);
  if(russianRoulette > continueProb)
    return bounces == 0 ? material.emission + glm::clamp(direct + env, 0.0f, 1.0f) : emission + direct + env;
  // incident illumination
  glm::vec3 incidentSample;
  float incidentPdf;
  Sample(incidentSample, incidentPdf, direction, norm, material);
  float incidentPdfX = getPdfByLight(pos, incidentSample, lightArea);
  float directPdfX = getPdfByMaterial(directDir, direction, norm, material);
  float directRatio = directPdf == 0 ? 0 : directPdf * directPdf / (directPdf * directPdf + directPdfX * directPdfX);
  float incidentRatio = incidentPdf * incidentPdf / (incidentPdf * incidentPdf + incidentPdfX * incidentPdfX);
  incident = SampleRay(pos, incidentSample, x, y, sample, true, incidentRatio, bounces + 1) *
             getBRDF(material, hit, norm, direction, incidentSample) * 
             std::max(glm::dot(norm, incidentSample), 0.0f) /
             incidentPdf;
  if(bounces == 0)
    return glm::clamp(material.emission + (directRatio * direct + env + incident) / continueProb, 0.0f, 1.0f);
  else
    return emission + (directRatio * direct + env + incident) / continueProb;
}
}  // namespace sparks
