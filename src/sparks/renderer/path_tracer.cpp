#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

namespace sparks {
PathTracer::PathTracer(const RendererSettings *render_settings,
                       const Scene *scene) {
  render_settings_ = render_settings;
  scene_ = scene;
  rd = std::mt19937(std::random_device()());
  uniform = std::uniform_real_distribution<float>(0, 1);
}

/*
glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample) const {
  glm::vec3 throughput{1.0f};
  glm::vec3 radiance{0.0f};
  HitRecord hit_record;
  const int max_bounce = render_settings_->num_bounces;
  std::mt19937 rd(sample ^ x ^ y);
  for (int i = 0; i < max_bounce; i++) {
    auto t = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit_record);
    if (t > 0.0f) {
      auto &material =
          scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
      if (material.material_type == MATERIAL_TYPE_EMISSION) {
        radiance += throughput * material.emission * material.emission_strength;
        break;
      } else {
        throughput *=
            material.albedo_color *
            glm::vec3{scene_->GetTextures()[material.albedo_texture_id].Sample(
                hit_record.tex_coord)};
        origin = hit_record.position;
        direction = scene_->GetEnvmapLightDirection();
        radiance += throughput * scene_->GetEnvmapMinorColor();
        throughput *=
            std::max(glm::dot(direction, hit_record.normal), 0.0f) * 2.0f;
        if (scene_->TraceRay(origin, direction, 1e-3f, 1e4f, nullptr) < 0.0f) {
          radiance += throughput * scene_->GetEnvmapMajorColor();
        }
        break;
      }
    } else {
      radiance += throughput * glm::vec3{scene_->SampleEnvmap(direction)};
      break;
    }
  }
  return radiance;
}
*/

// sample from light source, returns sample result(position on light source), pdf, normal vector and light
void PathTracer::SampleFromLight(glm::vec3 &res, glm::vec3 &norm, float &pdf) {
  const std::vector<Entity> &lightSources = scene_->GetEntities(); // all of them shall be meshes
  std::vector<float> prefixSum;
  std::vector<std::tuple<glm::vec3, glm::vec3, glm::vec3>> triangles; 
  std::vector<Material> materials;
  float totalArea = 0;
  // calculate total area
  for (int i = 0; i < lightSources.size(); i++) {
    const Entity &object = lightSources[i];
    if (object.GetMaterial().material_type != MATERIAL_TYPE_EMISSION)
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
  // uniformly sample
  float sample = uniform(rd) * totalArea;
  int totalTriangles = triangles.size();
  int triangleFound = std::lower_bound(prefixSum.begin(), prefixSum.end(), sample) - prefixSum.begin();
  float t = uniform(rd);
  float x = 1 - sqrt(1 - t); 
  float y = uniform(rd) * (1 - x);
  glm::vec3 v[3];
  std::cerr << triangleFound << " " << triangles.size() << std::endl;
  std::tie(v[0], v[1], v[2]) = triangles[triangleFound];
  res = x * v[0] + y * v[1] + (1 - x - y) * v[2];
  norm = glm::normalize(glm::cross(v[2] - v[0], v[1] - v[0]));
  pdf = 1 / totalArea;
}

// sample from cosine(among the hemisphere)
void PathTracer::SampleFromCosine(glm::vec3 &res, float &pdf, glm::vec3 localNorm) {
  float t1 = uniform(rd), t2 = uniform(rd);
  float theta = acos(sqrt(1 - t1)), phi = 2 * PI * t2;
  glm::vec3 localRes(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
  glm::mat4 local2World = glm::rotate(glm::mat4(1.0f), glm::dot(glm::vec3(0, 0, 1), localNorm), glm::cross(glm::vec3(0, 0, 1), localNorm));
  res = glm::vec3(local2World * glm::vec4(localRes, 1.0f));
  pdf = cos(theta) * sin(theta) / PI;
}

// path trace main algorithm
glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample,
                                bool moreBounces,
                                float weight) {           
  glm::vec3 emission{0.0f}, direct{0.0f}, incident{0.0f};
  HitRecord hit;
  float intersection = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit);
  if(intersection <= 0.0f)
    return glm::vec3{0.0f};
  glm::vec3 pos = hit.position;
  glm::vec3 norm = hit.geometry_normal;
  Material material = scene_->GetEntity(hit.hit_entity_id).GetMaterial();
  // emission
  emission = material.emission * material.emission_strength;
  if(!moreBounces)
    return emission * weight; // only direct illumination is effected by MIS
  // light source sampling
  glm::vec3 directSample;
  float directPdf;
  glm::vec3 directNorm;
  SampleFromLight(directSample, directNorm, directPdf);
  glm::vec3 directDir = glm::normalize(directSample - pos);
  direct = SampleRay(pos, directDir, x, y, sample, false, 1.0) * 
           material.albedo_color *
           glm::dot(directNorm, -directSample) *
           glm::dot(directNorm, norm) /
           glm::dot(directSample - pos, directSample - pos) /
           directPdf;
  float russianRoulette = uniform(rd);
  if(russianRoulette > continueProb)
    return emission + direct;
  // incident illumination
  glm::vec3 incidentSample;
  float incidentPdf;
  SampleFromCosine(incidentSample, incidentPdf, norm);
  float combineRatio = directPdf * directPdf / (directPdf * directPdf + incidentPdf * incidentPdf);
  incident = SampleRay(pos, -incidentSample, x, y, sample, true, 1 - combineRatio) * 
             material.albedo_color * 
             glm::dot(norm, incidentSample) /
             incidentPdf / 
             continueProb;
  return emission + combineRatio * direct + incident;
}
}  // namespace sparks
