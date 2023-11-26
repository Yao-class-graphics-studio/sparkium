#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <memory>

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

float PathTracer::getPdfByLight(glm::vec3 pos, glm::vec3 sample, float area) {
  HitRecord tmp;
  float t = scene_->TraceRay(pos, sample, 1e-3f, 1e4f, &tmp);
  if(t > 0.0f && scene_->GetEntity(tmp.hit_entity_id).GetMaterial().material_type == MATERIAL_TYPE_EMISSION)
    return 1 / area * t * t / glm::dot(tmp.geometry_normal, -sample);
  else  
    return 0.0f;
}

static float PowerHeuristic(float pdf1, float pdf2) {
  return pdf1 * pdf1 / (pdf1 * pdf1 + pdf2 * pdf2);
}

// path trace main algorithm
glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample,
                                int bounces) {           
  glm::vec3 emission{0.0f}, direct{0.0f}, env{0.0f}, incident{0.0f};
  HitRecord hit;
  float intersection = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit);
  auto FirstBounceClamp = [bounces](const glm::vec3 &L) -> glm::vec3 {
    return bounces == 0 ? glm::clamp(L, 0.0f, 1.0f) : L;
  };
  if(intersection <= 0.0f)
    return FirstBounceClamp(glm::vec3(scene_->SampleEnvmap(direction)));

  glm::vec3 pos = hit.position + 3e-5f * hit.geometry_normal;
  glm::vec3 norm = hit.geometry_normal;
  const Material &material = scene_->GetEntity(hit.hit_entity_id).GetMaterial();
  std::unique_ptr<BSDF> bsdf(material.ComputeBSDF(hit, scene_));

  //sample incident
  glm::vec3 incidentSample;
  float incidentPdf;
  BxDFType sampledType;
  glm::vec3 brdf;

  float actualContinueProb = bounces > 2 ? continueProb : 1;
  float russianRoulette = uniform(rd);
  if (russianRoulette > actualContinueProb || bounces >= render_settings_->num_bounces) {
    incidentPdf = 0;
    sampledType = BxDFType(0);
    brdf = glm::vec3{0.0f};
    incidentSample = glm::vec3{0.0f};
    incident = glm::vec3{0.0f};
  } else {
    brdf = bsdf->Sample_f(-direction, incidentSample,
                          glm::vec2{uniform(rd), uniform(rd)}, incidentPdf,
                          BxDFType(BSDF_ALL), sampledType);
    if (incidentPdf > 0.0f)
      incident = SampleRay(hit.position + 3e-5f * incidentSample,
                           incidentSample, x, y, sample, bounces + 1) *
                 brdf *
                 std::max(((sampledType & BSDF_TRANSMISSION) ? -1 : 1) *
                              glm::dot(norm, incidentSample),
                          0.0f) /
                 (actualContinueProb * incidentPdf);
        //TODO: BSSRDF
    else
      incident = glm::vec3{0.0f};
  }
  if (sampledType & BSDF_SPECULAR) {
    return FirstBounceClamp(incident);
  }

  // emission
  emission = material.material_type == MATERIAL_TYPE_EMISSION
                    ? material.emission * material.emission_strength
                    : glm::vec3{0.0f};

  // direct illumination (from environment map)
  // diffuse reflection only
  if (hit.front_face) {
    env += bsdf->f(-direction, -direction, BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)) * scene_->GetEnvmapMinorColor();
    glm::vec3 envDirection = scene_->GetEnvmapLightDirection();
    if (scene_->TraceRay(pos, envDirection, 1e-3f, 1e4f, nullptr) < 0.0f) {
      env += bsdf->f(-direction, envDirection,
                     BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)) *
             scene_->GetEnvmapMajorColor() *
             std::max(glm::dot(envDirection, norm), 0.0f) * 2.0f;
    }
  }

  // direct illumination (from emitting entities)
  glm::vec3 directSample, directDir = glm::vec3(0, 0, 1);
  float directPdf = 0, lightArea;
  glm::vec3 directNorm;
  SampleFromLight(directSample, directNorm, lightArea, hit.hit_entity_id);
  if(lightArea != 0) {
    directDir = glm::normalize(directSample - pos);
    HitRecord tmpHit;
    float check = scene_->TraceRay(pos, directDir, 1e-3f, 1e4f, &tmpHit);
    if(glm::length(pos + check * directDir - directSample) < 0.01f) {
      directPdf = 1 / lightArea * glm::dot(directSample - pos, directSample - pos) / std::fabs(glm::dot(directNorm, -directDir));
      Material lightMaterial = scene_->GetEntity(tmpHit.hit_entity_id).GetMaterial();
      direct = lightMaterial.emission * lightMaterial.emission_strength * bsdf->f(-direction, directDir, BxDFType(BSDF_ALL)) *
              std::max(glm::dot(norm, directDir), 0.0f) /
              directPdf;

      float directPdfX = bsdf->Pdf(-direction, directDir);

      float directRatio =
          directPdf == 0 ? 0 : PowerHeuristic(directPdf, directPdfX);
      direct *= directRatio;
    }
  } else {
    directPdf = 0;
  }
  //float russianRoulette = uniform(rd);
  //if (russianRoulette > continueProb) {
  //  return bounces == 0
  //             ? material.emission + glm::clamp(direct + env, 0.0f, 1.0f)
  //             : emission + direct + env;
  //}
  // incident illumination
  
  float incidentPdfX = getPdfByLight(pos, incidentSample, lightArea);

  float incidentRatio =
      incidentPdf == 0 ? 0 : PowerHeuristic(incidentPdf, incidentPdfX);
  //float directRatio = directPdf == 0 ? 0 : directPdf / (directPdf + directPdfX);
  //float incidentRatio = incidentPdf / (incidentPdf + incidentPdfX);
  //printf_s("%f %f %f %f\n", directPdf, directPdfX, incidentPdf, incidentPdfX);
  incident *= incidentRatio;

  glm::vec3 L = emission + direct + env + incident;
  if (bounces == 0)
    L = glm::clamp(L, 0.0f, 1.0f);
  return L;
}
}  // namespace sparks
